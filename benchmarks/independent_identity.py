#!/usr/bin/env python3
"""Separately coded TxID-v1 identity oracle; standard library only.

This module reads the published specification, not TxID's implementation objects.
It deliberately has no imports from ``txid``. Its strict GTF subset extracts
transcript identity, not annotation classification, abundances or provenance
attributes. See docs/independent-identity.md for the validation boundary.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile
from dataclasses import dataclass, field
from typing import Iterable, Mapping, Sequence


REFERENCE_VERSION = "independent-exact-v1.0"
_FINGERPRINT = re.compile(r"sha256:[0-9a-f]{64}\Z")
_ATTRIBUTE = re.compile(r"[A-Za-z_][A-Za-z0-9_.:-]*")
_POSITIVE_INTEGER = re.compile(r"[0-9]+\Z")


class IndependentIdentityError(ValueError):
    """An input does not satisfy this verifier's documented identity contract."""


def _token(value: object, name: str) -> str:
    if not isinstance(value, str) or not value or any(c.isspace() or ord(c) < 32 or ord(c) == 127 for c in value):
        raise IndependentIdentityError(f"{name} must be a nonempty whitespace-free string")
    return value


def _exon_intervals(exons: Iterable[Sequence[int]]) -> tuple[tuple[int, int], ...]:
    intervals: list[tuple[int, int]] = []
    for pair in exons:
        if len(pair) != 2 or any(type(n) is not int for n in pair):
            raise IndependentIdentityError("each exon must contain exactly two integer coordinates")
        start, end = pair
        if start < 1 or end < start:
            raise IndependentIdentityError(f"invalid 1-based closed exon [{start}, {end}]")
        intervals.append((start, end))
    if not intervals:
        raise IndependentIdentityError("a transcript must contain at least one exon")
    intervals.sort()
    for left, right in zip(intervals, intervals[1:]):
        if right[0] <= left[1]:
            raise IndependentIdentityError(f"duplicate or overlapping exons: {left}, {right}")
        if right[0] == left[1] + 1:
            raise IndependentIdentityError(f"adjacent exons imply an empty intron: {left}, {right}")
    return tuple(intervals)


def compute_identity(
    *,
    assembly_fingerprint: str,
    contig: str,
    strand: str,
    exons: Iterable[Sequence[int]],
) -> dict[str, dict[str, str] | None]:
    """Compute complete canonical JSON, SHA-256 and 96-bit public SC/TF/SE IDs.

    ``contig`` must already be a primary name in the supplied assembly context.
    The supplied fingerprint is syntax-checked, not inferred from the annotation.
    Exons are 1-based closed; record order has no meaning.
    """
    if not isinstance(assembly_fingerprint, str) or not _FINGERPRINT.fullmatch(assembly_fingerprint):
        raise IndependentIdentityError("assembly_fingerprint must be sha256: followed by 64 lowercase hex digits")
    _token(contig, "contig")
    if strand not in ("+", "-"):
        raise IndependentIdentityError("transcript strand must be + or -")
    intervals = _exon_intervals(exons)
    low, high = intervals[0][0], intervals[-1][1]
    if len(intervals) == 1:
        families = {"form": ("SE1", {"start": low, "end": high})}
    else:
        gaps = [[a[1] + 1, b[0] - 1] for a, b in zip(intervals, intervals[1:])]
        if strand == "-":
            gaps.reverse()
        families = {
            "splice_chain": ("SC1", {"introns": gaps}),
            "form": ("TF1", {
                "introns": gaps,
                "tss": low if strand == "+" else high,
                "tes": high if strand == "+" else low,
            }),
        }
    result: dict[str, dict[str, str] | None] = {"splice_chain": None}
    for role, (family, coordinates) in families.items():
        obj = dict(coordinates, algorithm=family, assembly=assembly_fingerprint,
                   contig=contig, strand=strand)
        # Spell out ordering independently; default JSON separators contain spaces.
        ordered = {name: obj[name] for name in sorted(obj)}
        payload = json.dumps(ordered, ensure_ascii=False, separators=(",", ":"))
        digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()
        result[role] = {
            "canonical_json": payload,
            "full_digest": digest,
            "public_id": "txid:" + family + "." + digest[:24],
        }
    return result


@dataclass(frozen=True)
class IndependentTranscript:
    transcript_id: str
    gene_id: str | None
    contig: str
    strand: str
    exons: tuple[tuple[int, int], ...]

    def identity(self, assembly_fingerprint: str) -> dict[str, dict[str, str] | None]:
        return compute_identity(assembly_fingerprint=assembly_fingerprint,
                                contig=self.contig, strand=self.strand, exons=self.exons)


def _attributes(text: str) -> list[tuple[str, str]]:
    """Read quoted GTF values without splitting semicolons inside a string."""
    if text == ".":
        return []
    result = []
    index = 0
    while index < len(text):
        while index < len(text) and text[index].isspace():
            index += 1
        if index == len(text):
            break
        match = _ATTRIBUTE.match(text, index)
        if match is None:
            raise IndependentIdentityError("malformed GTF attribute name")
        key = match.group()
        index = match.end()
        if index == len(text) or not text[index].isspace():
            raise IndependentIdentityError(f"attribute {key!r} must be followed by whitespace and a value")
        while index < len(text) and text[index].isspace():
            index += 1
        if index == len(text):
            raise IndependentIdentityError(f"attribute {key!r} has no value")
        if text[index] == '"':
            index += 1
            value: list[str] = []
            while index < len(text) and text[index] != '"':
                character = text[index]
                index += 1
                if character == "\\":
                    if index == len(text) or text[index] not in {'"', "\\", "n", "r", "t"}:
                        raise IndependentIdentityError(f"unsupported or truncated escape in attribute {key!r}")
                    character = {'n': '\n', 'r': '\r', 't': '\t'}.get(text[index], text[index])
                    index += 1
                value.append(character)
            if index == len(text):
                raise IndependentIdentityError(f"unterminated quoted attribute {key!r}")
            index += 1
            parsed = "".join(value)
        else:
            start = index
            while index < len(text) and not text[index].isspace() and text[index] != ";":
                if text[index] in {'"', "\\", "="}:
                    raise IndependentIdentityError(f"malformed unquoted value for attribute {key!r}")
                index += 1
            parsed = text[start:index]
            if not parsed:
                raise IndependentIdentityError(f"attribute {key!r} has no value")
        while index < len(text) and text[index].isspace():
            index += 1
        if index < len(text):
            if text[index] != ";":
                raise IndependentIdentityError(f"missing semicolon after attribute {key!r}")
            index += 1
        result.append((key, parsed))
    for key in ("gene_id", "transcript_id"):
        values = [value for name, value in result if name == key]
        if len(values) > 1:
            raise IndependentIdentityError(f"duplicate {key} attribute")
        if values and (not values[0] or any(ord(c) < 32 or ord(c) == 127 for c in values[0])):
            raise IndependentIdentityError(f"empty or control-containing {key}")
    return result


@dataclass
class _TranscriptRows:
    contig: str
    strand: str
    first_line: int
    genes: set[str] = field(default_factory=set)
    exons: list[tuple[int, int]] = field(default_factory=list)
    bounds: tuple[int, int] | None = None


def _context(
    aliases: Mapping[str, str] | None,
    lengths: Mapping[str, int] | None,
) -> tuple[dict[str, str], dict[str, int] | None]:
    names = dict(aliases or {})
    sizes = None if lengths is None else dict(lengths)
    if sizes is not None:
        if not sizes:
            raise IndependentIdentityError("contig_lengths must not be empty")
        for name, length in sizes.items():
            _token(name, "primary contig")
            if type(length) is not int or length < 1:
                raise IndependentIdentityError(f"invalid reference contig length for {name!r}")
    if names and sizes is None:
        raise IndependentIdentityError("explicit contig aliases require contig_lengths to validate primary targets")
    for alias, target in names.items():
        _token(alias, "contig alias")
        _token(target, "alias target")
        if sizes is None or target not in sizes:
            raise IndependentIdentityError(f"alias {alias!r} targets unknown primary contig {target!r}")
        if alias in sizes and alias != target:
            raise IndependentIdentityError(f"alias {alias!r} conflicts with a primary contig")
    return names, sizes


def read_gtf_models(
    path: str | Path,
    *,
    contig_aliases: Mapping[str, str] | None = None,
    contig_lengths: Mapping[str, int] | None = None,
) -> dict[str, IndependentTranscript]:
    """Parse a plain/gzip GTF independently and return transcript-ID keyed models.

    Transcript meta-features are optional, but any present bounds must exactly
    span their exons. Duplicate/overlapping/adjacent exons and conflicting IDs,
    contigs, strands, genes or transcript rows are errors. A supplied length map
    restricts the reference domain and enables upper-bound validation. Explicit
    aliases require such a map; no contig-name heuristic is used.
    """
    aliases, lengths = _context(contig_aliases, contig_lengths)
    source = Path(path)
    with source.open("rb") as handle:
        compressed = handle.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    records: dict[str, _TranscriptRows] = {}
    with opener(source, "rt", encoding="utf-8", newline=None) as handle:
        for number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            try:
                columns = raw.rstrip("\n").split("\t")
                if len(columns) != 9:
                    raise IndependentIdentityError("GTF records must have exactly nine tab-separated columns")
                contig, origin, feature, start_text, end_text, score, strand, phase, attrs = columns
                _token(contig, "contig")
                _token(origin, "source")
                _token(feature, "feature")
                if not _POSITIVE_INTEGER.fullmatch(start_text) or not _POSITIVE_INTEGER.fullmatch(end_text):
                    raise IndependentIdentityError("GTF coordinates must be positive decimal integers")
                start, end = int(start_text), int(end_text)
                if start < 1 or end < start:
                    raise IndependentIdentityError("invalid 1-based closed GTF coordinates")
                if score != "." and not math.isfinite(float(score)):
                    raise IndependentIdentityError("GTF score must be finite or '.'")
                if phase not in (".", "0", "1", "2"):
                    raise IndependentIdentityError("GTF phase must be '.', 0, 1 or 2")
                if strand not in ("+", "-", "."):
                    raise IndependentIdentityError("GTF strand must be '+', '-' or '.'")
                attributes = dict(_attributes(attrs))
                primary = aliases.get(contig, contig)
                if lengths is not None:
                    if primary not in lengths:
                        raise IndependentIdentityError(f"unknown contig {contig!r} in selected reference context")
                    if end > lengths[primary]:
                        raise IndependentIdentityError(f"record end exceeds reference length for {primary!r}")
                if feature not in ("transcript", "exon"):
                    continue
                if strand == ".":
                    raise IndependentIdentityError("transcript/exon strand must be + or -")
                transcript_id = attributes.get("transcript_id")
                if not transcript_id:
                    raise IndependentIdentityError("transcript/exon record lacks transcript_id")
                state = records.setdefault(transcript_id, _TranscriptRows(primary, strand, number))
                if state.contig != primary or state.strand != strand:
                    raise IndependentIdentityError(f"transcript {transcript_id!r} mixes contigs or strands")
                if "gene_id" in attributes:
                    state.genes.add(attributes["gene_id"])
                    if len(state.genes) > 1:
                        raise IndependentIdentityError(f"transcript {transcript_id!r} has conflicting gene_id values")
                if feature == "exon":
                    state.exons.append((start, end))
                else:
                    if state.bounds is not None:
                        raise IndependentIdentityError(f"duplicate transcript meta-feature for {transcript_id!r}")
                    state.bounds = (start, end)
            except (ValueError, OverflowError) as error:
                raise IndependentIdentityError(f"{source}:{number}: {error}") from error
    if not records:
        raise IndependentIdentityError(f"{source}: no transcript/exon models found")
    result = {}
    for transcript_id, state in sorted(records.items()):
        try:
            intervals = _exon_intervals(state.exons)
            if state.bounds is not None and state.bounds != (intervals[0][0], intervals[-1][1]):
                raise IndependentIdentityError("transcript meta-feature bounds differ from its exon span")
        except IndependentIdentityError as error:
            raise IndependentIdentityError(f"{source}:{state.first_line}: transcript {transcript_id!r}: {error}") from error
        result[transcript_id] = IndependentTranscript(
            transcript_id, next(iter(state.genes), None), state.contig, state.strand, intervals)
    return result


def sha256_file(path: str | Path) -> str:
    """Return SHA-256 of exact file bytes, including gzip headers if present."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _source_record(path: str | Path, role: str) -> dict[str, str | int]:
    file = Path(path)
    return {"role": role, "path": str(file), "bytes": file.stat().st_size, "sha256": sha256_file(file)}


def verify_conformance(path: str | Path) -> dict[str, object]:
    """Check canonical bytes and both digests in the existing published vectors."""
    document = json.loads(Path(path).read_text(encoding="utf-8"))
    failures = []
    identities = 0
    for vector in document["vectors"]:
        actual = compute_identity(assembly_fingerprint=document["assembly_fingerprint"],
                                  contig=vector["contig"], strand=vector["strand"], exons=vector["exons"])
        for role in ("form", "splice_chain"):
            expected = vector[role]
            identities += expected is not None
            if actual[role] != expected:
                failures.append({"vector": vector["name"], "role": role,
                                 "expected": expected, "observed": actual[role]})
    root = Path(__file__).resolve().parents[1]
    return {
        "schema": "txid.independent-conformance.v1",
        "calculator_version": REFERENCE_VERSION,
        "python_version": sys.version.split()[0],
        "scope": "Separate standard-library implementation of published exact-v1 rules; no TxID imports. Supplied assembly fingerprint; no independent FASTA validation.",
        "assembly_fingerprint": document["assembly_fingerprint"],
        "vectors": len(document["vectors"]), "identities_checked": identities,
        "valid": not failures, "failures": failures,
        "inputs": [_source_record(path, "golden_vectors"),
                   _source_record(Path(__file__), "independent_calculator"),
                   _source_record(root / "docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md", "specification"),
                   _source_record(root / "schemas/canonical-object.schema.json", "canonical_schema")],
    }


def _aliases_from_tsv(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["alias", "primary"]:
            raise IndependentIdentityError("alias TSV header must be exactly: alias<TAB>primary")
        result = {}
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise IndependentIdentityError(f"{path}:{number}: alias TSV must have two columns")
            if row["alias"] in result:
                raise IndependentIdentityError(f"{path}:{number}: duplicate contig alias")
            result[row["alias"]] = row["primary"]
    return result


def _emit(report: dict[str, object], output: Path | None) -> None:
    payload = json.dumps(report, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    if output is None:
        sys.stdout.write(payload)
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", dir=output.parent,
                                         prefix=".independent-identity-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        # An atomic, exclusive link never replaces a previous result or input.
        os.link(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    conform = commands.add_parser("conformance", help="independently verify existing exact-v1 golden vectors")
    conform.add_argument("--vectors", type=Path, required=True)
    conform.add_argument("--output", type=Path, help="new report path; existing files are never overwritten")
    gtf = commands.add_parser("gtf", help="calculate transcript keys from a strict plain/gzip GTF")
    gtf.add_argument("--input", type=Path, required=True)
    gtf.add_argument("--assembly-fingerprint", required=True)
    gtf.add_argument("--contig-lengths", type=Path, help="JSON object mapping primary names to positive lengths")
    gtf.add_argument("--contig-aliases", type=Path, help="explicit alias/primary TSV; requires --contig-lengths")
    gtf.add_argument("--output", type=Path, help="new report path; default stdout")
    args = parser.parse_args(argv)
    try:
        if args.command == "conformance":
            report = verify_conformance(args.vectors)
        else:
            lengths = None if args.contig_lengths is None else json.loads(args.contig_lengths.read_text(encoding="utf-8"))
            aliases = None if args.contig_aliases is None else _aliases_from_tsv(args.contig_aliases)
            models = read_gtf_models(args.input, contig_aliases=aliases, contig_lengths=lengths)
            rows = []
            for model in models.values():
                rows.append({"transcript_id": model.transcript_id, "gene_id": model.gene_id,
                             "contig": model.contig, "strand": model.strand, "exons": model.exons,
                             **model.identity(args.assembly_fingerprint)})
            sources = [_source_record(args.input, "input_gtf"), _source_record(Path(__file__), "independent_calculator")]
            for source, role in ((args.contig_lengths, "contig_lengths"), (args.contig_aliases, "contig_aliases")):
                if source is not None:
                    sources.append(_source_record(source, role))
            report = {"schema": "txid.independent-gtf-identities.v1", "calculator_version": REFERENCE_VERSION,
                      "python_version": sys.version.split()[0],
                      "assembly_fingerprint": args.assembly_fingerprint,
                      "reference_validation": "supplied contig names and lengths" if lengths is not None else "fingerprint supplied; contig lengths unverified",
                      "scope": "Strict GTF exon structures only; no FASTA fingerprint recalculation, annotation classification or nonidentity attribute round trip.",
                      "transcripts": len(rows), "models": rows, "inputs": sources, "valid": True}
        _emit(report, args.output)
        return 0 if report["valid"] else 1
    except (OSError, ValueError, KeyError, TypeError, EOFError) as error:
        print(f"independent-identity: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
