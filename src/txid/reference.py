"""Reference FASTA and explicit contig alias handling."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import TextIO

from .errors import ReferenceError
from .models import ContigRecord, SequenceCollection, TranscriptModel


def _open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("r", encoding="utf-8", newline=None)


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def load_aliases(path: str | Path | None, primary_names: set[str]) -> tuple[tuple[str, str], ...]:
    if path is None:
        return ()
    aliases: dict[str, str] = {}
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if line_number == 1 and [item.lower() for item in columns[:2]] == ["alias", "primary"]:
                continue
            if len(columns) != 2 or not all(columns):
                raise ReferenceError(
                    f"alias table line {line_number}: expected non-empty alias and primary columns"
                )
            alias, primary = columns
            if primary not in primary_names:
                raise ReferenceError(
                    f"alias table line {line_number}: target {primary!r} is absent from FASTA"
                )
            if alias in primary_names and alias != primary:
                raise ReferenceError(
                    f"alias table line {line_number}: primary contig {alias!r} cannot alias {primary!r}"
                )
            previous = aliases.get(alias)
            if previous is not None and previous != primary:
                raise ReferenceError(
                    f"alias {alias!r} maps to both {previous!r} and {primary!r}"
                )
            aliases[alias] = primary
    return tuple(sorted(aliases.items()))


def read_fasta(
    path: str | Path, *, assembly_name: str, alias_path: str | Path | None = None
) -> SequenceCollection:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ReferenceError(f"FASTA line {line_number}: empty sequence name")
                if name in sequences:
                    raise ReferenceError(f"FASTA line {line_number}: duplicate sequence name {name!r}")
                sequences[name] = []
                current = name
            elif current is None:
                raise ReferenceError(f"FASTA line {line_number}: sequence occurs before a header")
            else:
                sequence = "".join(line.split()).upper()
                if not sequence.isascii():
                    raise ReferenceError(f"FASTA line {line_number}: sequence must be ASCII")
                sequences[current].append(sequence)
    if not sequences:
        raise ReferenceError("reference FASTA contains no sequences")

    contigs: list[ContigRecord] = []
    fingerprint_rows: list[str] = []
    for name in sorted(sequences):
        sequence = "".join(sequences[name])
        if not sequence:
            raise ReferenceError(f"FASTA sequence {name!r} is empty")
        sequence_digest = hashlib.sha256(sequence.encode("ascii")).hexdigest()
        contigs.append(ContigRecord(name, len(sequence), sequence_digest))
        fingerprint_rows.append(f"{name}\t{len(sequence)}\t{sequence_digest}\n")
    fingerprint = "sha256:" + hashlib.sha256(
        "".join(fingerprint_rows).encode("utf-8")
    ).hexdigest()
    aliases = load_aliases(alias_path, set(sequences))
    return SequenceCollection(
        name=assembly_name,
        fingerprint=fingerprint,
        contigs=tuple(contigs),
        aliases=aliases,
    )


def resolve_and_validate(
    transcript: TranscriptModel, reference: SequenceCollection
) -> TranscriptModel:
    contig = reference.resolve(transcript.contig)
    if transcript.end > contig.length:
        raise ReferenceError(
            f"transcript {transcript.original_transcript_id!r} ends at {transcript.end}, "
            f"beyond contig {contig.name!r} length {contig.length}"
        )
    return transcript.with_contig(contig.name)

