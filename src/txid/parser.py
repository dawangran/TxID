"""GTF and GFF3 parsing into the shared TxID transcript model."""

from __future__ import annotations

import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Iterable, TextIO
from urllib.parse import unquote

from .errors import AnnotationParseError
from .models import Exon, TranscriptModel, freeze_attributes


def _open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("r", encoding="utf-8", newline=None)


def _split_gtf_fields(text: str, line_number: int) -> list[str]:
    fields: list[str] = []
    current: list[str] = []
    quoted = False
    escaped = False
    for char in text:
        if escaped:
            current.append(char)
            escaped = False
        elif char == "\\" and quoted:
            current.append(char)
            escaped = True
        elif char == '"':
            current.append(char)
            quoted = not quoted
        elif char == ";" and not quoted:
            fields.append("".join(current).strip())
            current = []
        else:
            current.append(char)
    if quoted or escaped:
        raise AnnotationParseError(f"line {line_number}: unterminated quoted GTF attribute")
    if "".join(current).strip():
        fields.append("".join(current).strip())
    return fields


def parse_gtf_attribute_items(
    text: str, line_number: int = 0
) -> list[tuple[str, str]]:
    """Parse ordered GTF attributes without collapsing repeated keys."""

    result: list[tuple[str, str]] = []
    for field in _split_gtf_fields(text, line_number):
        if not field:
            continue
        parts = field.split(None, 1)
        if len(parts) != 2:
            raise AnnotationParseError(
                f"line {line_number}: malformed GTF attribute {field!r}"
            )
        key, raw = parts
        raw = raw.strip()
        if raw.startswith('"'):
            if len(raw) < 2 or not raw.endswith('"'):
                raise AnnotationParseError(
                    f"line {line_number}: malformed quoted GTF value for {key!r}"
                )
            raw = raw[1:-1]
            value_chars: list[str] = []
            escaped = False
            for char in raw:
                if escaped:
                    value_chars.append(char)
                    escaped = False
                elif char == "\\":
                    escaped = True
                else:
                    value_chars.append(char)
            if escaped:
                raise AnnotationParseError(
                    f"line {line_number}: trailing escape in GTF value for {key!r}"
                )
            value = "".join(value_chars)
        else:
            value = raw
        result.append((key, value))
    return result


def parse_gtf_attributes(text: str, line_number: int = 0) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, value in parse_gtf_attribute_items(text, line_number):
        if key in result:
            raise AnnotationParseError(
                f"line {line_number}: duplicate GTF attribute key {key!r}"
            )
        result[key] = value
    return result


def parse_gff3_attributes(text: str, line_number: int = 0) -> dict[str, str]:
    if text == "." or not text:
        return {}
    result: dict[str, str] = {}
    for field in text.split(";"):
        if not field:
            continue
        if "=" not in field:
            raise AnnotationParseError(
                f"line {line_number}: malformed GFF3 attribute {field!r}"
            )
        key, value = field.split("=", 1)
        if not key:
            raise AnnotationParseError(f"line {line_number}: empty GFF3 attribute key")
        decoded_key = unquote(key)
        if decoded_key in result:
            raise AnnotationParseError(
                f"line {line_number}: duplicate GFF3 attribute key {decoded_key!r}"
            )
        result[decoded_key] = unquote(value)
    return result


def _detect_format(lines: Iterable[str], path: str | Path | None = None) -> str:
    if path is not None:
        suffixes = Path(path).suffixes
        useful = suffixes[-2] if suffixes and suffixes[-1] == ".gz" and len(suffixes) > 1 else suffixes[-1] if suffixes else ""
        if useful.lower() in {".gff", ".gff3"}:
            return "gff3"
        if useful.lower() == ".gtf":
            return "gtf"
    for line in lines:
        if line.startswith("##gff-version"):
            return "gff3"
        if line.startswith("#") or not line.strip():
            continue
        columns = line.rstrip("\n").split("\t")
        if len(columns) == 9 and "=" in columns[8] and 'transcript_id "' not in columns[8]:
            return "gff3"
        return "gtf"
    raise AnnotationParseError("annotation contains no feature records")


def parse_annotation(path: str | Path, fmt: str | None = None) -> list[TranscriptModel]:
    with _open_text(path) as handle:
        return parse_annotation_lines(handle, fmt=fmt, path=path)


def parse_annotation_lines(
    lines: Iterable[str], *, fmt: str | None = None, path: str | Path | None = None
) -> list[TranscriptModel]:
    materialized = list(lines)
    annotation_format = (fmt or _detect_format(materialized, path)).lower()
    if annotation_format not in {"gtf", "gff3"}:
        raise AnnotationParseError(f"unsupported annotation format {annotation_format!r}")

    records: list[tuple[int, list[str], dict[str, str]]] = []
    for line_number, raw_line in enumerate(materialized, 1):
        line = raw_line.rstrip("\r\n")
        if not line or line.startswith("#"):
            continue
        columns = line.split("\t")
        if len(columns) != 9:
            raise AnnotationParseError(
                f"line {line_number}: expected 9 tab-delimited columns, found {len(columns)}"
            )
        try:
            start = int(columns[3])
            end = int(columns[4])
        except ValueError as error:
            raise AnnotationParseError(
                f"line {line_number}: start/end must be integers"
            ) from error
        if start < 1 or end < start:
            raise AnnotationParseError(
                f"line {line_number}: invalid 1-based closed interval {start}-{end}"
            )
        if columns[6] not in {"+", "-", "."}:
            raise AnnotationParseError(
                f"line {line_number}: unsupported strand {columns[6]!r}"
            )
        attrs = (
            parse_gtf_attributes(columns[8], line_number)
            if annotation_format == "gtf"
            else parse_gff3_attributes(columns[8], line_number)
        )
        records.append((line_number, columns, attrs))

    transcript_meta: dict[str, tuple[list[str], dict[str, str]]] = {}
    transcript_gene: dict[str, str | None] = {}
    gene_ids: set[str] = set()
    if annotation_format == "gff3":
        for line_number, columns, attrs in records:
            feature = columns[2].lower()
            if feature == "gene":
                gene_id = attrs.get("ID")
                if gene_id:
                    gene_ids.add(gene_id)
            elif feature in {"mrna", "transcript", "lnc_rna", "ncrna", "rrna", "trna"}:
                transcript_id = attrs.get("ID")
                if not transcript_id:
                    raise AnnotationParseError(
                        f"line {line_number}: GFF3 transcript feature has no ID"
                    )
                transcript_meta[transcript_id] = (columns, attrs)
                parents = [item for item in attrs.get("Parent", "").split(",") if item]
                transcript_gene[transcript_id] = parents[0] if len(parents) == 1 else None

    exon_groups: dict[str, list[tuple[int, list[str], dict[str, str]]]] = defaultdict(list)
    for line_number, columns, attrs in records:
        if columns[2].lower() != "exon":
            if annotation_format == "gtf" and columns[2].lower() in {"transcript", "mrna"}:
                transcript_id = attrs.get("transcript_id")
                if not transcript_id:
                    raise AnnotationParseError(
                        f"line {line_number}: GTF transcript feature has no transcript_id"
                    )
                transcript_meta[transcript_id] = (columns, attrs)
                transcript_gene[transcript_id] = attrs.get("gene_id")
            continue
        if annotation_format == "gtf":
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                raise AnnotationParseError(
                    f"line {line_number}: GTF exon has no transcript_id"
                )
            exon_groups[transcript_id].append((line_number, columns, attrs))
            gene_id = attrs.get("gene_id")
            previous = transcript_gene.get(transcript_id)
            if previous is not None and gene_id is not None and previous != gene_id:
                raise AnnotationParseError(
                    f"line {line_number}: transcript {transcript_id!r} has conflicting gene_id values"
                )
            transcript_gene[transcript_id] = previous or gene_id
        else:
            parents = [item for item in attrs.get("Parent", "").split(",") if item]
            if not parents:
                raise AnnotationParseError(
                    f"line {line_number}: GFF3 exon has no Parent"
                )
            for transcript_id in parents:
                exon_groups[transcript_id].append((line_number, columns, attrs))
                if transcript_id not in transcript_gene:
                    gene_id = attrs.get("gene_id")
                    transcript_gene[transcript_id] = gene_id

    models: list[TranscriptModel] = []
    for transcript_id, exon_records in exon_groups.items():
        contigs = {record[1][0] for record in exon_records}
        strands = {record[1][6] for record in exon_records}
        if len(contigs) != 1 or len(strands) != 1 or "." in strands:
            line_numbers = ", ".join(str(record[0]) for record in exon_records)
            raise AnnotationParseError(
                f"transcript {transcript_id!r} has inconsistent contig/strand at lines {line_numbers}"
            )
        exons = tuple(
            Exon(
                start=int(columns[3]),
                end=int(columns[4]),
                source=columns[1],
                score=columns[5],
                phase=columns[7],
                attributes=freeze_attributes(attrs),
            )
            for _, columns, attrs in exon_records
        )
        meta = transcript_meta.get(transcript_id)
        if meta is None:
            attributes: dict[str, str] = {}
            source = exon_records[0][1][1]
        else:
            attributes = meta[1]
            source = meta[0][1]
            if meta[0][0] not in contigs or meta[0][6] not in strands:
                raise AnnotationParseError(
                    f"transcript feature {transcript_id!r} conflicts with its exons"
                )
        models.append(
            TranscriptModel(
                contig=next(iter(contigs)),
                strand=next(iter(strands)),
                exons=exons,
                original_transcript_id=transcript_id,
                original_gene_id=transcript_gene.get(transcript_id),
                source=source,
                source_format=annotation_format,
                attributes=freeze_attributes(attributes),
            )
        )
    if not models:
        raise AnnotationParseError("annotation contains no exon-derived transcript models")
    return sorted(
        models,
        key=lambda model: (
            model.contig,
            model.start,
            model.end,
            model.strand,
            model.original_transcript_id,
        ),
    )


def annotation_fingerprint(transcripts: Iterable[TranscriptModel]) -> str:
    rows = [
        {
            "contig": transcript.contig,
            "exons": [[exon.start, exon.end] for exon in transcript.exons],
            "gene_id": transcript.original_gene_id,
            "strand": transcript.strand,
            "transcript_id": transcript.original_transcript_id,
        }
        for transcript in transcripts
    ]
    rows.sort(
        key=lambda row: (
            row["contig"],
            row["exons"],
            row["strand"],
            row["gene_id"] or "",
            row["transcript_id"],
        )
    )
    payload = json.dumps(
        rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return "sha256:" + hashlib.sha256(payload).hexdigest()
