"""High-level workflows shared by the CLI and tests."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .classify import build_reference_index, classify
from .errors import RegistryError, TxIDError
from .identity import identify
from .models import Assignment, IdentityBundle, TranscriptModel
from .parser import annotation_fingerprint, parse_annotation
from .reference import file_sha256, read_fasta, resolve_and_validate
from .registry import Registry
from .writers import ensure_distinct_paths, write_gtf, write_mapping


@dataclass(frozen=True, slots=True)
class ImportResult:
    import_id: int
    assignments: tuple[Assignment, ...]
    created: bool
    input_checksum: str


def _stage_annotation(
    path: str | Path,
    reference,
    *,
    fmt: str | None = None,
) -> list[tuple[TranscriptModel, IdentityBundle]]:
    parsed = parse_annotation(path, fmt=fmt)
    resolved = [resolve_and_validate(item, reference) for item in parsed]
    staged = [(item, identify(item, reference.fingerprint)) for item in resolved]
    return sorted(
        staged,
        key=lambda item: (
            item[0].contig,
            item[0].start,
            item[0].end,
            item[0].strand,
            item[1].form.public_id,
            item[0].original_transcript_id,
        ),
    )


def initialize_registry(
    database: str | Path,
    *,
    fasta: str | Path,
    assembly_name: str,
    annotation: str | Path,
    annotation_name: str,
    alias_path: str | Path | None = None,
    annotation_format: str | None = None,
) -> dict[str, Any]:
    reference = read_fasta(fasta, assembly_name=assembly_name, alias_path=alias_path)
    staged = _stage_annotation(annotation, reference, fmt=annotation_format)
    fingerprint = annotation_fingerprint(item[0] for item in staged)
    registry = Registry.create(
        database,
        reference,
        fasta_checksum=file_sha256(fasta),
        alias_checksum=file_sha256(alias_path) if alias_path else None,
        software_version=__version__,
    )
    try:
        with registry.transaction():
            annotation_id = registry.add_annotation(
                name=annotation_name,
                fingerprint=fingerprint,
                input_checksum=file_sha256(annotation),
                transcripts=staged,
            )
    except Exception:
        registry.close()
        target = Path(database)
        if target.exists():
            target.unlink()
        raise
    finally:
        if registry.connection:
            registry.close()
    return {
        "registry": str(database),
        "assembly": assembly_name,
        "assembly_fingerprint": reference.fingerprint,
        "annotation": annotation_name,
        "annotation_fingerprint": fingerprint,
        "reference_transcripts": len(staged),
        "schema_version": 1,
    }


def add_annotation_context(
    database: str | Path,
    *,
    annotation: str | Path,
    annotation_name: str,
    annotation_format: str | None = None,
) -> dict[str, Any]:
    with Registry(database) as registry:
        reference = registry.reference()
        staged = _stage_annotation(annotation, reference, fmt=annotation_format)
        fingerprint = annotation_fingerprint(item[0] for item in staged)
        with registry.transaction():
            annotation_id = registry.add_annotation(
                name=annotation_name,
                fingerprint=fingerprint,
                input_checksum=file_sha256(annotation),
                transcripts=staged,
            )
    return {
        "annotation_id": annotation_id,
        "annotation": annotation_name,
        "annotation_fingerprint": fingerprint,
        "reference_transcripts": len(staged),
    }


def import_file(
    database: str | Path,
    *,
    input_path: str | Path,
    sample: str,
    tool: str,
    annotation_name: str,
    output_gtf: str | Path,
    mapping_path: str | Path,
    annotation_format: str | None = None,
    fuzzy_splice_tolerance: int | None = None,
    fuzzy_end_tolerance: int | None = None,
) -> ImportResult:
    if not sample or not tool:
        raise TxIDError("sample and tool must be non-empty")
    ensure_distinct_paths(input_path, [output_gtf, mapping_path])
    if (fuzzy_splice_tolerance is None) != (fuzzy_end_tolerance is None):
        raise TxIDError("fuzzy mode requires both splice and end tolerances")
    checksum = file_sha256(input_path)
    with Registry(database) as registry:
        reference = registry.reference()
        annotation_row = registry.annotation(annotation_name)
        staged = _stage_annotation(input_path, reference, fmt=annotation_format)
        references = registry.references(int(annotation_row["id"]))
        reference_index = build_reference_index(references)
        provisional = []
        for transcript, identities in staged:
            result = classify(transcript, identities, reference_index)
            provisional.append(
                (
                    transcript,
                    identities,
                    result.label,
                    result.output_gene_id,
                    result.output_transcript_id,
                    result.gene_candidates,
                )
            )
        import_id, assignments, created = registry.store_import(
            annotation_id=int(annotation_row["id"]),
            sample=sample,
            tool=tool,
            source_label=str(input_path),
            input_checksum=checksum,
            options={
                "annotation_format": annotation_format,
                "fuzzy_splice_tolerance": fuzzy_splice_tolerance,
                "fuzzy_end_tolerance": fuzzy_end_tolerance,
            },
            software_version=__version__,
            provisional=provisional,
            annotation_name=annotation_name,
        )
        if fuzzy_splice_tolerance is not None and fuzzy_end_tolerance is not None:
            registry.run_fuzzy(
                splice_tolerance=fuzzy_splice_tolerance,
                end_tolerance=fuzzy_end_tolerance,
                software_version=__version__,
            )
            assignments = registry.load_assignments(import_id, annotation_name)
    write_gtf(output_gtf, assignments)
    write_mapping(mapping_path, assignments, sample=sample, tool=tool)
    return ImportResult(import_id, tuple(assignments), created, checksum)


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return slug or "unnamed"


def batch_import(
    database: str | Path,
    *,
    manifest: str | Path,
    output_dir: str | Path,
    fuzzy_splice_tolerance: int | None = None,
    fuzzy_end_tolerance: int | None = None,
) -> list[dict[str, Any]]:
    with Path(manifest).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"input", "sample", "tool", "annotation_name"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise TxIDError(
                "batch manifest must have tab-delimited columns: input, sample, tool, annotation_name"
            )
        rows = [dict(row) for row in reader]
    if not rows:
        raise TxIDError("batch manifest has no data rows")
    for index, row in enumerate(rows, 2):
        if not all(row.get(key) for key in required):
            raise TxIDError(f"batch manifest line {index} has an empty required field")
        if not Path(row["input"]).is_file():
            raise TxIDError(f"batch manifest line {index}: input does not exist: {row['input']}")
    rows.sort(key=lambda row: (row["annotation_name"], row["sample"], row["tool"], row["input"]))
    destination = Path(output_dir)
    results: list[dict[str, Any]] = []
    for row in rows:
        checksum = file_sha256(row["input"])
        stem = ".".join(
            [
                _slug(row["sample"]),
                _slug(row["tool"]),
                _slug(row["annotation_name"]),
                checksum.split(":", 1)[1][:10],
            ]
        )
        output_gtf = destination / f"{stem}.txid.gtf"
        mapping = destination / f"{stem}.mapping.tsv"
        result = import_file(
            database,
            input_path=row["input"],
            sample=row["sample"],
            tool=row["tool"],
            annotation_name=row["annotation_name"],
            output_gtf=output_gtf,
            mapping_path=mapping,
            annotation_format=row.get("format") or None,
            fuzzy_splice_tolerance=fuzzy_splice_tolerance,
            fuzzy_end_tolerance=fuzzy_end_tolerance,
        )
        results.append(
            {
                "input": row["input"],
                "sample": row["sample"],
                "tool": row["tool"],
                "annotation_name": row["annotation_name"],
                "import_id": result.import_id,
                "created": result.created,
                "transcripts": len(result.assignments),
                "output_gtf": str(output_gtf),
                "mapping": str(mapping),
            }
        )
    return results
