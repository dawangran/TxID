"""High-level workflows shared by the CLI and tests."""

from __future__ import annotations

import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .classify import ReferenceIndex, build_reference_index, classify
from .errors import RegistryError, TxIDError
from .identity import identify
from .models import Assignment, IdentityBundle, SequenceCollection, TranscriptModel
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


def _validate_fuzzy_options(
    splice_tolerance: int | None, end_tolerance: int | None
) -> None:
    if (splice_tolerance is None) != (end_tolerance is None):
        raise TxIDError("fuzzy mode requires both splice and end tolerances")
    if splice_tolerance is not None and splice_tolerance < 0:
        raise TxIDError("fuzzy splice tolerance must be non-negative")
    if end_tolerance is not None and end_tolerance < 0:
        raise TxIDError("fuzzy end tolerance must be non-negative")


def _stage_annotation(
    path: str | Path,
    reference,
    *,
    fmt: str | None = None,
) -> list[tuple[TranscriptModel, IdentityBundle]]:
    parsed = parse_annotation(path, fmt=fmt)
    contig_map = reference.contig_map()
    alias_map = reference.alias_map()
    staged = []
    for item in parsed:
        resolved = resolve_and_validate(
            item,
            reference,
            contig_map=contig_map,
            alias_map=alias_map,
        )
        staged.append((resolved, identify(resolved, reference.fingerprint)))
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


def _import_with_registry(
    registry: Registry,
    *,
    reference: SequenceCollection,
    annotation_cache: dict[str, tuple[int, ReferenceIndex]],
    input_path: str | Path,
    input_checksum: str,
    sample: str,
    tool: str,
    annotation_name: str,
    annotation_format: str | None,
    fuzzy_splice_tolerance: int | None,
    fuzzy_end_tolerance: int | None,
) -> ImportResult:
    annotation_id, reference_index = _annotation_import_context(
        registry, annotation_cache, annotation_name
    )

    existing = registry.find_import(annotation_id, sample, tool, input_checksum)
    if existing is not None:
        return ImportResult(
            existing,
            tuple(registry.load_assignments(existing, annotation_name)),
            False,
            input_checksum,
        )

    staged = _stage_annotation(input_path, reference, fmt=annotation_format)
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
        annotation_id=annotation_id,
        sample=sample,
        tool=tool,
        source_label=str(input_path),
        input_checksum=input_checksum,
        options={
            "annotation_format": annotation_format,
            "fuzzy_splice_tolerance": fuzzy_splice_tolerance,
            "fuzzy_end_tolerance": fuzzy_end_tolerance,
        },
        software_version=__version__,
        provisional=provisional,
        annotation_name=annotation_name,
    )
    return ImportResult(import_id, tuple(assignments), created, input_checksum)


def _annotation_import_context(
    registry: Registry,
    annotation_cache: dict[str, tuple[int, ReferenceIndex]],
    annotation_name: str,
) -> tuple[int, ReferenceIndex]:
    context = annotation_cache.get(annotation_name)
    if context is None:
        annotation_row = registry.annotation(annotation_name)
        annotation_id = int(annotation_row["id"])
        context = (
            annotation_id,
            build_reference_index(registry.references(annotation_id)),
        )
        annotation_cache[annotation_name] = context
    return context


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
    _validate_fuzzy_options(fuzzy_splice_tolerance, fuzzy_end_tolerance)
    checksum = file_sha256(input_path)
    with Registry(database) as registry:
        reference = registry.reference()
        result = _import_with_registry(
            registry,
            reference=reference,
            annotation_cache={},
            input_path=input_path,
            input_checksum=checksum,
            sample=sample,
            tool=tool,
            annotation_name=annotation_name,
            annotation_format=annotation_format,
            fuzzy_splice_tolerance=fuzzy_splice_tolerance,
            fuzzy_end_tolerance=fuzzy_end_tolerance,
        )
        if fuzzy_splice_tolerance is not None and fuzzy_end_tolerance is not None:
            registry.run_fuzzy(
                splice_tolerance=fuzzy_splice_tolerance,
                end_tolerance=fuzzy_end_tolerance,
                software_version=__version__,
            )
            result = ImportResult(
                result.import_id,
                tuple(registry.load_assignments(result.import_id, annotation_name)),
                result.created,
                checksum,
            )
    write_gtf(output_gtf, result.assignments)
    write_mapping(mapping_path, result.assignments, sample=sample, tool=tool)
    return result


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._")
    return slug or "unnamed"


def _sample_name_from_path(path: str | Path) -> str:
    name = Path(path).name
    if name.lower().endswith(".gz"):
        name = name[:-3]
    for suffix in (".gff3", ".gff", ".gtf"):
        if name.lower().endswith(suffix):
            name = name[: -len(suffix)]
            break
    if not name:
        raise TxIDError(
            f"cannot derive a sample name from input {path!s}; supply --samples"
        )
    return name


def batch_import(
    database: str | Path,
    *,
    manifest: str | Path,
    output_dir: str | Path,
    fuzzy_splice_tolerance: int | None = None,
    fuzzy_end_tolerance: int | None = None,
) -> list[dict[str, Any]]:
    _validate_fuzzy_options(fuzzy_splice_tolerance, fuzzy_end_tolerance)
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
    return _batch_import_rows(
        database,
        rows=rows,
        output_dir=output_dir,
        fuzzy_splice_tolerance=fuzzy_splice_tolerance,
        fuzzy_end_tolerance=fuzzy_end_tolerance,
    )


def multi_import(
    database: str | Path,
    *,
    inputs: list[str | Path],
    tool: str,
    annotation_name: str,
    output_dir: str | Path,
    samples: list[str] | None = None,
    annotation_format: str | None = None,
    fuzzy_splice_tolerance: int | None = None,
    fuzzy_end_tolerance: int | None = None,
) -> list[dict[str, Any]]:
    _validate_fuzzy_options(fuzzy_splice_tolerance, fuzzy_end_tolerance)
    if not inputs:
        raise TxIDError("multi-add requires at least one input")
    if not tool or not annotation_name:
        raise TxIDError("tool and annotation name must be non-empty")
    if samples is not None and len(samples) != len(inputs):
        raise TxIDError(
            "--samples must contain exactly one sample name for each --input path"
        )

    sample_names = (
        [_sample_name_from_path(path) for path in inputs]
        if samples is None
        else list(samples)
    )
    rows: list[dict[str, str]] = []
    for index, (input_path, sample) in enumerate(zip(inputs, sample_names), 1):
        if not sample:
            raise TxIDError(f"multi-add sample {index} is empty")
        path = Path(input_path)
        if not path.is_file():
            raise TxIDError(f"multi-add input {index} does not exist: {path}")
        rows.append(
            {
                "input": str(path),
                "sample": sample,
                "tool": tool,
                "annotation_name": annotation_name,
                "format": annotation_format or "",
            }
        )

    return _batch_import_rows(
        database,
        rows=rows,
        output_dir=output_dir,
        fuzzy_splice_tolerance=fuzzy_splice_tolerance,
        fuzzy_end_tolerance=fuzzy_end_tolerance,
    )


def _batch_import_rows(
    database: str | Path,
    *,
    rows: list[dict[str, str]],
    output_dir: str | Path,
    fuzzy_splice_tolerance: int | None,
    fuzzy_end_tolerance: int | None,
) -> list[dict[str, Any]]:
    _validate_fuzzy_options(fuzzy_splice_tolerance, fuzzy_end_tolerance)
    rows = [dict(row) for row in rows]
    rows.sort(key=lambda row: (row["annotation_name"], row["sample"], row["tool"], row["input"]))
    destination = Path(output_dir)
    prepared: list[tuple[dict[str, str], str, Path, Path]] = []
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
        ensure_distinct_paths(row["input"], [output_gtf, mapping])
        prepared.append((row, checksum, output_gtf, mapping))

    def summarize(
        row: dict[str, str], result: ImportResult, output_gtf: Path, mapping: Path
    ) -> dict[str, Any]:
        return {
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

    results: list[dict[str, Any]] = []
    fuzzy_enabled = fuzzy_splice_tolerance is not None
    pending: list[tuple[dict[str, str], ImportResult, Path, Path]] = []
    with Registry(database) as registry:
        reference = registry.reference()
        annotation_cache: dict[str, tuple[int, ReferenceIndex]] = {}
        for row, checksum, output_gtf, mapping in prepared:
            result = _import_with_registry(
                registry,
                reference=reference,
                annotation_cache=annotation_cache,
                input_path=row["input"],
                input_checksum=checksum,
                sample=row["sample"],
                tool=row["tool"],
                annotation_name=row["annotation_name"],
                annotation_format=row.get("format") or None,
                fuzzy_splice_tolerance=fuzzy_splice_tolerance,
                fuzzy_end_tolerance=fuzzy_end_tolerance,
            )
            if fuzzy_enabled:
                pending.append((row, result, output_gtf, mapping))
            else:
                write_gtf(output_gtf, result.assignments)
                write_mapping(
                    mapping,
                    result.assignments,
                    sample=row["sample"],
                    tool=row["tool"],
                )
                results.append(summarize(row, result, output_gtf, mapping))

        if fuzzy_enabled:
            assert fuzzy_splice_tolerance is not None
            assert fuzzy_end_tolerance is not None
            registry.run_fuzzy(
                splice_tolerance=fuzzy_splice_tolerance,
                end_tolerance=fuzzy_end_tolerance,
                software_version=__version__,
            )
            for row, stored, output_gtf, mapping in pending:
                result = ImportResult(
                    stored.import_id,
                    tuple(
                        registry.load_assignments(
                            stored.import_id, row["annotation_name"]
                        )
                    ),
                    stored.created,
                    stored.input_checksum,
                )
                write_gtf(output_gtf, result.assignments)
                write_mapping(
                    mapping,
                    result.assignments,
                    sample=row["sample"],
                    tool=row["tool"],
                )
                results.append(summarize(row, result, output_gtf, mapping))
    return results
