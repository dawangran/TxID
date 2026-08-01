#!/usr/bin/env python3
"""Validate publication benchmark dataset and reference manifests."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


DATASET_COLUMNS = {
    "dataset_id",
    "sample_id",
    "replicate",
    "cohort",
    "platform",
    "read_type",
    "input_kind",
    "input_path",
    "reference_id",
    "annotation_id",
    "truth_set",
    "public_accession",
    "source_url",
    "sha256",
    "enabled",
}
REFERENCE_COLUMNS = {
    "reference_id",
    "assembly",
    "annotation_id",
    "annotation_name",
    "fasta",
    "annotation",
    "truth_gtf",
    "source_url",
    "license",
    "fasta_sha256",
    "annotation_sha256",
    "enabled",
}
TRUE = {"1", "true", "yes"}
FALSE = {"0", "false", "no"}


def _read(path: Path, required: set[str]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        missing = sorted(required - fields)
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{path}: manifest is empty")
    return rows


def _enabled(row: dict[str, str], path: Path, number: int) -> bool:
    value = row["enabled"].strip().lower()
    if value not in TRUE | FALSE:
        raise ValueError(f"{path}:{number}: enabled must be true or false")
    return value in TRUE


def _duplicates(values: list[str]) -> list[str]:
    return sorted(key for key, count in Counter(values).items() if count > 1)


def _expand(value: str, *, repo: Path, outdir: Path) -> Path:
    return Path(value.replace("{repo}", str(repo)).replace("{outdir}", str(outdir)))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(
    datasets_path: Path,
    references_path: Path,
    *,
    repo: Path,
    outdir: Path,
    allow_deferred_generated: bool,
) -> dict[str, object]:
    datasets = _read(datasets_path, DATASET_COLUMNS)
    references = _read(references_path, REFERENCE_COLUMNS)
    duplicate_datasets = _duplicates([row["dataset_id"] for row in datasets])
    if duplicate_datasets:
        raise ValueError(f"duplicate dataset_id values: {', '.join(duplicate_datasets)}")
    duplicate_contexts = _duplicates(
        [f"{row['reference_id']}::{row['annotation_id']}" for row in references]
    )
    if duplicate_contexts:
        raise ValueError(
            "duplicate reference/annotation contexts: " + ", ".join(duplicate_contexts)
        )

    enabled_references: dict[tuple[str, str], dict[str, str]] = {}
    for number, row in enumerate(references, 2):
        if not _enabled(row, references_path, number):
            continue
        key = (row["reference_id"], row["annotation_id"])
        enabled_references[key] = row
        for column, checksum_column in (
            ("fasta", "fasta_sha256"),
            ("annotation", "annotation_sha256"),
        ):
            raw = row[column]
            if raw in {"", "NA"}:
                raise ValueError(f"{references_path}:{number}: missing {column}")
            path = _expand(raw, repo=repo, outdir=outdir)
            deferred = allow_deferred_generated and row["source_url"] == "generated"
            if not path.is_file() and not deferred:
                raise ValueError(f"{references_path}:{number}: file not found: {path}")
            expected = row[checksum_column]
            if path.is_file() and expected not in {"generated", "NA", ""}:
                observed = _sha256(path)
                if observed.lower() != expected.lower():
                    raise ValueError(
                        f"{references_path}:{number}: checksum mismatch for {path}"
                    )

    enabled_datasets = []
    for number, row in enumerate(datasets, 2):
        if not _enabled(row, datasets_path, number):
            continue
        enabled_datasets.append(row)
        key = (row["reference_id"], row["annotation_id"])
        if key not in enabled_references:
            raise ValueError(
                f"{datasets_path}:{number}: missing enabled reference context {key}"
            )
        if row["input_kind"] not in {"generated", "local", "local_required"}:
            raise ValueError(
                f"{datasets_path}:{number}: unsupported input_kind {row['input_kind']!r}"
            )
        raw_path = row["input_path"]
        if raw_path in {"", "NA"}:
            raise ValueError(f"{datasets_path}:{number}: missing input_path")
        path = _expand(raw_path, repo=repo, outdir=outdir)
        deferred = allow_deferred_generated and row["input_kind"] == "generated"
        if not path.is_file() and not deferred:
            raise ValueError(f"{datasets_path}:{number}: input not found: {path}")
        expected = row["sha256"]
        if path.is_file() and expected not in {"generated", "NA", ""}:
            if _sha256(path).lower() != expected.lower():
                raise ValueError(
                    f"{datasets_path}:{number}: checksum mismatch for {path}"
                )
        if row["input_kind"] != "generated":
            if row["public_accession"] in {"", "NA"}:
                raise ValueError(
                    f"{datasets_path}:{number}: public/local data require an accession"
                )
            if row["source_url"] in {"", "NA"}:
                raise ValueError(
                    f"{datasets_path}:{number}: public/local data require a source URL"
                )
            if row["sha256"] in {"", "NA", "generated"}:
                raise ValueError(
                    f"{datasets_path}:{number}: staged public data require SHA-256"
                )

    if not enabled_datasets:
        raise ValueError("no enabled datasets")
    return {
        "dataset_manifest": str(datasets_path.resolve()),
        "reference_manifest": str(references_path.resolve()),
        "enabled_datasets": sorted(row["dataset_id"] for row in enabled_datasets),
        "enabled_reference_contexts": sorted(
            f"{reference_id}::{annotation_id}"
            for reference_id, annotation_id in enabled_references
        ),
        "status": "valid",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-deferred-generated", action="store_true")
    args = parser.parse_args()
    result = validate(
        args.datasets,
        args.references,
        repo=args.repo.resolve(),
        outdir=args.outdir.resolve(),
        allow_deferred_generated=args.allow_deferred_generated,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
