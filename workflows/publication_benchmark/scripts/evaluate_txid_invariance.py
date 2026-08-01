#!/usr/bin/env python3
"""Compare TxID orders/releases and emit a presence-matrix consequence table."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def _read(path: Path) -> list[dict[str, str]]:
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", newline="", encoding="utf-8")
    else:
        handle = path.open("r", newline="", encoding="utf-8")
    with handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _mapping(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], dict[str, str]]:
    result = {}
    for row in rows:
        key = (row["sample"], row["tool"], row["original_transcript_id"])
        if key in result:
            raise ValueError(f"duplicate observation key: {key}")
        result[key] = row
    return result


def _payload_bytes(path: Path) -> bytes:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return handle.read()
    return path.read_bytes()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(_payload_bytes(path))
    return digest.hexdigest()


def _incremental_metrics(
    snapshot_dir: Path, final_catalog: Path
) -> dict[str, object]:
    paths = sorted(snapshot_dir.glob("*.catalog.tsv"))
    if not paths:
        raise ValueError(f"no incremental catalog snapshots in {snapshot_dir}")
    previous: dict[str, dict[str, object]] = {}
    disappearances = 0
    metadata_changes = 0
    classification_regressions = 0
    for path in paths:
        current_rows = _read(path)
        current: dict[str, dict[str, object]] = {}
        for row in current_rows:
            form_id = row["form_id"]
            state = current.setdefault(
                form_id,
                {
                    "splice_chain_id": row["splice_chain_id"],
                    "classifications": set(),
                },
            )
            if state["splice_chain_id"] != row["splice_chain_id"]:
                raise ValueError(
                    f"form_id has conflicting splice chains in {path}: "
                    f"{form_id}"
                )
            classifications = state["classifications"]
            assert isinstance(classifications, set)
            classifications.update(
                filter(None, row["classifications"].split(","))
            )
        disappearances += len(set(previous) - set(current))
        for form_id in set(previous) & set(current):
            old = previous[form_id]
            new = current[form_id]
            if old["splice_chain_id"] != new["splice_chain_id"]:
                metadata_changes += 1
            old_classes = old["classifications"]
            new_classes = new["classifications"]
            assert isinstance(old_classes, set)
            assert isinstance(new_classes, set)
            if not old_classes <= new_classes:
                classification_regressions += 1
        previous = current
    return {
        "incremental_snapshots": len(paths),
        "incremental_identifier_disappearances": disappearances,
        "incremental_structural_metadata_changes": metadata_changes,
        "incremental_classification_regressions": classification_regressions,
        "incremental_final_catalog_bytes_identical": (
            _payload_bytes(paths[-1]) == _payload_bytes(final_catalog)
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sorted-observations", type=Path, required=True)
    parser.add_argument("--shuffled-observations", type=Path, required=True)
    parser.add_argument("--release2-observations", type=Path, required=True)
    parser.add_argument("--sorted-catalog", type=Path, required=True)
    parser.add_argument("--shuffled-catalog", type=Path, required=True)
    parser.add_argument("--sorted-snapshot-dir", type=Path, required=True)
    parser.add_argument(
        "--exclude-matrix-tool",
        action="append",
        default=[],
        help="Tool observations retained for identity checks but excluded from the sample matrix.",
    )
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--summary-tsv", type=Path, required=True)
    parser.add_argument("--summary-json", type=Path, required=True)
    args = parser.parse_args()
    sorted_map = _mapping(_read(args.sorted_observations))
    shuffled_map = _mapping(_read(args.shuffled_observations))
    release2_map = _mapping(_read(args.release2_observations))
    if set(sorted_map) != set(shuffled_map) or set(sorted_map) != set(release2_map):
        raise ValueError("TxID runs do not contain identical observation keys")
    keys = sorted(sorted_map)
    order_form_changes = sum(
        sorted_map[key]["txid_form"] != shuffled_map[key]["txid_form"] for key in keys
    )
    order_partition_changes = sum(
        sorted_map[key]["txid_sc"] != shuffled_map[key]["txid_sc"] for key in keys
    )
    release_form_changes = sum(
        sorted_map[key]["txid_form"] != release2_map[key]["txid_form"] for key in keys
    )
    release_class_changes = sum(
        sorted_map[key]["classification"] != release2_map[key]["classification"]
        for key in keys
    )
    matrix_map = {
        key: row
        for key, row in sorted_map.items()
        if row["tool"] not in set(args.exclude_matrix_tool)
    }
    samples = sorted({key[0] for key in matrix_map})
    forms = sorted({row["txid_form"] for row in matrix_map.values()})
    values: dict[tuple[str, str], int] = defaultdict(int)
    for key, row in matrix_map.items():
        values[(key[0], row["txid_form"])] += 1
    args.matrix.parent.mkdir(parents=True, exist_ok=True)
    with args.matrix.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", *forms])
        for sample in samples:
            writer.writerow([sample, *[values[(sample, form)] for form in forms]])
    raw_columns = len(
        {
            (row["tool"], row["original_transcript_id"])
            for row in matrix_map.values()
        }
    )
    nonzero = sum(value > 0 for value in values.values())
    catalog_hashes = {
        "sorted": _sha256(args.sorted_catalog),
        "shuffled": _sha256(args.shuffled_catalog),
    }
    incremental = _incremental_metrics(
        args.sorted_snapshot_dir, args.sorted_catalog
    )
    summary = {
        "observations": len(keys),
        "matrix_observations": len(matrix_map),
        "samples": len(samples),
        "matrix_excluded_tools": sorted(set(args.exclude_matrix_tool)),
        "raw_upstream_columns": raw_columns,
        "txid_form_columns": len(forms),
        "matrix_nonzero": nonzero,
        "matrix_nonzero_fraction": (
            nonzero / (len(samples) * len(forms)) if samples and forms else 0.0
        ),
        "order_form_identifier_changes": order_form_changes,
        "order_splice_identifier_changes": order_partition_changes,
        "annotation_release_form_identifier_changes": release_form_changes,
        "annotation_release_classification_changes": release_class_changes,
        "catalog_sha256": catalog_hashes,
        "catalog_bytes_identical": (
            _payload_bytes(args.sorted_catalog)
            == _payload_bytes(args.shuffled_catalog)
        ),
        **incremental,
    }
    rows = [
        {
            "metric": key,
            "value": (
                json.dumps(value, sort_keys=True)
                if isinstance(value, (dict, list))
                else value
            ),
        }
        for key, value in summary.items()
    ]
    with args.summary_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.summary_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    incremental_failed = (
        incremental["incremental_identifier_disappearances"]
        or incremental["incremental_structural_metadata_changes"]
        or incremental["incremental_classification_regressions"]
        or not incremental["incremental_final_catalog_bytes_identical"]
    )
    if (
        order_form_changes
        or order_partition_changes
        or release_form_changes
        or incremental_failed
    ):
        raise SystemExit("TxID exact-identity invariance check failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
