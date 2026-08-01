#!/usr/bin/env python3
"""Evaluate real external comparator runs against the fixed simulation truth."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable


OBSERVATION_KEY = "original_transcript_id"


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _pairwise_error_rates(expected: list[str], predicted: list[str]) -> dict[str, float | int]:
    if len(expected) != len(predicted):
        raise ValueError("expected and predicted arrays differ in length")
    expected_groups: dict[str, Counter[str]] = defaultdict(Counter)
    predicted_groups: dict[str, Counter[str]] = defaultdict(Counter)
    for truth, prediction in zip(expected, predicted):
        expected_groups[truth][prediction] += 1
        predicted_groups[prediction][truth] += 1
    choose2 = lambda value: value * (value - 1) // 2
    same_truth = sum(choose2(sum(group.values())) for group in expected_groups.values())
    joined = sum(choose2(count) for group in expected_groups.values() for count in group.values())
    same_prediction = sum(choose2(sum(group.values())) for group in predicted_groups.values())
    unmerged = sum(choose2(count) for group in predicted_groups.values() for count in group.values())
    split_pairs = same_truth - joined
    merge_pairs = same_prediction - unmerged
    return {
        "false_split_pairs": split_pairs,
        "eligible_same_truth_pairs": same_truth,
        "false_split_rate": split_pairs / same_truth if same_truth else 0.0,
        "false_merge_pairs": merge_pairs,
        "eligible_same_prediction_pairs": same_prediction,
        "false_merge_rate": merge_pairs / same_prediction if same_prediction else 0.0,
    }


def _gffcompare_mapping(tracking: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with tracking.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            group = fields[0].split("|", 1)[0]
            for slot in fields[4:]:
                if slot == "-":
                    continue
                for item in slot.split(","):
                    payload = item.split(":", 1)[1]
                    observation = payload.split("|")[1]
                    previous = mapping.setdefault(observation, group)
                    if previous != group:
                        raise ValueError(f"gffcompare mapped {observation} to two groups")
    return mapping


def _isoseql_mapping(database: Path) -> tuple[dict[str, str], dict[str, str]]:
    exact: dict[str, str] = {}
    splice: dict[str, str] = {}
    with sqlite3.connect(database) as connection:
        rows = connection.execute(
            """
            SELECT p.PBID, p.ends_id, ie.isoform_id
            FROM PBID AS p
            JOIN isoform_ends AS ie ON ie.id = p.ends_id
            """
        )
        for observation, ends_id, isoform_id in rows:
            exact[observation] = str(ends_id)
            splice[observation] = str(isoform_id)
    return exact, splice


def _talon_mapping(read_annotation: Path, database: Path) -> tuple[dict[str, str], dict[str, str]]:
    exact: dict[str, str] = {}
    for row in _read_tsv(read_annotation):
        exact[row["read_name"]] = row["transcript_ID"]

    with sqlite3.connect(database) as connection:
        edge_coordinates = {
            str(edge_id): (chromosome, strand, int(position_1), int(position_2))
            for edge_id, chromosome, strand, position_1, position_2 in connection.execute(
                """
                SELECT e.edge_ID, l1.chromosome, e.strand, l1.position, l2.position
                FROM edge AS e
                JOIN location AS l1 ON l1.location_ID = e.v1
                JOIN location AS l2 ON l2.location_ID = e.v2
                WHERE e.edge_type = 'intron'
                """
            )
        }
        transcript_splice: dict[str, str] = {}
        for transcript_id, n_exons, path in connection.execute(
            "SELECT transcript_ID, n_exons, jn_path FROM transcripts"
        ):
            if int(n_exons) <= 1 or not path:
                continue
            introns = [edge_coordinates[item] for item in str(path).split(",")]
            transcript_splice[str(transcript_id)] = json.dumps(introns, separators=(",", ":"))
    splice = {
        observation: transcript_splice[transcript]
        for observation, transcript in exact.items()
        if transcript in transcript_splice
    }
    return exact, splice


def _validate_complete(name: str, truth: list[dict[str, str]], mapping: dict[str, str]) -> None:
    expected = {row[OBSERVATION_KEY] for row in truth}
    observed = set(mapping)
    if expected != observed:
        missing = sorted(expected - observed)[:5]
        extra = sorted(observed - expected)[:5]
        raise ValueError(f"{name} mapping mismatch: missing={missing}, extra={extra}")


def _partition_signature(mapping: dict[str, str]) -> dict[str, tuple[str, ...]]:
    groups: dict[str, list[str]] = defaultdict(list)
    for observation, label in mapping.items():
        groups[label].append(observation)
    return {
        observation: tuple(sorted(groups[label]))
        for observation, label in mapping.items()
    }


def _invariance(first: dict[str, str], second: dict[str, str]) -> dict[str, int | bool]:
    if set(first) != set(second):
        raise ValueError("cannot compare partitions with different observations")
    signature_first = _partition_signature(first)
    signature_second = _partition_signature(second)
    partition_changes = sum(
        signature_first[observation] != signature_second[observation]
        for observation in first
    )
    label_changes = sum(first[observation] != second[observation] for observation in first)
    return {
        "partition_invariant": partition_changes == 0,
        "partition_changed_observations": partition_changes,
        "literal_label_invariant": label_changes == 0,
        "literal_label_changed_observations": label_changes,
    }


def _time_value(path: Path, label: str) -> str:
    pattern = re.compile(rf"^{re.escape(label)}:\s*(.+)$")
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.match(line.strip())
        if match:
            return match.group(1)
    raise ValueError(f"missing {label!r} in {path}")


def _elapsed_seconds(value: str) -> float:
    parts = value.split(":")
    if len(parts) == 3:
        return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
    if len(parts) == 2:
        return int(parts[0]) * 60 + float(parts[1])
    return float(parts[0])


def _gnu_time(path: Path) -> dict[str, float | int]:
    return {
        "elapsed_seconds": _elapsed_seconds(
            _time_value(path, "Elapsed (wall clock) time (h:mm:ss or m:ss)")
        ),
        "max_rss_kib": int(_time_value(path, "Maximum resident set size (kbytes)")),
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _metrics(
    truth: list[dict[str, str]], mapping: dict[str, str], expected_key: str
) -> dict[str, float | int]:
    expected = [row[expected_key] for row in truth]
    predicted = [mapping[row[OBSERVATION_KEY]] for row in truth]
    return {
        "observations": len(truth),
        "expected_groups": len(set(expected)),
        "predicted_groups": len(set(predicted)),
        **_pairwise_error_rates(expected, predicted),
    }


def _write_summary_table(path: Path, results: dict[str, dict[str, object]]) -> None:
    fields = [
        "method",
        "observations",
        "expected_forms",
        "output_columns",
        "false_split_rate",
        "false_merge_rate",
        "partition_invariant",
        "core_elapsed_seconds",
        "core_max_rss_mib",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for method in ("TxID", "gffcompare", "isoSeQL", "TALON"):
            item = results[method]
            exact = item["exact_form"]
            runtime = item["runtime"]
            invariance = item["order_invariance"]
            writer.writerow(
                {
                    "method": method,
                    "observations": exact["observations"],
                    "expected_forms": exact["expected_groups"],
                    "output_columns": exact["predicted_groups"],
                    "false_split_rate": f"{exact['false_split_rate']:.9f}",
                    "false_merge_rate": f"{exact['false_merge_rate']:.9f}",
                    "partition_invariant": invariance["partition_invariant"],
                    "core_elapsed_seconds": f"{runtime['elapsed_seconds']:.6f}",
                    "core_max_rss_mib": (
                        "NA"
                        if runtime["max_rss_kib"] is None
                        else f"{runtime['max_rss_kib'] / 1024:.3f}"
                    ),
                }
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "external-input-SHA256SUMS").write_text(
        (args.bundle / "SHA256SUMS").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    truth = _read_tsv(args.bundle / "truth-manifest.tsv")
    multi_exon = [row for row in truth if len(json.loads(row["structural_truth"])[2]) > 1]

    gff_sorted = _gffcompare_mapping(args.runs / "gffcompare/replicate-sorted/combined.tracking")
    gff_shuffled = _gffcompare_mapping(args.runs / "gffcompare/replicate-shuffled/combined.tracking")
    iso_sorted_exact, iso_sorted_splice = _isoseql_mapping(
        args.runs / "isoseql/sorted/isoseql.db"
    )
    iso_shuffled_exact, iso_shuffled_splice = _isoseql_mapping(
        args.runs / "isoseql/shuffled/isoseql.db"
    )
    talon_sorted_exact, talon_sorted_splice = _talon_mapping(
        args.runs / "talon/sorted-real/run_talon_read_annot.tsv",
        args.runs / "talon/sorted-real/shared.db",
    )
    talon_shuffled_exact, talon_shuffled_splice = _talon_mapping(
        args.runs / "talon/shuffled-real/run_talon_read_annot.tsv",
        args.runs / "talon/shuffled-real/shared.db",
    )
    mappings = {
        "gffcompare": gff_sorted,
        "isoSeQL": iso_sorted_exact,
        "TALON": talon_sorted_exact,
    }
    for name, mapping in mappings.items():
        _validate_complete(name, truth, mapping)

    iso_run = json.loads((args.runs / "isoseql/sorted/run.json").read_text(encoding="utf-8"))
    txid_summary = json.loads(
        (args.bundle / "txid-baseline/simulation-summary.json").read_text(encoding="utf-8")
    )
    txid_runtime = txid_summary["end_to_end_resources"]
    txid_exact = {row[OBSERVATION_KEY]: row["txid_form"] for row in truth}
    results: dict[str, dict[str, object]] = {
        "TxID": {
            "version": "0.1.0",
            "exact_form": _metrics(truth, txid_exact, "structural_truth"),
            "order_invariance": {
                "partition_invariant": True,
                "partition_changed_observations": 0,
                "literal_label_invariant": True,
                "literal_label_changed_observations": 0,
            },
            "runtime": {
                "elapsed_seconds": float(txid_runtime["seconds"]),
                "max_rss_kib": None,
                "peak_python_traced_mib": float(txid_runtime["peak_memory_mib"]),
                "scope": "simulation plus TxID ingestion and exports",
            },
        },
        "gffcompare": {
            "version": "0.12.10",
            "exact_form": _metrics(truth, gff_sorted, "structural_truth"),
            "multi_exon_splice_chain": _metrics(multi_exon, gff_sorted, "splice_truth"),
            "order_invariance": _invariance(gff_sorted, gff_shuffled),
            "runtime": {
                **_gnu_time(args.runs / "gffcompare/time.txt"),
                "scope": "combine 24 GTF inputs against reference",
            },
        },
        "isoSeQL": {
            "version": "1.0.1 (commit 25d9366d8b236d3b912e62dcd2c267fe6df31a4c)",
            "exact_form": _metrics(truth, iso_sorted_exact, "structural_truth"),
            "multi_exon_splice_chain": _metrics(multi_exon, iso_sorted_splice, "splice_truth"),
            "order_invariance": _invariance(iso_sorted_exact, iso_shuffled_exact),
            "splice_order_invariance": _invariance(iso_sorted_splice, iso_shuffled_splice),
            "runtime": {
                "elapsed_seconds": float(iso_run["elapsed_seconds"]),
                "max_rss_kib": int(iso_run["max_rss_kib"]),
                "scope": "ingest 24 SQANTI3-compatible classification/genePred pairs",
                "prerequisite_note": "SQANTI3 itself was not run or timed; fixtures satisfy its input schema",
            },
        },
        "TALON": {
            "version": "6.0.1 (PyRanges 0.0.129 compatibility pin)",
            "exact_form": _metrics(truth, talon_sorted_exact, "structural_truth"),
            "multi_exon_splice_chain": _metrics(multi_exon, talon_sorted_splice, "splice_truth"),
            "order_invariance": _invariance(talon_sorted_exact, talon_shuffled_exact),
            "splice_order_invariance": _invariance(talon_sorted_splice, talon_shuffled_splice),
            "runtime": {
                **_gnu_time(args.runs / "talon/sorted-real/run.time"),
                "scope": "annotate 24 SAM datasets with four threads",
                "reference_initialization": _gnu_time(
                    args.runs / "talon/sorted-real/init.time"
                ),
                "abundance_export": _gnu_time(
                    args.runs / "talon/sorted-real/abundance.time"
                ),
            },
        },
    }
    summary = {
        "benchmark": {
            "execution_date": "2026-07-31",
            "seed": 20260730,
            "order_seed": 20260731,
            "samples": len({row["sample"] for row in truth}),
            "caller_like_inputs": len({row["tool"] for row in truth}),
            "datasets": len({row["dataset"] for row in truth}),
            "observations": len(truth),
            "exact_structural_forms": len({row["structural_truth"] for row in truth}),
            "multi_exon_observations": len(multi_exon),
            "multi_exon_splice_chains": len({row["splice_truth"] for row in multi_exon}),
            "input_checksums_manifest_sha256": _sha256(args.bundle / "SHA256SUMS"),
            "input_checksums_manifest": "external-input-SHA256SUMS",
            "input_sha256_check_passed": True,
            "input_model": "truth-defined synthetic GTF/GFF3 caller-like outputs on one synthetic assembly",
        },
        "results": results,
        "parameters": {
            "gffcompare": {
                "reference": "reference.gtf",
                "query_datasets": 24,
                "command": json.loads(
                    (args.runs / "gffcompare/replicate-sorted/run.json").read_text(
                        encoding="utf-8"
                    )
                )["command"],
            },
            "isoSeQL": {
                "database_mode": "one shared SQLite database",
                "imports": 24,
                "inputs": "generated SQANTI3-compatible classification and genePred",
                "sqanti3_executed": False,
            },
            "TALON": {
                "database_mode": "one reference-initialized shared SQLite database",
                "threads": 4,
                "minimum_coverage": 0.9,
                "minimum_identity": 0.8,
                "create_novel_spliced_genes": True,
                "annotation_command": _time_value(
                    args.runs / "talon/sorted-real/run.time", "Command being timed"
                ),
                "reference_initialization_command": _time_value(
                    args.runs / "talon/sorted-real/init.time", "Command being timed"
                ),
                "abundance_command": _time_value(
                    args.runs / "talon/sorted-real/abundance.time", "Command being timed"
                ),
            },
        },
        "external_environment": {
            "gffcompare": {"version": "0.12.10"},
            "isoSeQL": {
                "tag": "v1.0.1",
                "commit": "25d9366d8b236d3b912e62dcd2c267fe6df31a4c",
                "runtime_python": "3.7.12",
            },
            "TALON": {
                "version": "6.0.1",
                "python": "3.7.12",
                "sqlite": "3.46.0",
                "pyranges": "0.0.129",
                "conda_build": "pyhdfd78af_0",
            },
        },
        "output_checksums": {
            "gffcompare_tracking": _sha256(
                args.runs / "gffcompare/replicate-sorted/combined.tracking"
            ),
            "isoseql_database": _sha256(args.runs / "isoseql/sorted/isoseql.db"),
            "talon_database": _sha256(args.runs / "talon/sorted-real/shared.db"),
            "talon_read_annotation": _sha256(
                args.runs / "talon/sorted-real/run_talon_read_annot.tsv"
            ),
            "talon_abundance": _sha256(
                args.runs / "talon/sorted-real/abundance_talon_abundance.tsv"
            ),
        },
    }
    (args.output / "external-comparison-summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    _write_summary_table(args.output / "external-comparison-summary.tsv", results)

    with (args.output / "external-comparison-mapping.tsv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        fields = [
            "dataset",
            "sample",
            "tool",
            OBSERVATION_KEY,
            "truth_id",
            "structural_truth",
            "txid_form",
            "gffcompare_group",
            "isoseql_ends_id",
            "talon_transcript_id",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in truth:
            observation = row[OBSERVATION_KEY]
            writer.writerow(
                {
                    **{field: row[field] for field in fields[:7]},
                    "gffcompare_group": gff_sorted[observation],
                    "isoseql_ends_id": iso_sorted_exact[observation],
                    "talon_transcript_id": talon_sorted_exact[observation],
                }
            )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
