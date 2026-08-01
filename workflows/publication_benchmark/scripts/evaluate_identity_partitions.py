#!/usr/bin/env python3
"""Evaluate identity partitions on fixed GTF models using raw exon structures."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import sqlite3
import sys
from collections import Counter, defaultdict
from pathlib import Path

from benchmark_lib import clustered_ratio_bootstrap


Member = tuple[str, str, str]
Signature = tuple[object, ...]


def _load_parser(repo: Path):
    sys.path.insert(0, str((repo / "src").resolve()))
    from txid.parser import parse_annotation

    return parse_annotation


def _form_signature(model) -> Signature:
    return (
        model.contig,
        model.strand,
        tuple((exon.start, exon.end) for exon in model.exons),
    )


def _splice_signature(model) -> Signature:
    if model.exon_count == 1:
        return ("single-exon", *_form_signature(model))
    return ("multi-exon", model.contig, model.strand, tuple(model.introns))


def _choose_two(value: int) -> int:
    return value * (value - 1) // 2


def _partition_metrics(
    clusters: dict[str, set[Member]],
    form_by_member: dict[Member, Signature],
    splice_by_member: dict[Member, Signature],
    *,
    bootstrap_replicates: int = 2000,
    bootstrap_seed: int = 20260731,
    allow_missing: bool = False,
) -> dict[str, object]:
    member_cluster: dict[Member, str] = {}
    for cluster, members in clusters.items():
        if not members:
            raise ValueError(f"identity cluster {cluster!r} is empty")
        for member in members:
            if member not in form_by_member:
                raise ValueError(f"identity output contains unknown model: {member}")
            if member in member_cluster:
                raise ValueError(f"model occurs in multiple identity clusters: {member}")
            member_cluster[member] = cluster
    missing = set(form_by_member) - set(member_cluster)
    if missing and not allow_missing:
        example = sorted(missing)[0]
        raise ValueError(
            f"identity output omits {len(missing)} fixed models; first={example}"
        )
    evaluated_forms = {
        member: signature
        for member, signature in form_by_member.items()
        if member in member_cluster
    }

    forms_to_clusters: dict[Signature, set[str]] = defaultdict(set)
    form_cluster_counts: dict[Signature, Counter[str]] = defaultdict(Counter)
    exact_mix = 0
    splice_mix = 0
    end_variant_mix = 0
    false_merge_pairs = 0
    merge_strata = []
    for cluster, members in clusters.items():
        form_counts = Counter(form_by_member[member] for member in members)
        splice_keys = {splice_by_member[member] for member in members}
        for signature, count in form_counts.items():
            forms_to_clusters[signature].add(cluster)
            form_cluster_counts[signature][cluster] = count
        if len(form_counts) > 1:
            exact_mix += 1
            cluster_false_merges = _choose_two(len(members)) - sum(
                _choose_two(count) for count in form_counts.values()
            )
            false_merge_pairs += cluster_false_merges
            if len(splice_keys) == 1:
                end_variant_mix += 1
        else:
            cluster_false_merges = 0
        merge_strata.append(
            (cluster_false_merges, _choose_two(len(members)))
        )
        if len(splice_keys) > 1:
            splice_mix += 1

    false_split_pairs = 0
    split_strata = []
    for signature, cluster_ids in forms_to_clusters.items():
        counts = form_cluster_counts[signature]
        eligible_pairs = _choose_two(sum(counts.values()))
        group_false_splits = eligible_pairs - sum(
            _choose_two(count) for count in counts.values()
        )
        false_split_pairs += group_false_splits
        split_strata.append((group_false_splits, eligible_pairs))
    eligible_same_target_pairs = sum(
        _choose_two(count)
        for count in Counter(evaluated_forms.values()).values()
    )
    eligible_same_cluster_pairs = sum(
        _choose_two(len(members)) for members in clusters.values()
    )
    matrix_rows = len(
        {(member[0], member[1]) for member in evaluated_forms}
    )
    occupied_matrix_cells = sum(
        len({(member[0], member[1]) for member in members})
        for members in clusters.values()
    )
    possible_matrix_cells = matrix_rows * len(clusters)
    merge_interval = clustered_ratio_bootstrap(
        merge_strata,
        replicates=bootstrap_replicates,
        seed=bootstrap_seed,
    )
    split_interval = clustered_ratio_bootstrap(
        split_strata,
        replicates=bootstrap_replicates,
        seed=bootstrap_seed + 1,
    )

    return {
        "input_observations": len(form_by_member),
        "observations": len(evaluated_forms),
        "omitted_observations": len(missing),
        "observations_expected": len(form_by_member),
        "observations_emitted": len(evaluated_forms),
        "observations_missing": len(missing),
        "distinct_exact_forms": len(set(evaluated_forms.values())),
        "identity_clusters": len(clusters),
        "matrix_rows": matrix_rows,
        "matrix_columns": len(clusters),
        "matrix_occupied_cells": occupied_matrix_cells,
        "matrix_sparsity": (
            1.0 - occupied_matrix_cells / possible_matrix_cells
            if possible_matrix_cells
            else 0.0
        ),
        "matrix_columns_in_multiple_rows": sum(
            len({(member[0], member[1]) for member in members}) > 1
            for members in clusters.values()
        ),
        "clusters_mixing_exact_forms": exact_mix,
        "clusters_mixing_splice_structures": splice_mix,
        "clusters_merging_end_variants_only": end_variant_mix,
        "exact_forms_split_across_clusters": sum(
            len(cluster_ids) > 1 for cluster_ids in forms_to_clusters.values()
        ),
        "eligible_same_target_pairs": eligible_same_target_pairs,
        "eligible_same_cluster_pairs": eligible_same_cluster_pairs,
        "observation_pairs_false_merged": false_merge_pairs,
        "observation_pairs_false_split": false_split_pairs,
        "false_merge_rate": (
            false_merge_pairs / eligible_same_cluster_pairs
            if eligible_same_cluster_pairs
            else 0.0
        ),
        "false_split_rate": (
            false_split_pairs / eligible_same_target_pairs
            if eligible_same_target_pairs
            else 0.0
        ),
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "bootstrap_algorithm": (
            "ordinary resampling with replacement of independent structural units"
        ),
        "false_merge_bootstrap_unit": "identity cluster",
        "false_merge_ci95_lower": merge_interval["lower"],
        "false_merge_ci95_upper": merge_interval["upper"],
        "false_merge_bootstrap_effective_replicates": merge_interval[
            "effective_replicates"
        ],
        "false_split_bootstrap_unit": "exact structural target",
        "false_split_ci95_lower": split_interval["lower"],
        "false_split_ci95_upper": split_interval["upper"],
        "false_split_bootstrap_effective_replicates": split_interval[
            "effective_replicates"
        ],
    }


def _load_models(
    raw_inputs: list[str], parse_annotation
) -> tuple[
    dict[Member, Signature],
    dict[Member, Signature],
    dict[Path, tuple[str, str]],
]:
    form_by_member: dict[Member, Signature] = {}
    splice_by_member: dict[Member, Signature] = {}
    labels_by_path: dict[Path, tuple[str, str]] = {}
    for raw in raw_inputs:
        caller, dataset, raw_path = raw.split("=", 2)
        path = Path(raw_path).resolve()
        label = (caller, dataset)
        if path in labels_by_path:
            raise ValueError(f"duplicate model path: {path}")
        labels_by_path[path] = label
        for model in parse_annotation(path):
            member = (caller, dataset, model.original_transcript_id)
            if member in form_by_member:
                raise ValueError(f"duplicate fixed-model key: {member}")
            form_by_member[member] = _form_signature(model)
            splice_by_member[member] = _splice_signature(model)
    return form_by_member, splice_by_member, labels_by_path


def _txid_rows(path: Path) -> dict[Member, dict[str, str]]:
    rows: dict[Member, dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            member = (
                row["tool"],
                row["sample"],
                row["original_transcript_id"],
            )
            if member in rows:
                raise ValueError(f"duplicate TxID observation key: {member}")
            rows[member] = row
    return rows


def _clusters_from_rows(
    rows: dict[Member, dict[str, str]], field: str
) -> dict[str, set[Member]]:
    clusters: dict[str, set[Member]] = defaultdict(set)
    for member, row in rows.items():
        clusters[row[field]].add(member)
    return dict(clusters)


def _txid_clusters(path: Path) -> dict[str, set[Member]]:
    return _clusters_from_rows(_txid_rows(path), "txid_form")


def _txid_order_comparison(sorted_path: Path, shuffled_path: Path) -> dict[str, object]:
    sorted_rows = _txid_rows(sorted_path)
    shuffled_rows = _txid_rows(shuffled_path)
    if set(sorted_rows) != set(shuffled_rows):
        raise ValueError("TxID order runs contain different observations")
    fields = (
        "txid_form",
        "txid_sc",
        "txid_transcript_id",
        "txid_gene_id",
        "classification",
        "gene_candidates",
    )
    return {
        "observations": len(sorted_rows),
        "field_changes": {
            field: sum(
                sorted_rows[member][field] != shuffled_rows[member][field]
                for member in sorted_rows
            )
            for field in fields
        },
        "exact_form_partition": _order_comparison(
            _clusters_from_rows(sorted_rows, "txid_form"),
            _clusters_from_rows(shuffled_rows, "txid_form"),
        ),
        "gene_locus_partition": _order_comparison(
            _clusters_from_rows(sorted_rows, "txid_gene_id"),
            _clusters_from_rows(shuffled_rows, "txid_gene_id"),
        ),
    }


def _tracking_clusters(
    tracking: Path,
    run_json: Path,
    labels_by_path: dict[Path, tuple[str, str]],
) -> dict[str, set[Member]]:
    record = json.loads(run_json.read_text(encoding="utf-8"))
    ordered_labels = []
    for raw_path in record["input_order"]:
        path = Path(raw_path).resolve()
        if path not in labels_by_path:
            raise ValueError(f"gffcompare run contains undeclared input: {path}")
        ordered_labels.append(labels_by_path[path])
    clusters: dict[str, set[Member]] = {}
    with tracking.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) != 4 + len(ordered_labels):
                raise ValueError(
                    f"{tracking}:{line_number}: expected {4 + len(ordered_labels)} "
                    f"columns, found {len(fields)}"
                )
            cluster = fields[0].split("|", 1)[0]
            members: set[Member] = set()
            for query_index, (field, label) in enumerate(
                zip(fields[4:], ordered_labels), 1
            ):
                if field == "-":
                    continue
                for item in field.split(","):
                    prefix, separator, payload = item.partition(":")
                    if not separator or prefix != f"q{query_index}":
                        raise ValueError(
                            f"{tracking}:{line_number}: malformed query field {item!r}"
                        )
                    parts = payload.split("|")
                    if len(parts) < 2:
                        raise ValueError(
                            f"{tracking}:{line_number}: no transcript ID in {item!r}"
                        )
                    members.add((label[0], label[1], parts[1]))
            if cluster in clusters:
                raise ValueError(f"duplicate gffcompare cluster ID: {cluster}")
            clusters[cluster] = members
    return clusters


def _isoseql_clusters(
    database: Path, labels_by_path: dict[Path, tuple[str, str]]
) -> tuple[dict[str, set[Member]], dict[str, set[Member]]]:
    label_by_experiment = {
        f"{dataset}-{caller}": (caller, dataset)
        for caller, dataset in labels_by_path.values()
    }
    exact: dict[str, set[Member]] = defaultdict(set)
    junction: dict[str, set[Member]] = defaultdict(set)
    with sqlite3.connect(database) as connection:
        rows = connection.execute(
            """
            SELECT p.PBID, e.exp_name, p.ends_id, ie.isoform_id
            FROM PBID AS p
            JOIN exp AS e ON e.id = p.exp
            JOIN isoform_ends AS ie ON ie.id = p.ends_id
            ORDER BY p.id
            """
        ).fetchall()
    for transcript_id, experiment, ends_id, isoform_id in rows:
        if experiment not in label_by_experiment:
            raise ValueError(
                f"isoSeQL database contains undeclared experiment {experiment!r}"
            )
        caller, dataset = label_by_experiment[experiment]
        member = (caller, dataset, transcript_id)
        exact[str(ends_id)].add(member)
        junction[str(isoform_id)].add(member)
    return dict(exact), dict(junction)


def _pair_relations(clusters: dict[str, set[Member]]) -> set[tuple[Member, Member]]:
    pairs = set()
    for members in clusters.values():
        for left, right in itertools.combinations(sorted(members), 2):
            pairs.add((left, right))
    return pairs


def _member_clusters(clusters: dict[str, set[Member]]) -> dict[Member, str]:
    return {
        member: cluster for cluster, members in clusters.items() for member in members
    }


def _order_comparison(
    sorted_partition: dict[str, set[Member]],
    shuffled_partition: dict[str, set[Member]],
) -> dict[str, int]:
    sorted_members = _member_clusters(sorted_partition)
    shuffled_members = _member_clusters(shuffled_partition)
    sorted_member_set = set(sorted_members)
    shuffled_member_set = set(shuffled_members)
    shared_members = sorted_member_set & shuffled_member_set
    sorted_cluster_sizes = Counter(
        sorted_members[member] for member in shared_members
    )
    shuffled_cluster_sizes = Counter(
        shuffled_members[member] for member in shared_members
    )
    sorted_pair_count = sum(
        _choose_two(count) for count in sorted_cluster_sizes.values()
    )
    shuffled_pair_count = sum(
        _choose_two(count) for count in shuffled_cluster_sizes.values()
    )
    joint_counts = Counter(
        (sorted_members[member], shuffled_members[member])
        for member in shared_members
    )
    pairs_shared_by_both_partitions = sum(
        _choose_two(count) for count in joint_counts.values()
    )
    return {
        "members_sorted": len(sorted_member_set),
        "members_shuffled": len(shuffled_member_set),
        "members_shared": len(shared_members),
        "members_only_sorted": len(sorted_member_set - shuffled_member_set),
        "members_only_shuffled": len(shuffled_member_set - sorted_member_set),
        "pair_relation_changes": (
            sorted_pair_count
            + shuffled_pair_count
            - 2 * pairs_shared_by_both_partitions
        ),
        "cluster_identifier_changes": sum(
            sorted_members[member] != shuffled_members[member]
            for member in shared_members
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--model", action="append", required=True, help="caller=dataset=GTF"
    )
    parser.add_argument("--txid-observations", type=Path, required=True)
    parser.add_argument("--txid-shuffled-observations", type=Path)
    parser.add_argument("--gffcompare-sorted-tracking", type=Path, required=True)
    parser.add_argument("--gffcompare-sorted-run", type=Path, required=True)
    parser.add_argument("--gffcompare-shuffled-tracking", type=Path, required=True)
    parser.add_argument("--gffcompare-shuffled-run", type=Path, required=True)
    parser.add_argument("--isoseql-sorted-db", type=Path, required=True)
    parser.add_argument("--isoseql-shuffled-db", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260731)
    args = parser.parse_args()

    parse_annotation = _load_parser(args.repo)
    form_by_member, splice_by_member, labels_by_path = _load_models(
        args.model, parse_annotation
    )
    exact_sorted, junction_sorted = _isoseql_clusters(
        args.isoseql_sorted_db, labels_by_path
    )
    exact_shuffled, junction_shuffled = _isoseql_clusters(
        args.isoseql_shuffled_db, labels_by_path
    )
    exact_partitions = {
        "txid-sorted": _txid_clusters(args.txid_observations),
        "gffcompare-sorted": _tracking_clusters(
            args.gffcompare_sorted_tracking,
            args.gffcompare_sorted_run,
            labels_by_path,
        ),
        "gffcompare-shuffled": _tracking_clusters(
            args.gffcompare_shuffled_tracking,
            args.gffcompare_shuffled_run,
            labels_by_path,
        ),
        "isoseql-exact-sorted": exact_sorted,
        "isoseql-exact-shuffled": exact_shuffled,
    }
    evaluations = {}
    for name, value in exact_partitions.items():
        evaluations[name] = _partition_metrics(
            value,
            form_by_member,
            splice_by_member,
            bootstrap_replicates=args.bootstrap_replicates,
            bootstrap_seed=args.bootstrap_seed,
            allow_missing=name.startswith("gffcompare-"),
        )
    junction_partitions = {
        "gffcompare-splice-sorted": exact_partitions["gffcompare-sorted"],
        "gffcompare-splice-shuffled": exact_partitions["gffcompare-shuffled"],
        "isoseql-junction-sorted": junction_sorted,
        "isoseql-junction-shuffled": junction_shuffled,
    }
    evaluations.update(
        {
            name: _partition_metrics(
                value,
                splice_by_member,
                splice_by_member,
                allow_missing=name.startswith("gffcompare-"),
                bootstrap_replicates=args.bootstrap_replicates,
                bootstrap_seed=args.bootstrap_seed,
            )
            for name, value in junction_partitions.items()
        }
    )
    order = {
        "gffcompare_exact": _order_comparison(
            exact_partitions["gffcompare-sorted"],
            exact_partitions["gffcompare-shuffled"],
        ),
        "isoseql_exact_ends": _order_comparison(exact_sorted, exact_shuffled),
        "isoseql_common_junction": _order_comparison(
            junction_sorted, junction_shuffled
        ),
    }
    if args.txid_shuffled_observations is not None:
        order["txid"] = _txid_order_comparison(
            args.txid_observations, args.txid_shuffled_observations
        )

    rows = []
    for name, metrics in evaluations.items():
        system, run = name.rsplit("-", 1)
        row = {"system": system, "run": run}
        row.update(metrics)
        rows.append(row)
    rows.sort(key=lambda row: (row["system"], row["run"]))
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(
        json.dumps(
            {
                "target_definition": (
                    "exact equality of assembly-local contig, strand, and all "
                    "1-based closed exon intervals parsed independently from fixed GTFs"
                ),
                "splice_target_definition": (
                    "equality of assembly-local contig, strand, and ordered intron "
                    "chain; single-exon models retain their exact interval"
                ),
                "bootstrap": {
                    "confidence": 0.95,
                    "algorithm": (
                        "ordinary resampling with replacement of independent "
                        "structural units"
                    ),
                    "merge_unit": "identity cluster",
                    "split_unit": "exact structural target",
                    "replicates": args.bootstrap_replicates,
                    "seed": args.bootstrap_seed,
                },
                "evaluations": evaluations,
                "order_comparisons": order,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
