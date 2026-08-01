#!/usr/bin/env python3
"""Evaluate TALON read partitions against exact-form and splice-chain truth."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from benchmark_lib import clustered_ratio_bootstrap


Member = tuple[str, str]


def _read_truth(
    read_truth_path: Path, transcript_truth_path: Path
) -> tuple[dict[Member, str], dict[Member, tuple[object, ...]]]:
    transcript_splice: dict[str, tuple[object, ...]] = {}
    with transcript_truth_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            introns = json.loads(row["introns"])
            if introns:
                signature = (
                    "multi-exon",
                    row["contig"],
                    row["strand"],
                    tuple(tuple(value) for value in introns),
                )
            else:
                signature = (
                    "single-exon",
                    row["contig"],
                    row["strand"],
                    tuple(tuple(value) for value in json.loads(row["exons"])),
                )
            transcript_splice[row["truth_id"]] = signature

    exact_by_read: dict[Member, str] = {}
    splice_by_read: dict[Member, tuple[object, ...]] = {}
    with read_truth_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            member = (row["sample"], row["read_id"])
            truth_id = row["truth_id"]
            if member in exact_by_read:
                raise ValueError(f"duplicate read truth key: {member}")
            if truth_id not in transcript_splice:
                raise ValueError(f"read truth references unknown transcript: {truth_id}")
            exact_by_read[member] = truth_id
            splice_by_read[member] = transcript_splice[truth_id]
    return exact_by_read, splice_by_read


def _talon_clusters(database: Path) -> dict[str, set[Member]]:
    clusters: dict[str, set[Member]] = defaultdict(set)
    seen: set[Member] = set()
    with sqlite3.connect(database) as connection:
        rows = connection.execute(
            """
            SELECT dataset, read_name, transcript_ID
            FROM observed
            ORDER BY dataset, read_name, obs_ID
            """
        ).fetchall()
    for dataset, read_name, transcript_id in rows:
        member = (dataset, read_name)
        if member in seen:
            raise ValueError(f"TALON database repeats read observation: {member}")
        seen.add(member)
        clusters[str(transcript_id)].add(member)
    return dict(clusters)


def _partition_metrics(
    clusters: dict[str, set[Member]],
    truth_by_member: dict[Member, object],
    *,
    allow_missing_truth: bool = False,
    bootstrap_replicates: int = 2000,
    bootstrap_seed: int = 20260731,
) -> dict[str, object]:
    cluster_by_member: dict[Member, str] = {}
    for cluster_id, members in clusters.items():
        for member in members:
            if member not in truth_by_member:
                raise ValueError(f"TALON database contains unknown read: {member}")
            if member in cluster_by_member:
                raise ValueError(f"read occurs in multiple TALON clusters: {member}")
            cluster_by_member[member] = cluster_id
    missing = set(truth_by_member) - set(cluster_by_member)
    if missing and not allow_missing_truth:
        raise ValueError(
            f"TALON database omits {len(missing)} truth reads; first={sorted(missing)[0]}"
        )
    evaluated_truth = {
        member: truth_by_member[member] for member in cluster_by_member
    }

    false_merge_pairs = 0
    merge_strata = []
    truth_clusters: dict[object, Counter[str]] = defaultdict(Counter)
    mixed_clusters = 0
    for cluster_id, members in clusters.items():
        counts = Counter(evaluated_truth[member] for member in members)
        if len(counts) > 1:
            mixed_clusters += 1
        total = len(members)
        eligible_pairs = total * (total - 1) // 2
        cluster_false_merges = eligible_pairs - sum(
            count * (count - 1) // 2 for count in counts.values()
        )
        false_merge_pairs += cluster_false_merges
        merge_strata.append((cluster_false_merges, eligible_pairs))
        for truth_value, count in counts.items():
            truth_clusters[truth_value][cluster_id] += count

    false_split_pairs = 0
    split_strata = []
    split_truth_groups = 0
    for counts in truth_clusters.values():
        if len(counts) > 1:
            split_truth_groups += 1
        total = sum(counts.values())
        eligible_pairs = total * (total - 1) // 2
        group_false_splits = eligible_pairs - sum(
            count * (count - 1) // 2 for count in counts.values()
        )
        false_split_pairs += group_false_splits
        split_strata.append((group_false_splits, eligible_pairs))
    eligible_same_truth_pairs = sum(
        sum(counts.values()) * (sum(counts.values()) - 1) // 2
        for counts in truth_clusters.values()
    )
    eligible_same_cluster_pairs = sum(
        len(members) * (len(members) - 1) // 2
        for members in clusters.values()
    )
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
        "observations": len(evaluated_truth),
        "truth_observations_requested": len(truth_by_member),
        "truth_observations_missing": len(missing),
        "truth_groups": len(set(evaluated_truth.values())),
        "identity_clusters": len(clusters),
        "clusters_mixing_truth_groups": mixed_clusters,
        "truth_groups_split_across_clusters": split_truth_groups,
        "eligible_same_truth_pairs": eligible_same_truth_pairs,
        "eligible_same_cluster_pairs": eligible_same_cluster_pairs,
        "observation_pairs_false_merged": false_merge_pairs,
        "observation_pairs_false_split": false_split_pairs,
        "false_merge_rate": (
            false_merge_pairs / eligible_same_cluster_pairs
            if eligible_same_cluster_pairs
            else 0.0
        ),
        "false_split_rate": (
            false_split_pairs / eligible_same_truth_pairs
            if eligible_same_truth_pairs
            else 0.0
        ),
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "false_merge_bootstrap_unit": "TALON transcript cluster",
        "false_merge_ci95_lower": merge_interval["lower"],
        "false_merge_ci95_upper": merge_interval["upper"],
        "false_merge_bootstrap_effective_replicates": merge_interval[
            "effective_replicates"
        ],
        "false_split_bootstrap_unit": "truth transcript group",
        "false_split_ci95_lower": split_interval["lower"],
        "false_split_ci95_upper": split_interval["upper"],
        "false_split_bootstrap_effective_replicates": split_interval[
            "effective_replicates"
        ],
    }


def _pair_relations(
    clusters: dict[str, set[Member]]
) -> set[tuple[Member, Member]]:
    return {
        pair
        for members in clusters.values()
        for pair in itertools.combinations(sorted(members), 2)
    }


def _order_comparison(
    sorted_clusters: dict[str, set[Member]],
    shuffled_clusters: dict[str, set[Member]],
    *,
    allow_member_difference: bool = False,
) -> dict[str, int]:
    sorted_member_cluster = {
        member: cluster for cluster, members in sorted_clusters.items() for member in members
    }
    shuffled_member_cluster = {
        member: cluster
        for cluster, members in shuffled_clusters.items()
        for member in members
    }
    sorted_members = set(sorted_member_cluster)
    shuffled_members = set(shuffled_member_cluster)
    if sorted_members != shuffled_members and not allow_member_difference:
        raise ValueError("TALON order runs contain different read observations")
    common = sorted_members & shuffled_members
    sorted_common = {
        cluster: members & common
        for cluster, members in sorted_clusters.items()
        if members & common
    }
    shuffled_common = {
        cluster: members & common
        for cluster, members in shuffled_clusters.items()
        if members & common
    }
    return {
        "pair_relation_changes": len(
            _pair_relations(sorted_common) ^ _pair_relations(shuffled_common)
        ),
        "cluster_identifier_changes": sum(
            sorted_member_cluster[member] != shuffled_member_cluster[member]
            for member in common
        ),
        "common_observations": len(common),
        "observation_set_symmetric_difference": len(
            sorted_members ^ shuffled_members
        ),
    }


def _truth_subset(
    clusters: dict[str, set[Member]], truth_members: set[Member]
) -> tuple[dict[str, set[Member]], int]:
    database_members = {
        member for members in clusters.values() for member in members
    }
    restricted = {
        cluster: members & truth_members
        for cluster, members in clusters.items()
        if members & truth_members
    }
    return restricted, len(database_members - truth_members)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--read-truth", type=Path, required=True)
    parser.add_argument("--transcript-truth", type=Path, required=True)
    parser.add_argument("--sorted-db", type=Path, required=True)
    parser.add_argument("--shuffled-db", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260731)
    parser.add_argument(
        "--truth-subset",
        action="store_true",
        help="ignore database reads absent from the audited truth table",
    )
    parser.add_argument(
        "--allow-missing-truth-reads",
        action="store_true",
        help="evaluate observed truth reads and report truth reads omitted by TALON",
    )
    args = parser.parse_args()

    exact_truth, splice_truth = _read_truth(
        args.read_truth, args.transcript_truth
    )
    sorted_clusters = _talon_clusters(args.sorted_db)
    shuffled_clusters = _talon_clusters(args.shuffled_db)
    ignored_unlabelled = {"sorted": 0, "shuffled": 0}
    if args.truth_subset:
        truth_members = set(exact_truth)
        sorted_clusters, ignored_unlabelled["sorted"] = _truth_subset(
            sorted_clusters, truth_members
        )
        shuffled_clusters, ignored_unlabelled["shuffled"] = _truth_subset(
            shuffled_clusters, truth_members
        )
    evaluations = {}
    rows = []
    for order, clusters in (
        ("sorted", sorted_clusters),
        ("shuffled", shuffled_clusters),
    ):
        for target, truth in (
            ("exact-form", exact_truth),
            ("splice-chain", splice_truth),
        ):
            metrics = _partition_metrics(
                clusters,
                truth,
                allow_missing_truth=args.allow_missing_truth_reads,
                bootstrap_replicates=args.bootstrap_replicates,
                bootstrap_seed=args.bootstrap_seed,
            )
            metrics["ignored_unlabelled_database_reads"] = ignored_unlabelled[order]
            evaluations[f"{target}-{order}"] = metrics
            rows.append({"target": target, "run": order, **metrics})
    rows.sort(key=lambda row: (row["target"], row["run"]))
    order_comparison = _order_comparison(
        sorted_clusters,
        shuffled_clusters,
        allow_member_difference=args.allow_missing_truth_reads,
    )

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
                "evaluation_unit": "individual aligned reads",
                "exact_target": "truth transcript form",
                "splice_target": (
                    "ordered intron chain; single-exon truth retains its exact interval"
                ),
                "truth_subset": args.truth_subset,
                "allow_missing_truth_reads": args.allow_missing_truth_reads,
                "bootstrap": {
                    "confidence": 0.95,
                    "merge_unit": "TALON transcript cluster",
                    "split_unit": "truth transcript group",
                    "replicates": args.bootstrap_replicates,
                    "seed": args.bootstrap_seed,
                },
                "evaluations": evaluations,
                "order_comparison": order_comparison,
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
