#!/usr/bin/env python3
"""Evaluate caller transcript models against exact truth structures."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import Counter
from pathlib import Path

from benchmark_lib import binary_precision_recall_bootstrap


def _load_parser(repo: Path):
    sys.path.insert(0, str((repo / "src").resolve()))
    from txid.parser import parse_annotation

    return parse_annotation


def _form_key(model) -> tuple[object, ...]:
    return (
        model.contig,
        model.strand,
        tuple((exon.start, exon.end) for exon in model.exons),
    )


def _splice_key(model) -> tuple[object, ...] | None:
    if model.exon_count == 1:
        return None
    return model.contig, model.strand, tuple(model.introns)


def _rates(
    truth: set[object],
    predictions: list[object],
    *,
    bootstrap_replicates: int,
    bootstrap_seed: int,
) -> dict[str, object]:
    unique = set(predictions)
    matched = unique & truth
    precision = len(matched) / len(unique) if unique else 0.0
    recall = len(matched) / len(truth) if truth else 0.0
    result = {
        "truth": len(truth),
        "predicted_observations": len(predictions),
        "predicted_unique": len(unique),
        "matched_unique": len(matched),
        "precision": precision,
        "recall": recall,
        "f1": 2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0,
        "redundant_predictions": len(predictions) - len(unique),
    }
    result.update(
        binary_precision_recall_bootstrap(
            [
                value in unique
                for value in sorted(truth, key=repr)
            ],
            [
                value in truth
                for value in sorted(unique, key=repr)
            ],
            replicates=bootstrap_replicates,
            seed=bootstrap_seed,
        )
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--truth", type=Path, required=True)
    parser.add_argument("--input", action="append", required=True, help="caller=dataset=GTF")
    parser.add_argument(
        "--truth-scope",
        choices=("all_predictions", "truth_contigs"),
        default="all_predictions",
        help=(
            "Evaluate all predictions, or only predictions on contigs represented "
            "in the truth GTF (required for spike-in truth embedded in a host genome)."
        ),
    )
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=2000)
    parser.add_argument("--bootstrap-seed", type=int, default=20260731)
    args = parser.parse_args()
    parse_annotation = _load_parser(args.repo)
    truth_models = parse_annotation(args.truth)
    truth_contigs = {model.contig for model in truth_models}
    truth_forms = {_form_key(model) for model in truth_models}
    truth_splices = {
        key for model in truth_models if (key := _splice_key(model)) is not None
    }
    truth_by_splice = {
        _splice_key(model): model
        for model in truth_models
        if _splice_key(model) is not None
    }
    rows = []
    details = []
    for raw in args.input:
        caller, dataset, raw_path = raw.split("=", 2)
        path = Path(raw_path)
        all_models = parse_annotation(path)
        models = (
            [model for model in all_models if model.contig in truth_contigs]
            if args.truth_scope == "truth_contigs"
            else all_models
        )
        form_keys = [_form_key(model) for model in models]
        splice_keys = [
            key for model in models if (key := _splice_key(model)) is not None
        ]
        form_metrics = _rates(
            truth_forms,
            form_keys,
            bootstrap_replicates=args.bootstrap_replicates,
            bootstrap_seed=args.bootstrap_seed,
        )
        splice_metrics = _rates(
            truth_splices,
            splice_keys,
            bootstrap_replicates=args.bootstrap_replicates,
            bootstrap_seed=args.bootstrap_seed + 1,
        )
        end_errors = []
        for model in models:
            key = _splice_key(model)
            if key is None or key not in truth_by_splice:
                continue
            matched = truth_by_splice[key]
            end_errors.extend([abs(model.tss - matched.tss), abs(model.tes - matched.tes)])
        row = {
            "caller": caller,
            "dataset": dataset,
            "predicted_models": len(models),
            "predicted_models_total": len(all_models),
            "exact_precision": f"{form_metrics['precision']:.6f}",
            "exact_recall": f"{form_metrics['recall']:.6f}",
            "exact_f1": f"{form_metrics['f1']:.6f}",
            "exact_precision_ci95_lower": (
                f"{form_metrics['precision_ci95_lower']:.6f}"
            ),
            "exact_precision_ci95_upper": (
                f"{form_metrics['precision_ci95_upper']:.6f}"
            ),
            "exact_recall_ci95_lower": (
                f"{form_metrics['recall_ci95_lower']:.6f}"
            ),
            "exact_recall_ci95_upper": (
                f"{form_metrics['recall_ci95_upper']:.6f}"
            ),
            "exact_f1_ci95_lower": f"{form_metrics['f1_ci95_lower']:.6f}",
            "exact_f1_ci95_upper": f"{form_metrics['f1_ci95_upper']:.6f}",
            "splice_precision": f"{splice_metrics['precision']:.6f}",
            "splice_recall": f"{splice_metrics['recall']:.6f}",
            "splice_f1": f"{splice_metrics['f1']:.6f}",
            "splice_precision_ci95_lower": (
                f"{splice_metrics['precision_ci95_lower']:.6f}"
            ),
            "splice_precision_ci95_upper": (
                f"{splice_metrics['precision_ci95_upper']:.6f}"
            ),
            "splice_recall_ci95_lower": (
                f"{splice_metrics['recall_ci95_lower']:.6f}"
            ),
            "splice_recall_ci95_upper": (
                f"{splice_metrics['recall_ci95_upper']:.6f}"
            ),
            "splice_f1_ci95_lower": f"{splice_metrics['f1_ci95_lower']:.6f}",
            "splice_f1_ci95_upper": f"{splice_metrics['f1_ci95_upper']:.6f}",
            "median_end_error_bp": (
                f"{statistics.median(end_errors):.3f}" if end_errors else "NA"
            ),
            "redundant_exact_predictions": form_metrics["redundant_predictions"],
        }
        rows.append(row)
        details.append(
            {
                "caller": caller,
                "dataset": dataset,
                "path": str(path.resolve()),
                "truth_scope": args.truth_scope,
                "predicted_models_total": len(all_models),
                "predicted_models_evaluated": len(models),
                "exact": form_metrics,
                "splice_chain": splice_metrics,
                "end_errors_bp": {
                    "n": len(end_errors),
                    "median": statistics.median(end_errors) if end_errors else None,
                    "maximum": max(end_errors) if end_errors else None,
                },
                "single_exon_predictions": sum(model.exon_count == 1 for model in models),
                "exon_count_distribution": dict(
                    sorted(Counter(model.exon_count for model in models).items())
                ),
            }
        )
    rows.sort(key=lambda row: (row["caller"], row["dataset"]))
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    args.output_json.write_text(
        json.dumps(
            {
                "truth_gtf": str(args.truth.resolve()),
                "truth_scope": args.truth_scope,
                "truth_contigs": sorted(truth_contigs),
                "truth_exact_forms": len(truth_forms),
                "truth_multi_exon_splice_chains": len(truth_splices),
                "bootstrap": {
                    "confidence": 0.95,
                    "exact_form_unit": "truth or predicted exact transcript form",
                    "splice_chain_unit": "truth or predicted splice chain",
                    "replicates": args.bootstrap_replicates,
                    "seed": args.bootstrap_seed,
                },
                "evaluations": details,
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
