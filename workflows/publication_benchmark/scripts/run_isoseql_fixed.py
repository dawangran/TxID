#!/usr/bin/env python3
"""Run isoSeQL over actual SQANTI3 outputs in a declared dataset order."""

from __future__ import annotations

import argparse
import csv
import random
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


def _parse(values: list[str]) -> list[tuple[str, str, Path, Path]]:
    rows = []
    for value in values:
        fields = value.split("=", 3)
        if len(fields) != 4:
            raise ValueError(
                f"invalid --input {value!r}; expected caller=dataset=classification=genePred"
            )
        rows.append((fields[0], fields[1], Path(fields[2]), Path(fields[3])))
    return rows


def _configs(root: Path, caller: str, dataset: str) -> tuple[Path, Path]:
    sample = root / f"{dataset}-{caller}.sample.csv"
    experiment = root / f"{dataset}-{caller}.experiment.csv"
    with sample.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sample_name", "tissue", "disease", "age", "sex"])
        writer.writerow([dataset, "benchmark", "NA", "NA", "NA"])
    with experiment.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "RIN",
                "date",
                "platform",
                "method",
                "vMap",
                "vReference",
                "vAnnot",
                "vLima",
                "vCCS",
                "vIsoseq3",
                "vCupcake",
                "vSQANTI",
                "exp_name",
            ]
        )
        writer.writerow(
            [
                "NA",
                "2026-07-31",
                "benchmark",
                caller,
                "minimap2-2.31",
                "manifest",
                "manifest",
                "NA",
                "NA",
                "NA",
                "NA",
                "6.0.1",
                f"{dataset}-{caller}",
            ]
        )
    return sample, experiment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--order", required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument(
        "--source-audit",
        type=Path,
        help="Audit JSON for any explicitly disclosed compatibility source tree.",
    )
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    rows = sorted(_parse(args.input), key=lambda row: (row[1], row[0]))
    seed = None
    if args.order.startswith("shuffled-"):
        seed = int(args.order.split("-", 1)[1])
        random.Random(seed).shuffle(rows)
    elif args.order != "sorted":
        raise ValueError(f"unsupported order: {args.order}")
    if args.database.exists():
        raise FileExistsError(f"refusing to overwrite database: {args.database}")
    commands = []
    for index, (caller, dataset, classification, genepred) in enumerate(rows, 1):
        sample_config, experiment_config = _configs(args.work_dir, caller, dataset)
        commands.append(
            run_command(
                [
                    "python",
                    str(args.source / "isoSeQL_run.py"),
                    "--classif",
                    classification,
                    "--genePred",
                    genepred,
                    "--sampleConfig",
                    sample_config,
                    "--expConfig",
                    experiment_config,
                    "--db",
                    args.database,
                ],
                log_prefix=args.work_dir / f"{index:03d}-{dataset}-{caller}",
            )
        )
    write_json(
        args.run_json,
        {
            "stage": "isoseql_fixed_models",
            "version_contract": "1.0.1",
            "source_commit": "25d9366d8b236d3b912e62dcd2c267fe6df31a4c",
            "source_audit": (
                None if args.source_audit is None else file_records([args.source_audit])[0]
            ),
            "sqanti3_included": True,
            "order": args.order,
            "order_seed": seed,
            "input_order": [
                {"caller": caller, "dataset": dataset}
                for caller, dataset, _, _ in rows
            ],
            "commands": commands,
            "inputs": file_records(
                [
                    path
                    for _, _, classification, genepred in rows
                    for path in (classification, genepred)
                ]
                + ([] if args.source_audit is None else [args.source_audit])
            ),
            "outputs": file_records([args.database]),
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
