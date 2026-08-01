#!/usr/bin/env python3
"""Run isoSeQL ingestion reproducibly over a prepared comparison bundle."""

from __future__ import annotations

import argparse
import csv
import json
import random
import resource
import subprocess
import time
from pathlib import Path


ORDER_SEED = 20260731


def _datasets(bundle: Path, order: str) -> list[dict[str, str]]:
    with (bundle / "inputs-manifest.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda row: row["dataset"])
    if order == "shuffled":
        random.Random(ORDER_SEED).shuffle(rows)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--isoseql-source", type=Path, required=True)
    parser.add_argument("--python", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--order", choices=("sorted", "shuffled"), default="sorted")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    database = args.output / "isoseql.db"
    if database.exists():
        raise SystemExit(f"refusing to overwrite existing database: {database}")

    rows = _datasets(args.bundle.resolve(), args.order)
    commands: list[dict[str, object]] = []
    started = time.perf_counter()
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    for index, row in enumerate(rows, start=1):
        input_dir = Path(row["isoseql_dir"])
        command = [
            str(args.python.resolve()),
            str((args.isoseql_source / "isoSeQL_run.py").resolve()),
            "--classif",
            str(input_dir / "classification.tsv"),
            "--genePred",
            str(input_dir / "models.genePred"),
            "--sampleConfig",
            str(input_dir / "sample.config"),
            "--expConfig",
            str(input_dir / "experiment.config"),
            "--db",
            str(database.resolve()),
        ]
        dataset_started = time.perf_counter()
        result = subprocess.run(command, text=True, capture_output=True, check=False)
        elapsed = time.perf_counter() - dataset_started
        log_path = args.output / f"{index:02d}-{row['dataset']}.log"
        log_path.write_text(
            result.stdout + "\n--- stderr ---\n" + result.stderr,
            encoding="utf-8",
        )
        commands.append(
            {
                "dataset": row["dataset"],
                "command": command,
                "elapsed_seconds": elapsed,
                "returncode": result.returncode,
                "log": log_path.name,
            }
        )
        if result.returncode != 0:
            (args.output / "commands.json").write_text(
                json.dumps(commands, indent=2) + "\n", encoding="utf-8"
            )
            raise SystemExit(
                f"isoSeQL failed for {row['dataset']} with exit {result.returncode}; "
                f"see {log_path}"
            )

    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    summary = {
        "tool": "isoSeQL",
        "order": args.order,
        "order_seed": ORDER_SEED if args.order == "shuffled" else None,
        "datasets": len(rows),
        "elapsed_seconds": time.perf_counter() - started,
        "child_user_seconds": usage_after.ru_utime - usage_before.ru_utime,
        "child_system_seconds": usage_after.ru_stime - usage_before.ru_stime,
        "max_rss_kib": usage_after.ru_maxrss,
        "database": str(database.resolve()),
        "commands": commands,
    }
    (args.output / "run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in summary.items() if key != "commands"}, indent=2))


if __name__ == "__main__":
    main()
