#!/usr/bin/env python3
"""Run gffcompare reproducibly over a prepared comparison bundle."""

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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--gffcompare", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--order", choices=("sorted", "shuffled"), default="sorted")
    args = parser.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    prefix = args.output / "combined"
    tracking = prefix.with_suffix(".tracking")
    if tracking.exists():
        raise SystemExit(f"refusing to overwrite existing output: {tracking}")

    with (args.bundle / "inputs-manifest.tsv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows.sort(key=lambda row: row["dataset"])
    if args.order == "shuffled":
        random.Random(ORDER_SEED).shuffle(rows)

    command = [
        str(args.gffcompare.resolve()),
        "-r",
        str((args.bundle / "reference.gtf").resolve()),
        "-o",
        str(prefix.resolve()),
        *[row["query_gtf"] for row in rows],
    ]
    usage_before = resource.getrusage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    elapsed = time.perf_counter() - started
    usage_after = resource.getrusage(resource.RUSAGE_CHILDREN)
    (args.output / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (args.output / "stderr.log").write_text(result.stderr, encoding="utf-8")
    summary = {
        "tool": "gffcompare",
        "order": args.order,
        "order_seed": ORDER_SEED if args.order == "shuffled" else None,
        "datasets": len(rows),
        "dataset_order": [row["dataset"] for row in rows],
        "elapsed_seconds": elapsed,
        "child_user_seconds": usage_after.ru_utime - usage_before.ru_utime,
        "child_system_seconds": usage_after.ru_stime - usage_before.ru_stime,
        "max_rss_kib": usage_after.ru_maxrss,
        "returncode": result.returncode,
        "command": command,
    }
    (args.output / "run.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    if result.returncode != 0:
        raise SystemExit(f"gffcompare failed with exit {result.returncode}")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
