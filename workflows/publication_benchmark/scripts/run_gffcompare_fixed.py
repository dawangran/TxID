#!/usr/bin/env python3
"""Run gffcompare on fixed caller GTFs in a declared input order."""

from __future__ import annotations

import argparse
import random
import shutil
import subprocess
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--input", action="append", type=Path, required=True)
    parser.add_argument("--order", required=True)
    parser.add_argument("--tracking", type=Path, required=True)
    parser.add_argument("--combined-gtf", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    inputs = sorted(args.input, key=lambda path: str(path))
    seed = None
    if args.order.startswith("shuffled-"):
        seed = int(args.order.split("-", 1)[1])
        random.Random(seed).shuffle(inputs)
    elif args.order != "sorted":
        raise ValueError(f"unsupported order: {args.order}")
    prefix = args.work_dir / "combined"
    record = run_command(
        [
            "gffcompare",
            "-r",
            str(args.reference),
            "-o",
            str(prefix),
            *[str(path) for path in inputs],
        ],
        log_prefix=args.work_dir / "gffcompare",
    )
    source_tracking = prefix.with_suffix(".tracking")
    source_combined = Path(f"{prefix}.combined.gtf")
    if not source_tracking.is_file() or not source_combined.is_file():
        raise FileNotFoundError("gffcompare did not emit tracking and combined GTF")
    args.tracking.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_tracking, args.tracking)
    shutil.copyfile(source_combined, args.combined_gtf)
    version = subprocess.run(
        ["gffcompare", "--version"], text=True, capture_output=True, check=False
    )
    record.update(
        {
            "stage": "gffcompare_fixed_models",
            "version": (version.stdout or version.stderr).strip(),
            "order": args.order,
            "order_seed": seed,
            "input_order": [str(path.resolve()) for path in inputs],
            "inputs": file_records([args.reference, *inputs]),
            "outputs": file_records([args.tracking, args.combined_gtf]),
        }
    )
    write_json(args.run_json, record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
