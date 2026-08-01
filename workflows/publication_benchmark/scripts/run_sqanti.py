#!/usr/bin/env python3
"""Run full SQANTI3 QC preprocessing needed by isoSeQL."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


REQUIRED_CLASSIFICATION = {
    "isoform",
    "structural_category",
    "associated_gene",
    "associated_transcript",
    "ref_exons",
    "exons",
    "subcategory",
    "all_canonical",
    "FL",
    "length",
}


def _find_one(root: Path, pattern: str) -> Path:
    paths = sorted(path for path in root.rglob(pattern) if path.is_file())
    if len(paths) != 1:
        raise FileNotFoundError(f"expected one {pattern} under {root}, found {paths}")
    return paths[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--isoforms", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    parser.add_argument("--classification", type=Path, required=True)
    parser.add_argument("--genepred", type=Path, required=True)
    parser.add_argument("--corrected-gtf", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    command = [
        "python",
        str(args.source / "sqanti3_qc.py"),
        "--isoforms",
        str(args.isoforms),
        "--refGTF",
        str(args.annotation),
        "--refFasta",
        str(args.reference),
        "--output",
        args.prefix,
        "--dir",
        str(args.work_dir),
        "--report",
        "skip",
        "--cpus",
        str(args.threads),
    ]
    record = run_command(
        command,
        log_prefix=args.work_dir / "sqanti3",
        env={"MPLCONFIGDIR": str((args.work_dir / "mpl").resolve())},
    )
    source_classification = _find_one(args.work_dir, "*_classification.txt")
    source_genepred = _find_one(args.work_dir, "*_corrected.genePred")
    source_gtf = _find_one(args.work_dir, "*_corrected.gtf")
    with source_classification.open(newline="", encoding="utf-8") as handle:
        fields = set(csv.DictReader(handle, delimiter="\t").fieldnames or [])
    missing = sorted(REQUIRED_CLASSIFICATION - fields)
    if missing:
        raise ValueError(
            "SQANTI3 output is not directly compatible with isoSeQL 1.0.1; "
            f"missing columns: {', '.join(missing)}"
        )
    for source, destination in (
        (source_classification, args.classification),
        (source_genepred, args.genepred),
        (source_gtf, args.corrected_gtf),
    ):
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
    record.update(
        {
            "stage": "sqanti3_preprocessing",
            "version_contract": "6.0.1",
            "inputs": file_records([args.isoforms, args.annotation, args.reference]),
            "outputs": file_records(
                [args.classification, args.genepred, args.corrected_gtf]
            ),
        }
    )
    write_json(args.run_json, record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
