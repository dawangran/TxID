#!/usr/bin/env python3
"""Run one real transcript caller and normalize its primary GTF output path."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


def _find_gtf(root: Path, patterns: tuple[str, ...]) -> Path:
    candidates = []
    for pattern in patterns:
        candidates.extend(root.rglob(pattern))
    candidates = sorted(
        {path.resolve() for path in candidates if path.is_file() and path.stat().st_size},
        key=lambda path: (len(path.parts), str(path)),
    )
    if not candidates:
        raise FileNotFoundError(f"no non-empty caller GTF found below {root}")
    return candidates[0]


def _version(argv: list[str], environment: dict[str, str] | None = None) -> str:
    import os

    result = subprocess.run(
        argv,
        text=True,
        capture_output=True,
        check=False,
        env=dict(os.environ, **dict(environment or {})),
    )
    return (result.stdout or result.stderr).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caller", choices=("isoquant", "flair", "stringtie"), required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--read-type", required=True)
    parser.add_argument("--output-gtf", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_gtf.parent.mkdir(parents=True, exist_ok=True)

    if args.caller == "isoquant":
        isoquant_entry = Path(__file__).with_name("isoquant_entry.py")
        command = [
            sys.executable,
            str(isoquant_entry),
            "--reference",
            str(args.reference),
            "--genedb",
            str(args.annotation),
            "--complete_genedb",
            "--bam",
            str(args.bam),
            "--data_type",
            args.read_type if args.read_type in {"pacbio_ccs", "pacbio", "nanopore"} else "nanopore",
            "--output",
            str(args.work_dir),
            "--prefix",
            args.dataset,
            "--threads",
            str(args.threads),
            "--force",
        ]
        environment = {
            "NUMBA_CACHE_DIR": str((args.work_dir / "numba-cache").resolve()),
            "MPLCONFIGDIR": str((args.work_dir / "matplotlib-cache").resolve()),
            "XDG_CONFIG_HOME": str((args.work_dir / "xdg-config").resolve()),
            "XDG_CACHE_HOME": str((args.work_dir / "xdg-cache").resolve()),
            "TXID_ISOQUANT_CONFIG_DIR": str(
                (args.work_dir / "isoquant-config").resolve()
            ),
        }
        patterns = ("*.transcript_models.gtf", "*.transcript_models.gff3")
        version = _version(
            [sys.executable, str(isoquant_entry), "--version"], environment
        )
    elif args.caller == "flair":
        prefix = args.work_dir / args.dataset
        command = [
            "flair",
            "transcriptome",
            "-b",
            str(args.bam),
            "-g",
            str(args.reference),
            "-f",
            str(args.annotation),
            "-o",
            str(prefix),
            "--threads",
            str(args.threads),
        ]
        environment = {
            "MPLCONFIGDIR": str((args.work_dir / "matplotlib-cache").resolve()),
            "XDG_CONFIG_HOME": str((args.work_dir / "xdg-config").resolve()),
            "XDG_CACHE_HOME": str((args.work_dir / "xdg-cache").resolve()),
        }
        patterns = ("*.isoforms.gtf", "*.gtf")
        version = _version(["flair", "--version"])
    else:
        direct_output = args.work_dir / f"{args.dataset}.gtf"
        command = [
            "stringtie",
            "-L",
            "-G",
            str(args.annotation),
            "-p",
            str(args.threads),
            "-o",
            str(direct_output),
            str(args.bam),
        ]
        environment = {}
        patterns = (f"{args.dataset}.gtf", "*.gtf")
        version = _version(["stringtie", "--version"])

    record = run_command(
        command,
        log_prefix=args.work_dir / args.caller,
        env=environment,
    )
    discovered = _find_gtf(args.work_dir, patterns)
    if discovered != args.output_gtf.resolve():
        shutil.copyfile(discovered, args.output_gtf)
    if args.output_gtf.stat().st_size == 0:
        raise RuntimeError(f"{args.caller} emitted an empty GTF")
    record.update(
        {
            "stage": "transcript_discovery",
            "caller": args.caller,
            "dataset": args.dataset,
            "version": version,
            "compatibility_entry": (
                "isoquant_entry.py redirects IsoQuant 3.13.0 mutable config "
                "from its hard-coded home path; analysis code is unchanged"
                if args.caller == "isoquant"
                else None
            ),
            "normalized_from": str(discovered),
            "inputs": file_records([args.reference, args.annotation, args.bam]),
            "outputs": file_records([args.output_gtf]),
        }
    )
    write_json(args.run_json, record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
