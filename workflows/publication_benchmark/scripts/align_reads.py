#!/usr/bin/env python3
"""Align one long-read dataset once and record a reusable sorted BAM."""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import tempfile
import time
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


def _preset(read_type: str, sirv: bool) -> list[str]:
    if read_type == "pacbio_ccs":
        args = ["-ax", "splice:hq", "-uf"]
    elif read_type in {"pacbio", "pacbio_clr"}:
        args = ["-ax", "splice:hq", "-uf"]
    elif read_type == "ont_directrna":
        args = ["-ax", "splice", "-uf", "-k14"]
    elif read_type in {"ont", "ont_cdna"}:
        args = ["-ax", "splice"]
    else:
        raise ValueError(f"unsupported read_type: {read_type}")
    args.extend(["--secondary=no", "--MD"])
    if sirv:
        args.append("--splice-flank=no")
    return args


def _version(command: list[str]) -> str:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    return (result.stdout or result.stderr).strip()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--reads", type=Path, required=True)
    parser.add_argument("--read-type", required=True)
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--bai", type=Path, required=True)
    parser.add_argument("--flagstat", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--sirv", action="store_true")
    args = parser.parse_args()
    args.bam.parent.mkdir(parents=True, exist_ok=True)
    minimap_argv = [
        "minimap2",
        *_preset(args.read_type, args.sirv),
        "-t",
        str(args.threads),
        str(args.reference),
        str(args.reads),
    ]
    samtools_argv = [
        "samtools",
        "sort",
        "-@",
        str(args.threads),
        "-o",
        str(args.bam),
    ]
    started_at = time.time()
    started = time.perf_counter()
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    with tempfile.NamedTemporaryFile(
        "wb", dir=args.bam.parent, prefix=".minimap2.", suffix=".sam", delete=False
    ) as temporary:
        temporary_path = Path(temporary.name)
        align = subprocess.run(
            minimap_argv,
            stdout=temporary,
            stderr=subprocess.PIPE,
            check=False,
        )
    (args.bam.parent / "minimap2.stderr.log").write_bytes(align.stderr)
    if align.returncode != 0:
        temporary_path.unlink(missing_ok=True)
        raise SystemExit(f"minimap2 failed with exit {align.returncode}")
    with temporary_path.open("rb") as sam_handle:
        sort = subprocess.run(
            samtools_argv,
            stdin=sam_handle,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    temporary_path.unlink(missing_ok=True)
    (args.bam.parent / "samtools-sort.stdout.log").write_bytes(sort.stdout)
    (args.bam.parent / "samtools-sort.stderr.log").write_bytes(sort.stderr)
    if sort.returncode != 0:
        raise SystemExit(f"samtools sort failed with exit {sort.returncode}")
    index_record = run_command(
        ["samtools", "index", "-@", str(args.threads), str(args.bam), str(args.bai)],
        log_prefix=args.bam.parent / "samtools-index",
    )
    flagstat_record = run_command(
        ["samtools", "flagstat", "-@", str(args.threads), str(args.bam)],
        log_prefix=args.bam.parent / "samtools-flagstat",
    )
    flagstat_stdout = Path(str(args.bam.parent / "samtools-flagstat") + ".stdout.log")
    args.flagstat.write_text(flagstat_stdout.read_text(encoding="utf-8"), encoding="utf-8")
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    record = {
        "stage": "alignment",
        "started_unix": started_at,
        "ended_unix": time.time(),
        "wall_seconds": time.perf_counter() - started,
        "child_user_seconds": after.ru_utime - before.ru_utime,
        "child_system_seconds": after.ru_stime - before.ru_stime,
        "max_rss_kib": after.ru_maxrss,
        "minimap2_version": _version(["minimap2", "--version"]),
        "samtools_version": _version(["samtools", "--version"]).splitlines()[0],
        "minimap2_argv": minimap_argv,
        "samtools_sort_argv": samtools_argv,
        "index": index_record,
        "flagstat": flagstat_record,
        "inputs": file_records([args.reference, args.reads]),
        "outputs": file_records([args.bam, args.bai, args.flagstat]),
        "environment": {
            key: os.environ[key]
            for key in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV")
            if key in os.environ
        },
    }
    write_json(args.run_json, record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
