#!/usr/bin/env python3
"""Fetch one external source repository at an exact revision."""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

from benchmark_lib import run_command, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository", required=True)
    parser.add_argument(
        "--cache",
        type=Path,
        help=(
            "optional local Git repository used for retrieval while retaining "
            "--repository as the declared upstream origin"
        ),
    )
    parser.add_argument("--revision", required=True)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--record", type=Path, required=True)
    args = parser.parse_args()
    if args.target.exists():
        raise FileExistsError(f"refusing to overwrite source checkout: {args.target}")
    retrieval_source = str(args.cache.resolve()) if args.cache else args.repository
    if args.cache is not None and not args.cache.exists():
        raise FileNotFoundError(f"source cache does not exist: {args.cache}")
    args.target.parent.mkdir(parents=True, exist_ok=True)
    staging = args.target.with_name(f".{args.target.name}.incomplete-{os.getpid()}")
    if staging.exists():
        raise FileExistsError(f"source staging path already exists: {staging}")
    commands = [
        run_command(
            ["git", "clone", "--no-checkout", retrieval_source, staging],
            log_prefix=args.record.parent / f"{args.target.name}-clone",
        ),
        run_command(
            ["git", "-C", staging, "checkout", "--detach", args.revision],
            log_prefix=args.record.parent / f"{args.target.name}-checkout",
        ),
    ]
    commit = subprocess.run(
        ["git", "-C", staging, "rev-parse", "HEAD"],
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    if len(args.revision) == 40 and commit != args.revision:
        raise RuntimeError(f"checked out {commit}, expected {args.revision}")
    os.replace(staging, args.target)
    write_json(
        args.record,
        {
            "repository": args.repository,
            "retrieval_source": retrieval_source,
            "used_local_cache": args.cache is not None,
            "requested_revision": args.revision,
            "resolved_commit": commit,
            "commands": commands,
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
