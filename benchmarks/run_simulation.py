#!/usr/bin/env python3
"""Run the reproducible TxID synthetic benchmark from a source checkout."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

try:
    from txid.simulation import run_simulation
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from txid.simulation import run_simulation


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path("benchmarks/results"))
    parser.add_argument("--work-dir", type=Path, help="retain generated FASTA, annotations, and registry")
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--transcripts", type=int, default=120)
    parser.add_argument("--samples", type=int, default=6)
    parser.add_argument("--tools", type=int, default=4, choices=range(1, 5))
    parser.add_argument("--artifact-rate", type=float, default=0.05)
    args = parser.parse_args()
    summary = run_simulation(
        args.output_dir,
        seed=args.seed,
        transcript_count=args.transcripts,
        sample_count=args.samples,
        tool_count=args.tools,
        artifact_rate=args.artifact_rate,
        work_dir=args.work_dir,
    )
    json.dump(summary, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

