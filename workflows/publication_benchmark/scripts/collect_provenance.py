#!/usr/bin/env python3
"""Collect rule execution records into one compact publication provenance bundle."""

from __future__ import annotations

import argparse
import csv
import json
import platform
import sys
from pathlib import Path

from benchmark_lib import file_records, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="append", type=Path, default=[])
    parser.add_argument("--manifest", action="append", type=Path, default=[])
    parser.add_argument("--artifact", action="append", type=Path, default=[])
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-tsv", type=Path, required=True)
    args = parser.parse_args()
    records = []
    rows = []
    for path in sorted(set(args.record), key=str):
        value = json.loads(path.read_text(encoding="utf-8"))
        records.append({"path": str(path.resolve()), "record": value})
        rows.append(
            {
                "record": str(path),
                "stage": value.get("stage", "unknown"),
                "tool": value.get("caller", value.get("tool", "")),
                "dataset": value.get("dataset", ""),
                "order": value.get("order", ""),
                "wall_seconds": value.get("wall_seconds", value.get("commands", [{}])[-1].get("wall_seconds", "") if value.get("commands") else ""),
                "max_rss_kib": value.get("max_rss_kib", value.get("commands", [{}])[-1].get("max_rss_kib", "") if value.get("commands") else ""),
            }
        )
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="", encoding="utf-8") as handle:
        fields = [
            "record",
            "stage",
            "tool",
            "dataset",
            "order",
            "wall_seconds",
            "max_rss_kib",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    write_json(
        args.output_json,
        {
            "schema": "txid-publication-provenance-1",
            "python": sys.version,
            "platform": platform.platform(),
            "manifests": file_records(args.manifest),
            "evaluated_artifacts": file_records(args.artifact),
            "records": records,
            "limitations": [
                "Tier-0 smoke data are generated integration fixtures, not biological evidence.",
                "Peak RSS fields are comparable only for records with identical stage boundaries.",
            ],
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
