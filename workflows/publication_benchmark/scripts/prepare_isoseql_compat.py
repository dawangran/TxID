#!/usr/bin/env python3
"""Stage an audited isoSeQL source tree that retains duplicate observations."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import time
from pathlib import Path


SOURCE_FILES = ("isoSeQL_run.py", "isoSeQL_db.py", "isoSeQL_parse.py")
ORIGINAL = (
    "c.execute('INSERT INTO ends_counts(ends_id, exp, read_count) VALUES (?,?,?)', "
    "(isoEndID, expID, classif[iso].count))"
)
REPLACEMENT = (
    "c.execute('INSERT OR IGNORE INTO ends_counts(ends_id, exp, read_count) "
    "VALUES (?,?,?)', (isoEndID, expID, classif[iso].count))"
)


def _file_record(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return {
        "path": str(path.resolve()),
        "bytes": size,
        "sha256": digest.hexdigest(),
    }


def prepare(source: Path, output: Path) -> dict[str, object]:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite compatibility source: {output}")
    missing = [name for name in SOURCE_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(
            f"isoSeQL source lacks required files: {', '.join(missing)}"
        )
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    temporary.mkdir(parents=True)
    try:
        for name in SOURCE_FILES:
            shutil.copyfile(source / name, temporary / name)
        database_module = temporary / "isoSeQL_db.py"
        text = database_module.read_text(encoding="utf-8")
        if text.count(ORIGINAL) != 1:
            raise ValueError(
                "isoSeQL duplicate-end insertion contract was not found exactly once"
            )
        database_module.write_text(
            text.replace(ORIGINAL, REPLACEMENT),
            encoding="utf-8",
            newline="",
        )
        os.replace(temporary, output)
    except Exception:
        if temporary.exists():
            shutil.rmtree(temporary)
        raise
    return {
        "schema": "txid.isoseql-compat-source.v1",
        "upstream_version": "isoSeQL 1.0.1",
        "upstream_commit": "25d9366d8b236d3b912e62dcd2c267fe6df31a4c",
        "compatibility_change": {
            "id": "duplicate-exact-ends-observation-retention",
            "reason": (
                "Unmodified isoSeQL raises a UNIQUE constraint error when one "
                "experiment contains multiple transcript names for the same exact "
                "end structure and the same read-count value."
            ),
            "operation": (
                "Change only the ends_counts insertion to INSERT OR IGNORE; PBID "
                "insertion remains unchanged, so every source transcript name is "
                "retained for identity-partition evaluation."
            ),
            "structural_queries_changed": False,
        },
        "source_files": [_file_record(source / name) for name in SOURCE_FILES],
        "output_files": [_file_record(output / name) for name in SOURCE_FILES],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    record = prepare(args.source.resolve(), args.output.resolve())
    record["started_unix"] = started
    record["ended_unix"] = time.time()
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.audit.with_name(f".{args.audit.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, args.audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
