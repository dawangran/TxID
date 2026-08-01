#!/usr/bin/env python3
"""Create deterministic contig-filtered GTFs with a checksum audit."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
from collections import Counter
from pathlib import Path


def _digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def filter_gtf(source: Path, target: Path, contigs: set[str]) -> dict[str, object]:
    if not contigs:
        raise ValueError("at least one contig is required")
    if target.exists():
        raise FileExistsError(f"refusing to overwrite subset: {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.{os.getpid()}.tmp")
    features: Counter[str] = Counter()
    comments = 0
    records = 0
    try:
        with (
            _open_text(source) as reader,
            temporary.open("wb") as raw_output,
            gzip.GzipFile(
                filename="", mode="wb", fileobj=raw_output, mtime=0
            ) as compressed,
            io.TextIOWrapper(compressed, encoding="utf-8", newline="\n") as writer,
        ):
            for line_number, raw_line in enumerate(reader, 1):
                if raw_line.startswith("#"):
                    writer.write(raw_line)
                    comments += 1
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"{source}:{line_number}: expected 9 GTF columns, "
                        f"found {len(fields)}"
                    )
                if fields[0] not in contigs:
                    continue
                writer.write(raw_line)
                records += 1
                features[fields[2]] += 1
        os.replace(temporary, target)
    except BaseException:
        if temporary.exists():
            temporary.unlink()
        raise
    if not records:
        raise ValueError(f"{source}: selected contigs contain no records")
    return {
        "source": str(source.resolve()),
        "source_bytes": source.stat().st_size,
        "source_sha256": _digest(source),
        "output": str(target.resolve()),
        "output_bytes": target.stat().st_size,
        "output_sha256": _digest(target),
        "comment_lines": comments,
        "records": records,
        "features": dict(sorted(features.items())),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="label=GTF_OR_GTF_GZ",
    )
    parser.add_argument("--contig", action="append", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    labels = set()
    for value in args.input:
        label, separator, raw_path = value.partition("=")
        if not separator or not label or not raw_path:
            raise ValueError(f"invalid --input {value!r}; expected label=GTF")
        if label in labels:
            raise ValueError(f"duplicate input label: {label}")
        labels.add(label)
        source = Path(raw_path)
        target = args.output_dir / f"{label}.gtf.gz"
        rows.append(
            {
                "label": label,
                **filter_gtf(source, target, set(args.contig)),
            }
        )
    audit = {
        "schema": "txid.gtf-contig-subset.v1",
        "contigs": sorted(set(args.contig)),
        "files": sorted(rows, key=lambda row: row["label"]),
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.audit.with_name(f".{args.audit.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, args.audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
