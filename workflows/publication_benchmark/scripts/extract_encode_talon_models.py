#!/usr/bin/env python3
"""Extract sample-specific TALON models from ENCODE merged annotation GTFs."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import re
import time
from collections import Counter
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


TRANSCRIPT_ID = re.compile(r'(?:^|;\s*)transcript_id\s+"([^"]+)"')
GENE_ID = re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')
KEPT_FEATURES = {"transcript", "exon"}


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


@contextmanager
def _open_output(path: Path, *, compressed: bool) -> Iterator[TextIO]:
    if not compressed:
        with path.open("w", encoding="utf-8", newline="") as handle:
            yield handle
        return
    with path.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, mtime=0
        ) as compressed_stream:
            with io.TextIOWrapper(
                compressed_stream, encoding="utf-8", newline=""
            ) as handle:
                yield handle


def _digest(path: Path, algorithm: str = "sha256") -> tuple[int, str]:
    digest = hashlib.new(algorithm)
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return size, digest.hexdigest()


def _payload_digest(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with _open_text(path) as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), ""):
            encoded = block.encode("utf-8")
            size += len(encoded)
            digest.update(encoded)
    return size, digest.hexdigest()


def _file_record(path: Path) -> dict[str, object]:
    size, sha256 = _digest(path)
    return {
        "path": str(path.resolve()),
        "bytes": size,
        "sha256": sha256,
    }


def extract(
    source: Path, output: Path, *, selected_source: str = "TALON"
) -> dict[str, object]:
    """Select TALON transcripts and all exons belonging to their transcript IDs."""

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    compressed = output.suffix == ".gz"
    feature_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    transcript_ids: set[str] = set()
    gene_ids: set[str] = set()
    transcript_lines = 0
    exon_lines = 0
    try:
        with _open_text(source) as reader:
            for line_number, raw_line in enumerate(reader, 1):
                if raw_line.startswith("#") or not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"{source}:{line_number}: expected 9 GTF columns, "
                        f"found {len(fields)}"
                    )
                source_counts[fields[1]] += 1
                feature_counts[fields[2]] += 1
                if fields[1] != selected_source or fields[2] != "transcript":
                    continue
                transcript_match = TRANSCRIPT_ID.search(fields[8])
                gene_match = GENE_ID.search(fields[8])
                if transcript_match is None or gene_match is None:
                    raise ValueError(
                        f"{source}:{line_number}: selected row lacks gene_id or "
                        "transcript_id"
                    )
                transcript_id = transcript_match.group(1)
                if transcript_id in transcript_ids:
                    raise ValueError(
                        f"{source}:{line_number}: duplicate transcript row for "
                        f"{transcript_id}"
                    )
                transcript_ids.add(transcript_id)
                gene_ids.add(gene_match.group(1))
        if not transcript_ids:
            raise ValueError(f"{source}: no {selected_source} transcript rows found")

        exon_transcript_ids: set[str] = set()
        with _open_text(source) as reader, _open_output(
            temporary, compressed=compressed
        ) as writer:
            for line_number, raw_line in enumerate(reader, 1):
                if raw_line.startswith("#") or not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"{source}:{line_number}: expected 9 GTF columns, "
                        f"found {len(fields)}"
                    )
                if fields[2] not in KEPT_FEATURES:
                    continue
                transcript_match = TRANSCRIPT_ID.search(fields[8])
                if transcript_match is None:
                    continue
                transcript_id = transcript_match.group(1)
                if transcript_id not in transcript_ids:
                    continue
                if fields[2] == "transcript":
                    if fields[1] != selected_source:
                        raise ValueError(
                            f"{source}:{line_number}: selected transcript ID "
                            f"{transcript_id} has a non-{selected_source} transcript row"
                        )
                    transcript_lines += 1
                else:
                    exon_transcript_ids.add(transcript_id)
                    exon_lines += 1
                writer.write(raw_line if raw_line.endswith("\n") else raw_line + "\n")
        missing_exons = sorted(transcript_ids - exon_transcript_ids)
        if missing_exons:
            raise ValueError(
                f"{source}: {len(missing_exons)} transcripts lack exon rows"
            )
        os.replace(temporary, output)
    finally:
        if temporary.exists():
            temporary.unlink()

    payload_bytes, payload_sha256 = _payload_digest(output)
    output_record = _file_record(output)
    output_record.update(
        {
            "compressed": compressed,
            "payload_bytes": payload_bytes,
            "payload_sha256": payload_sha256,
            "gzip_header_policy": (
                "RFC 1952 with mtime=0 and no original filename"
                if compressed
                else None
            ),
        }
    )
    return {
        "selector": {
            "transcript_source_column": selected_source,
            "exon_rule": "transcript_id belongs to a selected transcript",
            "features": sorted(KEPT_FEATURES),
        },
        "selected_lines": transcript_lines + exon_lines,
        "selected_transcript_lines": transcript_lines,
        "selected_exon_lines": exon_lines,
        "selected_transcripts": len(transcript_ids),
        "selected_genes": len(gene_ids),
        "input_feature_counts": dict(sorted(feature_counts.items())),
        "input_source_counts": dict(sorted(source_counts.items())),
        "input": _file_record(source),
        "output": output_record,
    }


def _parse_input(value: str) -> tuple[str, Path, Path]:
    fields = value.split("=", 2)
    if len(fields) != 3 or not all(fields):
        raise ValueError(
            f"invalid --input {value!r}; expected dataset=SOURCE_GTF=OUTPUT_GTF"
        )
    return fields[0], Path(fields[1]), Path(fields[2])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        action="append",
        required=True,
        help="dataset=SOURCE_GTF=OUTPUT_GTF[.gz]",
    )
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--source-column", default="TALON")
    args = parser.parse_args()
    started = time.time()
    records = []
    for dataset, source, output in sorted(
        (_parse_input(value) for value in args.input), key=lambda row: row[0]
    ):
        record = extract(source, output, selected_source=args.source_column)
        records.append({"dataset": dataset, **record})
    audit = {
        "schema": "txid.encode-talon-extraction.v1",
        "selection_contract": (
            "Select transcript rows whose GTF source column is TALON, then retain "
            "every exon row carrying one of those selected transcript_id values, "
            "regardless of exon source."
        ),
        "started_unix": started,
        "ended_unix": time.time(),
        "datasets": records,
    }
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.audit.with_name(f".{args.audit.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, args.audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
