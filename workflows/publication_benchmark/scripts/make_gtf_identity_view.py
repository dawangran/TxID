#!/usr/bin/env python3
"""Build an audited structure-and-identifier-only view of a GTF."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import os
import sys
import time
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


KEPT_FEATURES = {"transcript", "exon"}
IDENTIFIER_KEYS = ("gene_id", "transcript_id")


def _load_txid(repo: Path):
    source = str((repo / "src").resolve())
    if source not in sys.path:
        sys.path.insert(0, source)
    from txid.parser import parse_gtf_attribute_items
    from txid.writers import format_gtf_attributes

    return parse_gtf_attribute_items, format_gtf_attributes


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


def _digest(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
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


def project(source: Path, output: Path, *, repo: Path) -> dict[str, object]:
    parse_items, format_attributes = _load_txid(repo)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    compressed = output.suffix == ".gz"
    duplicate_keys: Counter[str] = Counter()
    input_features: Counter[str] = Counter()
    kept_features: Counter[str] = Counter()
    transcript_ids: set[str] = set()
    exon_transcript_ids: set[str] = set()
    try:
        with _open_text(source) as reader, _open_output(
            temporary, compressed=compressed
        ) as writer:
            writer.write(
                "# TxID benchmark identity view: coordinates and original "
                "gene_id/transcript_id only\n"
            )
            for line_number, raw_line in enumerate(reader, 1):
                if raw_line.startswith("#") or not raw_line.strip():
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(
                        f"{source}:{line_number}: expected 9 GTF columns, "
                        f"found {len(fields)}"
                    )
                feature = fields[2].lower()
                input_features[feature] += 1
                if feature not in KEPT_FEATURES:
                    continue
                values: dict[str, list[str]] = defaultdict(list)
                for key, value in parse_items(fields[8], line_number):
                    if values[key]:
                        duplicate_keys[key] += 1
                    values[key].append(value)
                identifiers = {}
                for key in IDENTIFIER_KEYS:
                    observed = values.get(key, [])
                    unique = sorted(set(observed))
                    if not observed:
                        raise ValueError(
                            f"{source}:{line_number}: {feature} row lacks {key}"
                        )
                    if len(unique) != 1:
                        raise ValueError(
                            f"{source}:{line_number}: conflicting duplicate {key} "
                            f"values: {unique}"
                        )
                    identifiers[key] = unique[0]
                transcript_id = identifiers["transcript_id"]
                if feature == "transcript":
                    if transcript_id in transcript_ids:
                        raise ValueError(
                            f"{source}:{line_number}: duplicate transcript row for "
                            f"{transcript_id}"
                        )
                    transcript_ids.add(transcript_id)
                else:
                    exon_transcript_ids.add(transcript_id)
                fields[8] = format_attributes(identifiers)
                writer.write("\t".join(fields) + "\n")
                kept_features[feature] += 1
        if not exon_transcript_ids:
            raise ValueError(f"{source}: identity view contains no exon rows")
        transcript_rows_missing_exons = transcript_ids - exon_transcript_ids
        if transcript_rows_missing_exons:
            raise ValueError(
                f"{source}: {len(transcript_rows_missing_exons)} transcript rows "
                "lack exon rows"
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
        "projection_contract": {
            "kept_features": sorted(KEPT_FEATURES),
            "kept_attributes": list(IDENTIFIER_KEYS),
            "coordinate_columns": "preserved byte-for-byte",
            "record_order": "preserved",
            "duplicate_nonidentifier_attributes": "counted in audit and omitted",
            "duplicate_identifier_attributes": (
                "accepted only when all values are identical"
            ),
        },
        "input_feature_counts": dict(sorted(input_features.items())),
        "kept_feature_counts": dict(sorted(kept_features.items())),
        "duplicate_attribute_occurrences": sum(duplicate_keys.values()),
        "duplicate_attribute_occurrences_by_key": dict(sorted(duplicate_keys.items())),
        "transcript_rows": len(transcript_ids),
        "transcripts_with_exons": len(exon_transcript_ids),
        "input": _file_record(source),
        "output": output_record,
    }


def _parse_input(value: str) -> tuple[str, Path, Path]:
    fields = value.split("=", 2)
    if len(fields) != 3 or not all(fields):
        raise ValueError(
            f"invalid --input {value!r}; expected label=SOURCE_GTF=OUTPUT_GTF[.gz]"
        )
    return fields[0], Path(fields[1]), Path(fields[2])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    records = []
    for label, source, output in sorted(
        (_parse_input(value) for value in args.input), key=lambda row: row[0]
    ):
        records.append(
            {
                "label": label,
                **project(source, output, repo=args.repo.resolve()),
            }
        )
    audit = {
        "schema": "txid.gtf-identity-view.v1",
        "started_unix": started,
        "ended_unix": time.time(),
        "views": records,
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
