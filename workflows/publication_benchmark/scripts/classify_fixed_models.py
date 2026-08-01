#!/usr/bin/env python3
"""Classify fixed transcript models without persisting registry observations."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import os
import time
from collections import Counter
from pathlib import Path

from txid.classify import build_reference_index, classify
from txid.registry import Registry
from txid.workflow import _stage_annotation


COLUMNS = [
    "sample",
    "tool",
    "original_transcript_id",
    "original_gene_id",
    "txid_transcript_id",
    "txid_gene_id",
    "txid_sc",
    "txid_form",
    "classification",
    "annotation_name",
    "gene_candidates",
]


def classify_models(
    *,
    database: Path,
    annotation_name: str,
    input_path: Path,
    sample: str,
    tool: str,
    output: Path,
) -> dict[str, object]:
    """Run the production parser/identity/classifier without allocating loci."""

    started = time.perf_counter()
    with Registry(database, read_only=True) as registry:
        reference = registry.reference()
        annotation = registry.annotation(annotation_name)
        references = registry.references(int(annotation["id"]))
    reference_index = build_reference_index(references)
    staged = _stage_annotation(input_path, reference)
    counts: Counter[str] = Counter()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    if temporary.exists():
        raise FileExistsError(
            f"refusing to overwrite temporary classification: {temporary}"
        )
    try:
        with temporary.open("wb") as raw_handle:
            with gzip.GzipFile(
                filename="",
                mode="wb",
                fileobj=raw_handle,
                compresslevel=6,
                mtime=0,
            ) as compressed_handle:
                with io.TextIOWrapper(
                    compressed_handle, encoding="utf-8", newline=""
                ) as handle:
                    writer = csv.DictWriter(
                        handle,
                        fieldnames=COLUMNS,
                        delimiter="\t",
                        lineterminator="\n",
                    )
                    writer.writeheader()
                    for transcript, identities in staged:
                        result = classify(
                            transcript, identities, reference_index
                        )
                        counts[result.label] += 1
                        writer.writerow(
                            {
                                "sample": sample,
                                "tool": tool,
                                "original_transcript_id": (
                                    transcript.original_transcript_id
                                ),
                                "original_gene_id": (
                                    transcript.original_gene_id
                                ),
                                "txid_transcript_id": (
                                    result.output_transcript_id
                                    or identities.form.public_id
                                ),
                                # This route deliberately does not allocate
                                # registry-managed GL accessions.
                                "txid_gene_id": result.output_gene_id or "",
                                "txid_sc": (
                                    identities.splice_chain.public_id
                                    if identities.splice_chain is not None
                                    else ""
                                ),
                                "txid_form": identities.form.public_id,
                                "classification": result.label,
                                "annotation_name": annotation_name,
                                "gene_candidates": ",".join(
                                    result.gene_candidates
                                ),
                            }
                        )
            raw_handle.flush()
            os.fsync(raw_handle.fileno())
        os.replace(temporary, output)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    return {
        "annotation": annotation_name,
        "classification_counts": dict(sorted(counts.items())),
        "database": str(database),
        "input": str(input_path),
        "mode": "classification_only_no_observation_persistence",
        "output": str(output),
        "reference_transcripts": len(references),
        "sample": sample,
        "tool": tool,
        "transcripts": len(staged),
        "wall_seconds": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--annotation-name", required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--tool", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    record = classify_models(
        database=args.database,
        annotation_name=args.annotation_name,
        input_path=args.input,
        sample=args.sample,
        tool=args.tool,
        output=args.output,
    )
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    args.audit.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(record, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
