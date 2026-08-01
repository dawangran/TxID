"""Deterministic, atomic TxID output writers."""

from __future__ import annotations

import csv
import json
import os
import tempfile
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import TextIO

from .errors import TxIDError
from .models import Assignment

MAPPING_COLUMNS = [
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
    "fuzzy_cluster",
    "fuzzy_bridge_status",
]


def _atomic_text(path: str | Path, writer: Callable[[TextIO], None]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            "w", encoding="utf-8", newline="", dir=target.parent, prefix=f".{target.name}.", delete=False
        ) as handle:
            temporary_name = handle.name
            writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, target)
    except Exception:
        if temporary_name is not None:
            temporary = Path(temporary_name)
            if temporary.exists():
                temporary.unlink()
        raise


def ensure_distinct_paths(input_path: str | Path, outputs: Iterable[str | Path]) -> None:
    source = Path(input_path).resolve()
    resolved = [Path(path).resolve() for path in outputs]
    if source in resolved:
        raise TxIDError("refusing to overwrite the input annotation in place")
    if len(set(resolved)) != len(resolved):
        raise TxIDError("output annotation and mapping table must use different paths")


def _gtf_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")


def format_gtf_attributes(attributes: dict[str, str]) -> str:
    preferred = ["gene_id", "transcript_id", "exon_number"]
    keys = [key for key in preferred if key in attributes]
    keys.extend(sorted(set(attributes) - set(keys)))
    return " ".join(f'{key} "{_gtf_escape(str(attributes[key]))}";' for key in keys)


def _upstream_attributes(values: dict[str, str], reserved: set[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    for key, value in values.items():
        if key in {"gene_id", "transcript_id"}:
            continue
        output_key = f"txid_upstream_{key}" if key in reserved or key in {"ID", "Parent"} else key
        while output_key in result or output_key in reserved:
            output_key = "txid_upstream_" + output_key
        result[output_key] = value
    return result


def _base_attributes(assignment: Assignment) -> dict[str, str]:
    identity = assignment.identities
    attrs = {
        "gene_id": assignment.output_gene_id,
        "transcript_id": assignment.output_transcript_id,
        "txid_annotation": assignment.annotation_name,
        "txid_classification": assignment.classification,
        "txid_form": identity.form.public_id,
        "txid_original_transcript_id": assignment.transcript.original_transcript_id,
    }
    if assignment.transcript.original_gene_id is not None:
        attrs["txid_original_gene_id"] = assignment.transcript.original_gene_id
    if identity.splice_chain is not None:
        attrs["txid_sc"] = identity.splice_chain.public_id
    if assignment.gene_candidates:
        attrs["txid_gene_candidates"] = ",".join(assignment.gene_candidates)
    if assignment.fuzzy_cluster is not None:
        attrs["txid_fuzzy_cluster"] = assignment.fuzzy_cluster
    if assignment.fuzzy_bridge_status is not None:
        attrs["txid_fuzzy_bridge_status"] = assignment.fuzzy_bridge_status
    return attrs


def write_gtf(path: str | Path, assignments: Iterable[Assignment]) -> None:
    ordered = sorted(
        assignments,
        key=lambda item: (
            item.transcript.contig,
            item.transcript.start,
            item.transcript.end,
            item.transcript.strand,
            item.identities.form.public_id,
            item.transcript.original_transcript_id,
        ),
    )

    def emit(handle: TextIO) -> None:
        handle.write("# txid rewritten annotation; coordinates are 1-based closed\n")
        for assignment in ordered:
            transcript = assignment.transcript
            base = _base_attributes(assignment)
            transcript_attrs = dict(base)
            transcript_attrs.update(
                _upstream_attributes(transcript.attribute_dict(), set(base))
            )
            handle.write(
                "\t".join(
                    [
                        transcript.contig,
                        "txid",
                        "transcript",
                        str(transcript.start),
                        str(transcript.end),
                        ".",
                        transcript.strand,
                        ".",
                        format_gtf_attributes(transcript_attrs),
                    ]
                )
                + "\n"
            )
            transcript_order = (
                list(transcript.exons)
                if transcript.strand == "+"
                else list(reversed(transcript.exons))
            )
            exon_number = {exon: str(index) for index, exon in enumerate(transcript_order, 1)}
            for exon in transcript.exons:
                exon_attrs = dict(base)
                exon_attrs["exon_number"] = exon_number[exon]
                exon_attrs.update(
                    _upstream_attributes(dict(exon.attributes), set(exon_attrs))
                )
                handle.write(
                    "\t".join(
                        [
                            transcript.contig,
                            exon.source if exon.source != "." else "txid",
                            "exon",
                            str(exon.start),
                            str(exon.end),
                            exon.score,
                            transcript.strand,
                            exon.phase,
                            format_gtf_attributes(exon_attrs),
                        ]
                    )
                    + "\n"
                )

    _atomic_text(path, emit)


def write_mapping(
    path: str | Path,
    assignments: Iterable[Assignment],
    *,
    sample: str,
    tool: str,
) -> None:
    ordered = sorted(assignments, key=lambda item: item.transcript.original_transcript_id)

    def emit(handle: TextIO) -> None:
        writer = csv.DictWriter(handle, fieldnames=MAPPING_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for assignment in ordered:
            writer.writerow(
                {
                    "sample": sample,
                    "tool": tool,
                    "original_transcript_id": assignment.transcript.original_transcript_id,
                    "original_gene_id": assignment.transcript.original_gene_id or "",
                    "txid_transcript_id": assignment.output_transcript_id,
                    "txid_gene_id": assignment.output_gene_id,
                    "txid_sc": assignment.identities.splice_chain.public_id if assignment.identities.splice_chain else "",
                    "txid_form": assignment.identities.form.public_id,
                    "classification": assignment.classification,
                    "annotation_name": assignment.annotation_name,
                    "gene_candidates": ",".join(assignment.gene_candidates),
                    "fuzzy_cluster": assignment.fuzzy_cluster or "",
                    "fuzzy_bridge_status": assignment.fuzzy_bridge_status or "",
                }
            )

    _atomic_text(path, emit)


def write_catalog(path: str | Path, rows: Iterable[dict[str, object]]) -> None:
    columns = [
        "form_id",
        "splice_chain_id",
        "gene_id",
        "classifications",
        "observation_count",
        "sample_count",
        "tool_count",
    ]

    def emit(handle: TextIO) -> None:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    _atomic_text(path, emit)


def write_json(path: str | Path, value: object) -> None:
    def emit(handle: TextIO) -> None:
        json.dump(value, handle, sort_keys=True, indent=2, ensure_ascii=False)
        handle.write("\n")

    _atomic_text(path, emit)
