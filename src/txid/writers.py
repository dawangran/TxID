"""Deterministic, atomic TxID output writers."""

from __future__ import annotations

import csv
import json
import os
import tempfile
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import Mapping, TextIO

from .errors import TxIDError
from .models import Assignment, Attributes

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
_GTF_ATTRIBUTE_RANK = {
    "gene_id": 0,
    "transcript_id": 1,
    "exon_number": 2,
}
_EXON_NUMBER = "exon_number"


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
    if "\\" not in value and '"' not in value and "\n" not in value:
        return value
    return value.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n")


def format_gtf_attributes(
    attributes: Mapping[str, str] | Iterable[tuple[str, str]],
) -> str:
    items = list(attributes.items() if isinstance(attributes, Mapping) else attributes)
    preferred: list[list[tuple[str, str]]] = [[], [], []]
    remaining: list[tuple[str, str]] = []
    for item in items:
        rank = _GTF_ATTRIBUTE_RANK.get(item[0])
        if rank is None:
            remaining.append(item)
        else:
            preferred[rank].append(item)
    remaining.sort(key=lambda item: item[0])
    ordered = preferred[0] + preferred[1] + preferred[2] + remaining
    return " ".join(
        [f'{key} "{_gtf_escape(str(value))}";' for key, value in ordered]
    )


def _upstream_attributes(
    values: Mapping[str, str] | Iterable[tuple[str, str]], reserved: set[str]
) -> Attributes:
    items = values.items() if isinstance(values, Mapping) else values
    result: list[tuple[str, str]] = []
    output_keys: dict[str, str] = {}
    owners: dict[str, str] = {}
    for key, value in items:
        if key in {"gene_id", "transcript_id"}:
            continue
        output_key = output_keys.get(key)
        if output_key is None:
            output_key = (
                f"txid_upstream_{key}"
                if key in reserved or key in {"ID", "Parent"}
                else key
            )
            while output_key in reserved or output_key in owners:
                output_key = "txid_upstream_" + output_key
            output_keys[key] = output_key
            owners[output_key] = key
        result.append((output_key, value))
    return tuple(result)


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
            base_items = tuple(base.items())
            base_keys = set(base)
            transcript_attrs = base_items + _upstream_attributes(
                transcript.attributes, base_keys
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
            exon_keys = base_keys | {_EXON_NUMBER}
            exon_count = transcript.exon_count
            for genomic_index, exon in enumerate(transcript.exons):
                exon_number = (
                    genomic_index + 1
                    if transcript.strand == "+"
                    else exon_count - genomic_index
                )
                exon_attrs = base_items + (
                    (_EXON_NUMBER, str(exon_number)),
                ) + _upstream_attributes(
                    exon.attributes, exon_keys
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
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(MAPPING_COLUMNS)
        for assignment in ordered:
            writer.writerow(
                (
                    sample,
                    tool,
                    assignment.transcript.original_transcript_id,
                    assignment.transcript.original_gene_id or "",
                    assignment.output_transcript_id,
                    assignment.output_gene_id,
                    assignment.identities.splice_chain.public_id
                    if assignment.identities.splice_chain
                    else "",
                    assignment.identities.form.public_id,
                    assignment.classification,
                    assignment.annotation_name,
                    ",".join(assignment.gene_candidates),
                    assignment.fuzzy_cluster or "",
                    assignment.fuzzy_bridge_status or "",
                )
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
