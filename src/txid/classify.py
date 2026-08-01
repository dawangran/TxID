"""Annotation-relative classification, kept separate from exact identity."""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from typing import Iterable

from .models import IdentityBundle, TranscriptModel


@dataclass(frozen=True, slots=True)
class ReferenceEntry:
    gene_id: str | None
    transcript_id: str
    contig: str
    strand: str
    start: int
    end: int
    form_id: str
    splice_chain_id: str | None


@dataclass(frozen=True, slots=True)
class Classification:
    label: str
    output_gene_id: str | None
    output_transcript_id: str
    gene_candidates: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class _GeneSpanIndex:
    starts: tuple[int, ...]
    spans: tuple[tuple[int, int, str], ...]
    prefix_max_ends: tuple[int, ...]

    @classmethod
    def build(
        cls, spans: Iterable[tuple[int, int, str]]
    ) -> "_GeneSpanIndex":
        ordered = tuple(sorted(spans, key=lambda item: (item[0], item[1], item[2])))
        maximum = 0
        prefix_max_ends = []
        for _, end, _ in ordered:
            maximum = max(maximum, end)
            prefix_max_ends.append(maximum)
        return cls(
            starts=tuple(start for start, _, _ in ordered),
            spans=ordered,
            prefix_max_ends=tuple(prefix_max_ends),
        )

    def overlapping(self, start: int, end: int) -> tuple[str, ...]:
        genes: set[str] = set()
        index = bisect_right(self.starts, end) - 1
        while index >= 0:
            if self.prefix_max_ends[index] < start:
                break
            _, span_end, gene_id = self.spans[index]
            if span_end >= start:
                genes.add(gene_id)
            index -= 1
        return tuple(sorted(genes))


@dataclass(frozen=True, slots=True)
class ReferenceIndex:
    """Reusable exact-form and aggregate-gene-span lookup index."""

    by_form: dict[str, tuple[ReferenceEntry, ...]]
    by_contig_strand: dict[tuple[str, str], _GeneSpanIndex]

    @classmethod
    def build(cls, references: Iterable[ReferenceEntry]) -> "ReferenceIndex":
        form_rows: dict[str, list[ReferenceEntry]] = {}
        gene_spans: dict[tuple[str, str, str], tuple[int, int]] = {}
        for reference in references:
            form_rows.setdefault(reference.form_id, []).append(reference)
            if reference.gene_id is None:
                continue
            key = (reference.gene_id, reference.contig, reference.strand)
            previous = gene_spans.get(key)
            if previous is None:
                gene_spans[key] = (reference.start, reference.end)
            else:
                gene_spans[key] = (
                    min(previous[0], reference.start),
                    max(previous[1], reference.end),
                )
        grouped_spans: dict[tuple[str, str], list[tuple[int, int, str]]] = {}
        for (gene_id, contig, strand), (start, end) in gene_spans.items():
            grouped_spans.setdefault((contig, strand), []).append(
                (start, end, gene_id)
            )
        return cls(
            by_form={
                form_id: tuple(
                    sorted(
                        rows,
                        key=lambda item: (
                            item.transcript_id,
                            item.gene_id or "",
                            item.contig,
                            item.strand,
                            item.start,
                            item.end,
                        ),
                    )
                )
                for form_id, rows in form_rows.items()
            },
            by_contig_strand={
                key: _GeneSpanIndex.build(spans)
                for key, spans in grouped_spans.items()
            },
        )

    def exact(self, form_id: str) -> tuple[ReferenceEntry, ...]:
        return self.by_form.get(form_id, ())

    def overlapping_genes(
        self, contig: str, strand: str, start: int, end: int
    ) -> tuple[str, ...]:
        index = self.by_contig_strand.get((contig, strand))
        return () if index is None else index.overlapping(start, end)


def build_reference_index(
    references: Iterable[ReferenceEntry] | ReferenceIndex,
) -> ReferenceIndex:
    if isinstance(references, ReferenceIndex):
        return references
    return ReferenceIndex.build(references)


def classify(
    transcript: TranscriptModel,
    identities: IdentityBundle,
    references: Iterable[ReferenceEntry] | ReferenceIndex,
) -> Classification:
    reference_index = build_reference_index(references)
    exact = reference_index.exact(identities.form.public_id)
    if exact:
        selected = min(
            exact,
            key=lambda item: (
                item.transcript_id != transcript.original_transcript_id,
                item.transcript_id,
                item.gene_id or "",
            ),
        )
        return Classification(
            label="known",
            output_gene_id=selected.gene_id,
            output_transcript_id=selected.transcript_id,
            gene_candidates=tuple(
                sorted({item.gene_id for item in exact if item.gene_id is not None})
            ),
        )

    overlapping_genes = reference_index.overlapping_genes(
        transcript.contig,
        transcript.strand,
        transcript.start,
        transcript.end,
    )
    if len(overlapping_genes) == 1:
        return Classification(
            label="novel_in_known_gene",
            output_gene_id=overlapping_genes[0],
            output_transcript_id=identities.form.public_id,
            gene_candidates=tuple(overlapping_genes),
        )
    if len(overlapping_genes) > 1:
        return Classification(
            label="ambiguous_gene",
            output_gene_id=None,
            output_transcript_id=identities.form.public_id,
            gene_candidates=tuple(overlapping_genes),
        )
    return Classification(
        label="new_locus",
        output_gene_id=None,
        output_transcript_id=identities.form.public_id,
    )
