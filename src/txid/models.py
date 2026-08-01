"""Shared, format-independent transcript data structures."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Iterable, Mapping

from .errors import AnnotationParseError

Attributes = tuple[tuple[str, str], ...]


def freeze_attributes(values: Mapping[str, str] | Iterable[tuple[str, str]]) -> Attributes:
    """Return attributes in stable key order while retaining the last duplicate."""
    return tuple(sorted(dict(values).items()))


@dataclass(frozen=True, slots=True)
class Exon:
    start: int
    end: int
    source: str = "."
    score: str = "."
    phase: str = "."
    attributes: Attributes = ()

    def __post_init__(self) -> None:
        if self.start < 1 or self.end < self.start:
            raise AnnotationParseError(
                f"invalid 1-based closed exon interval {self.start}-{self.end}"
            )


@dataclass(frozen=True, slots=True)
class TranscriptModel:
    contig: str
    strand: str
    exons: tuple[Exon, ...]
    original_transcript_id: str
    original_gene_id: str | None = None
    source: str = "."
    source_format: str = "gtf"
    attributes: Attributes = ()

    def __post_init__(self) -> None:
        if not self.contig:
            raise AnnotationParseError(
                f"transcript {self.original_transcript_id!r} has an empty contig"
            )
        if self.strand not in {"+", "-"}:
            raise AnnotationParseError(
                f"transcript {self.original_transcript_id!r} has unsupported strand {self.strand!r}"
            )
        if not self.original_transcript_id:
            raise AnnotationParseError("transcript has no identifier")
        if not self.exons:
            raise AnnotationParseError(
                f"transcript {self.original_transcript_id!r} has no exons"
            )
        ordered = tuple(sorted(self.exons, key=lambda exon: (exon.start, exon.end)))
        for left, right in zip(ordered, ordered[1:]):
            if left.end >= right.start:
                raise AnnotationParseError(
                    f"transcript {self.original_transcript_id!r} has overlapping exons "
                    f"{left.start}-{left.end} and {right.start}-{right.end}"
                )
        if ordered != self.exons:
            object.__setattr__(self, "exons", ordered)

    @property
    def start(self) -> int:
        return self.exons[0].start

    @property
    def end(self) -> int:
        return self.exons[-1].end

    @property
    def exon_count(self) -> int:
        return len(self.exons)

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end

    @property
    def tes(self) -> int:
        return self.end if self.strand == "+" else self.start

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        genomic = tuple(
            (left.end + 1, right.start - 1)
            for left, right in zip(self.exons, self.exons[1:])
        )
        return genomic if self.strand == "+" else tuple(reversed(genomic))

    def with_contig(self, primary_contig: str) -> "TranscriptModel":
        return replace(self, contig=primary_contig)

    def attribute_dict(self) -> dict[str, str]:
        return dict(self.attributes)


@dataclass(frozen=True, slots=True)
class ContigRecord:
    name: str
    length: int
    sequence_digest: str


@dataclass(frozen=True, slots=True)
class SequenceCollection:
    name: str
    fingerprint: str
    contigs: tuple[ContigRecord, ...]
    aliases: tuple[tuple[str, str], ...] = ()

    def contig_map(self) -> dict[str, ContigRecord]:
        return {contig.name: contig for contig in self.contigs}

    def alias_map(self) -> dict[str, str]:
        result = {contig.name: contig.name for contig in self.contigs}
        result.update(self.aliases)
        return result

    def resolve(self, name: str) -> ContigRecord:
        from .errors import ReferenceError

        primary = self.alias_map().get(name)
        if primary is None:
            raise ReferenceError(
                f"contig {name!r} is absent from reference {self.name!r}; "
                "supply an explicit validated alias table at registry initialization"
            )
        return self.contig_map()[primary]


@dataclass(frozen=True, slots=True)
class StructuralIdentity:
    family: str
    public_id: str
    public_digest: str
    full_digest: str
    canonical_json: str


@dataclass(frozen=True, slots=True)
class IdentityBundle:
    form: StructuralIdentity
    splice_chain: StructuralIdentity | None = None


@dataclass(frozen=True, slots=True)
class Assignment:
    transcript: TranscriptModel
    identities: IdentityBundle
    classification: str
    output_gene_id: str
    output_transcript_id: str
    annotation_name: str
    gene_candidates: tuple[str, ...] = ()
    fuzzy_cluster: str | None = None
    fuzzy_bridge_status: str | None = None
