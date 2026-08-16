"""Reference FASTA and explicit contig alias handling."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import TextIO

from .errors import ReferenceError
from .models import ContigRecord, SequenceCollection, TranscriptModel


def _open_text(path: str | Path) -> TextIO:
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline=None)
    return path.open("r", encoding="utf-8", newline=None)


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return "sha256:" + digest.hexdigest()


def load_aliases(path: str | Path | None, primary_names: set[str]) -> tuple[tuple[str, str], ...]:
    if path is None:
        return ()
    aliases: dict[str, str] = {}
    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\r\n")
            if not line or line.startswith("#"):
                continue
            columns = line.split("\t")
            if line_number == 1 and [item.lower() for item in columns[:2]] == ["alias", "primary"]:
                continue
            if len(columns) != 2 or not all(columns):
                raise ReferenceError(
                    f"alias table line {line_number}: expected non-empty alias and primary columns"
                )
            alias, primary = columns
            if primary not in primary_names:
                raise ReferenceError(
                    f"alias table line {line_number}: target {primary!r} is absent from FASTA"
                )
            if alias in primary_names and alias != primary:
                raise ReferenceError(
                    f"alias table line {line_number}: primary contig {alias!r} cannot alias {primary!r}"
                )
            previous = aliases.get(alias)
            if previous is not None and previous != primary:
                raise ReferenceError(
                    f"alias {alias!r} maps to both {previous!r} and {primary!r}"
                )
            aliases[alias] = primary
    return tuple(sorted(aliases.items()))


def read_fasta(
    path: str | Path, *, assembly_name: str, alias_path: str | Path | None = None
) -> SequenceCollection:
    contigs_by_name: dict[str, ContigRecord] = {}
    current: str | None = None
    current_length = 0
    current_digest = hashlib.sha256()
    digest_buffer: list[str] = []
    digest_buffer_length = 0

    def flush_digest_buffer() -> None:
        nonlocal digest_buffer_length
        if digest_buffer:
            current_digest.update("".join(digest_buffer).encode("ascii"))
            digest_buffer.clear()
            digest_buffer_length = 0

    def finish_contig() -> None:
        nonlocal current_length, current_digest
        if current is None:
            return
        flush_digest_buffer()
        contigs_by_name[current] = ContigRecord(
            current,
            current_length,
            current_digest.hexdigest(),
        )
        current_length = 0
        current_digest = hashlib.sha256()

    with _open_text(path) as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split(None, 1)[0]
                if not name:
                    raise ReferenceError(f"FASTA line {line_number}: empty sequence name")
                if name in contigs_by_name or name == current:
                    raise ReferenceError(f"FASTA line {line_number}: duplicate sequence name {name!r}")
                finish_contig()
                current = name
            elif current is None:
                raise ReferenceError(f"FASTA line {line_number}: sequence occurs before a header")
            else:
                sequence = "".join(line.split()).upper()
                if not sequence.isascii():
                    raise ReferenceError(f"FASTA line {line_number}: sequence must be ASCII")
                current_length += len(sequence)
                digest_buffer.append(sequence)
                digest_buffer_length += len(sequence)
                if digest_buffer_length >= 1024 * 1024:
                    flush_digest_buffer()
    finish_contig()
    if not contigs_by_name:
        raise ReferenceError("reference FASTA contains no sequences")

    contigs = [contigs_by_name[name] for name in sorted(contigs_by_name)]
    for contig in contigs:
        if contig.length == 0:
            raise ReferenceError(f"FASTA sequence {contig.name!r} is empty")
    fingerprint_rows = [
        f"{contig.name}\t{contig.length}\t{contig.sequence_digest}\n"
        for contig in contigs
    ]
    fingerprint = "sha256:" + hashlib.sha256(
        "".join(fingerprint_rows).encode("utf-8")
    ).hexdigest()
    aliases = load_aliases(alias_path, set(contigs_by_name))
    return SequenceCollection(
        name=assembly_name,
        fingerprint=fingerprint,
        contigs=tuple(contigs),
        aliases=aliases,
    )


def resolve_and_validate(
    transcript: TranscriptModel,
    reference: SequenceCollection,
    *,
    contig_map: dict[str, ContigRecord] | None = None,
    alias_map: dict[str, str] | None = None,
) -> TranscriptModel:
    if contig_map is None or alias_map is None:
        contig = reference.resolve(transcript.contig)
    else:
        primary = alias_map.get(transcript.contig)
        if primary is None:
            raise ReferenceError(
                f"contig {transcript.contig!r} is absent from reference {reference.name!r}; "
                "supply an explicit validated alias table at registry initialization"
            )
        contig = contig_map[primary]
    if transcript.end > contig.length:
        raise ReferenceError(
            f"transcript {transcript.original_transcript_id!r} ends at {transcript.end}, "
            f"beyond contig {contig.name!r} length {contig.length}"
        )
    return transcript if transcript.contig == contig.name else transcript.with_contig(contig.name)
