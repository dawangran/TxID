#!/usr/bin/env python3
"""Build SIRV transcript truth and assign public reads by unique full-length mapping."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import subprocess
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


ATTRIBUTE = re.compile(r'(?:^|;\s*)([^";\s]+)\s+"([^"]*)"')
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class Transcript:
    truth_id: str
    contig: str
    strand: str
    exons: tuple[tuple[int, int], ...]

    @property
    def introns(self) -> tuple[tuple[int, int], ...]:
        return tuple(
            (left[1] + 1, right[0] - 1)
            for left, right in zip(self.exons, self.exons[1:])
        )


@dataclass(frozen=True)
class Candidate:
    target: str
    score: int
    identity: float
    query_coverage: float
    target_coverage: float
    mapq: int


def _open_text(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def parse_truth_gtf(path: Path) -> list[Transcript]:
    records: dict[str, dict[str, object]] = {}
    with _open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}:{number}: expected 9 GTF fields")
            if fields[2].lower() != "exon":
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            transcript_id = attributes.get("transcript_id")
            if not transcript_id:
                raise ValueError(f"{path}:{number}: exon lacks transcript_id")
            contig, strand = fields[0], fields[6]
            if strand not in {"+", "-"}:
                raise ValueError(f"{path}:{number}: unsupported strand {strand!r}")
            start, end = int(fields[3]), int(fields[4])
            if start < 1 or end < start:
                raise ValueError(f"{path}:{number}: invalid exon interval")
            record = records.setdefault(
                transcript_id,
                {"contig": contig, "strand": strand, "exons": []},
            )
            if record["contig"] != contig or record["strand"] != strand:
                raise ValueError(f"{path}:{number}: inconsistent transcript context")
            exons = record["exons"]
            assert isinstance(exons, list)
            exons.append((start, end))
    if not records:
        raise ValueError(f"{path}: no exon records")
    transcripts = []
    for truth_id, record in records.items():
        exons = tuple(sorted(set(record["exons"])))
        for left, right in zip(exons, exons[1:]):
            if left[1] >= right[0]:
                raise ValueError(f"{path}: overlapping exons in {truth_id}")
        transcripts.append(
            Transcript(
                truth_id=truth_id,
                contig=str(record["contig"]),
                strand=str(record["strand"]),
                exons=exons,
            )
        )
    return sorted(transcripts, key=lambda value: value.truth_id)


def _load_selected_contigs(reference: Path, names: set[str]) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    current: str | None = None
    with _open_text(reference) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].strip().split()[0]
                current = name if name in names else None
                if current is not None:
                    if current in sequences:
                        raise ValueError(f"{reference}: duplicate FASTA contig {current}")
                    sequences[current] = []
            elif current is not None:
                sequences[current].append(line.strip())
    missing = sorted(names - set(sequences))
    if missing:
        raise ValueError(
            f"{reference}: missing truth contigs: {', '.join(missing[:5])}"
        )
    return {name: "".join(parts).upper() for name, parts in sequences.items()}


def build_transcriptome(
    reference: Path,
    truth_gtf: Path,
    transcript_fasta: Path,
    transcript_truth: Path,
) -> list[Transcript]:
    transcripts = parse_truth_gtf(truth_gtf)
    contigs = _load_selected_contigs(
        reference, {transcript.contig for transcript in transcripts}
    )
    transcript_fasta.parent.mkdir(parents=True, exist_ok=True)
    with transcript_fasta.open("w", encoding="utf-8", newline="\n") as fasta:
        for transcript in transcripts:
            contig = contigs[transcript.contig]
            sequence = "".join(
                contig[start - 1 : end] for start, end in transcript.exons
            )
            if transcript.strand == "-":
                sequence = sequence.translate(COMPLEMENT)[::-1]
            fasta.write(f">{transcript.truth_id}\n")
            for offset in range(0, len(sequence), 80):
                fasta.write(sequence[offset : offset + 80] + "\n")
    transcript_truth.parent.mkdir(parents=True, exist_ok=True)
    with transcript_truth.open("w", newline="", encoding="utf-8") as handle:
        fields = ["truth_id", "contig", "strand", "exons", "introns", "length"]
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for transcript in transcripts:
            writer.writerow(
                {
                    "truth_id": transcript.truth_id,
                    "contig": transcript.contig,
                    "strand": transcript.strand,
                    "exons": json.dumps(transcript.exons, separators=(",", ":")),
                    "introns": json.dumps(
                        transcript.introns, separators=(",", ":")
                    ),
                    "length": sum(end - start + 1 for start, end in transcript.exons),
                }
            )
    return transcripts


def _parse_tags(values: list[str]) -> dict[str, str]:
    tags = {}
    for value in values:
        parts = value.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def parse_paf(path: Path) -> dict[str, list[Candidate]]:
    by_query_target: dict[tuple[str, str], Candidate] = {}
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{path}:{number}: malformed PAF")
            query, target = fields[0], fields[5]
            query_length, query_start, query_end = map(
                int, (fields[1], fields[2], fields[3])
            )
            target_length, target_start, target_end = map(
                int, (fields[6], fields[7], fields[8])
            )
            matches, alignment_length, mapq = map(
                int, (fields[9], fields[10], fields[11])
            )
            if query_length <= 0 or target_length <= 0 or alignment_length <= 0:
                raise ValueError(f"{path}:{number}: non-positive PAF length")
            tags = _parse_tags(fields[12:])
            score = int(tags.get("AS", matches))
            candidate = Candidate(
                target=target,
                score=score,
                identity=matches / alignment_length,
                query_coverage=(query_end - query_start) / query_length,
                target_coverage=(target_end - target_start) / target_length,
                mapq=mapq,
            )
            key = (query, target)
            previous = by_query_target.get(key)
            if previous is None or (
                candidate.score,
                candidate.identity,
                candidate.query_coverage,
                candidate.target_coverage,
            ) > (
                previous.score,
                previous.identity,
                previous.query_coverage,
                previous.target_coverage,
            ):
                by_query_target[key] = candidate
    by_query: dict[str, list[Candidate]] = defaultdict(list)
    for (query, _), candidate in by_query_target.items():
        by_query[query].append(candidate)
    return dict(by_query)


def assign_candidates(
    by_query: dict[str, list[Candidate]],
    *,
    min_identity: float,
    min_query_coverage: float,
    min_target_coverage: float,
    min_score_margin: int,
) -> tuple[list[tuple[str, str]], Counter[str]]:
    accepted = []
    counts: Counter[str] = Counter()
    for query in sorted(by_query):
        candidates = sorted(
            by_query[query],
            key=lambda value: (
                value.score,
                value.identity,
                value.query_coverage,
                value.target_coverage,
                value.mapq,
                value.target,
            ),
            reverse=True,
        )
        best = candidates[0]
        if best.identity < min_identity:
            counts["below_identity"] += 1
            continue
        if best.query_coverage < min_query_coverage:
            counts["below_query_coverage"] += 1
            continue
        if best.target_coverage < min_target_coverage:
            counts["below_target_coverage"] += 1
            continue
        if len(candidates) > 1 and best.score - candidates[1].score < min_score_margin:
            counts["ambiguous_score_margin"] += 1
            continue
        accepted.append((query, best.target))
        counts["accepted"] += 1
    counts["paf_queries"] = len(by_query)
    return accepted, counts


def _input_spec(value: str) -> tuple[str, str, Path]:
    parts = value.split("=", 2)
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "--input must be DATASET_ID=READ_TYPE=FASTQ"
        )
    return parts[0], parts[1], Path(parts[2])


def run(args: argparse.Namespace) -> dict[str, object]:
    started = time.time()
    transcripts = build_transcriptome(
        args.reference,
        args.truth_gtf,
        args.transcript_fasta,
        args.transcript_truth,
    )
    args.work_dir.mkdir(parents=True, exist_ok=True)
    truth_rows = []
    datasets = []
    for dataset, read_type, fastq in sorted(args.input):
        preset = "map-hifi" if read_type == "pacbio_ccs" else "map-ont"
        paf = args.work_dir / f"{dataset}.sirv-transcriptome.paf"
        stderr = args.work_dir / f"{dataset}.minimap2.stderr.log"
        command = [
            args.minimap2,
            "-t",
            str(args.threads),
            "-x",
            preset,
            "--secondary=yes",
            "-N",
            "20",
            str(args.transcript_fasta),
            str(fastq),
        ]
        command_started = time.time()
        with paf.open("w", encoding="utf-8") as stdout, stderr.open(
            "w", encoding="utf-8"
        ) as error:
            result = subprocess.run(
                command, stdout=stdout, stderr=error, check=False, text=True
            )
        if result.returncode != 0:
            raise RuntimeError(
                f"minimap2 failed for {dataset} with exit {result.returncode}; "
                f"see {stderr}"
            )
        min_identity = (
            args.pacbio_min_identity
            if read_type == "pacbio_ccs"
            else args.ont_min_identity
        )
        assignments, counts = assign_candidates(
            parse_paf(paf),
            min_identity=min_identity,
            min_query_coverage=args.min_query_coverage,
            min_target_coverage=args.min_target_coverage,
            min_score_margin=args.min_score_margin,
        )
        truth_rows.extend(
            {"read_id": read_id, "sample": dataset, "truth_id": truth_id}
            for read_id, truth_id in assignments
        )
        datasets.append(
            {
                "dataset_id": dataset,
                "read_type": read_type,
                "fastq": str(fastq.resolve()),
                "fastq_bytes": fastq.stat().st_size,
                "paf": str(paf.resolve()),
                "paf_bytes": paf.stat().st_size,
                "paf_sha256": _sha256(paf),
                "stderr": str(stderr.resolve()),
                "argv": command,
                "returncode": result.returncode,
                "wall_seconds": time.time() - command_started,
                "thresholds": {
                    "min_identity": min_identity,
                    "min_query_coverage": args.min_query_coverage,
                    "min_target_coverage": args.min_target_coverage,
                    "min_score_margin": args.min_score_margin,
                },
                "assignment_counts": dict(sorted(counts.items())),
            }
        )
    truth_rows.sort(key=lambda row: (row["sample"], row["read_id"]))
    args.read_truth.parent.mkdir(parents=True, exist_ok=True)
    with args.read_truth.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["read_id", "sample", "truth_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(truth_rows)
    record = {
        "schema": "txid.sirv-read-truth.v1",
        "method": (
            "unique best minimap2 alignment to GTF-derived SIRV transcript "
            "sequences after identity, query-coverage, target-coverage, and "
            "alignment-score-margin filtering"
        ),
        "reference": str(args.reference.resolve()),
        "truth_gtf": str(args.truth_gtf.resolve()),
        "truth_transcripts": len(transcripts),
        "accepted_read_truth_rows": len(truth_rows),
        "transcript_fasta": {
            "path": str(args.transcript_fasta.resolve()),
            "sha256": _sha256(args.transcript_fasta),
        },
        "transcript_truth": {
            "path": str(args.transcript_truth.resolve()),
            "sha256": _sha256(args.transcript_truth),
        },
        "read_truth": {
            "path": str(args.read_truth.resolve()),
            "sha256": _sha256(args.read_truth),
        },
        "datasets": datasets,
        "started_unix": started,
        "ended_unix": time.time(),
    }
    _write_json(args.run_json, record)
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--truth-gtf", type=Path, required=True)
    parser.add_argument("--input", type=_input_spec, action="append", required=True)
    parser.add_argument("--transcript-fasta", type=Path, required=True)
    parser.add_argument("--transcript-truth", type=Path, required=True)
    parser.add_argument("--read-truth", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--pacbio-min-identity", type=float, default=0.97)
    parser.add_argument("--ont-min-identity", type=float, default=0.85)
    parser.add_argument("--min-query-coverage", type=float, default=0.90)
    parser.add_argument("--min-target-coverage", type=float, default=0.90)
    parser.add_argument("--min-score-margin", type=int, default=50)
    parser.add_argument("--minimap2", default="minimap2")
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
