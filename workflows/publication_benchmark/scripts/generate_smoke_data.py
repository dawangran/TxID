#!/usr/bin/env python3
"""Generate deterministic sequence-backed truth and full-length smoke reads."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import random
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Model:
    transcript_id: str
    gene_id: str
    strand: str
    exons: tuple[tuple[int, int], ...]
    annotation_release: int


MODELS = (
    Model("truth.g1.main", "truth.g1", "+", ((1001, 1200), (1501, 1700), (2101, 2350)), 1),
    Model("truth.g1.ends", "truth.g1", "+", ((951, 1200), (1501, 1700), (2101, 2400)), 2),
    Model("truth.g1.alt", "truth.g1", "+", ((1001, 1200), (1601, 1700), (2101, 2350)), 2),
    Model("truth.g2.single", "truth.g2", "+", ((4001, 4650),), 1),
    Model("truth.g3.main", "truth.g3", "-", ((7001, 7200), (7601, 7800), (8201, 8500)), 1),
    Model("truth.g3.ends", "truth.g3", "-", ((6951, 7200), (7601, 7800), (8201, 8550)), 2),
    Model("truth.g4.main", "truth.g4", "+", ((10501, 10750), (11201, 11400)), 1),
    Model("truth.g5.main", "truth.g5", "-", ((14001, 14300), (14801, 15100)), 1),
)
SAMPLE_COUNTS = {
    "smoke-pb-r1": {
        model.transcript_id: 12 + (3 if model.transcript_id.endswith("main") else 0)
        for model in MODELS
    },
    "smoke-pb-r2": {
        model.transcript_id: (
            18
            if model.transcript_id in {"truth.g1.ends", "truth.g3.ends"}
            else 10
        )
        for model in MODELS
        if model.transcript_id != "truth.g1.alt"
    },
}


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _attributes(values: dict[str, str]) -> str:
    return " ".join(f'{key} "{value}";' for key, value in values.items())


def _gtf(models: tuple[Model, ...]) -> str:
    by_gene: dict[str, list[Model]] = {}
    for model in models:
        by_gene.setdefault(model.gene_id, []).append(model)
    lines = []
    for gene_id in sorted(by_gene):
        gene_models = by_gene[gene_id]
        strand = gene_models[0].strand
        start = min(model.exons[0][0] for model in gene_models)
        end = max(model.exons[-1][1] for model in gene_models)
        lines.append(
            "\t".join(
                [
                    "chrSmoke",
                    "txid-smoke",
                    "gene",
                    str(start),
                    str(end),
                    ".",
                    strand,
                    ".",
                    _attributes({"gene_id": gene_id, "gene_name": gene_id}),
                ]
            )
        )
        for model in sorted(gene_models, key=lambda item: item.transcript_id):
            base = {
                "gene_id": gene_id,
                "transcript_id": model.transcript_id,
                "truth_id": model.transcript_id,
            }
            lines.append(
                "\t".join(
                    [
                        "chrSmoke",
                        "txid-smoke",
                        "transcript",
                        str(model.exons[0][0]),
                        str(model.exons[-1][1]),
                        ".",
                        model.strand,
                        ".",
                        _attributes(base),
                    ]
                )
            )
            ordered = (
                model.exons if model.strand == "+" else tuple(reversed(model.exons))
            )
            exon_number = {exon: str(index) for index, exon in enumerate(ordered, 1)}
            for exon in model.exons:
                lines.append(
                    "\t".join(
                        [
                            "chrSmoke",
                            "txid-smoke",
                            "exon",
                            str(exon[0]),
                            str(exon[1]),
                            ".",
                            model.strand,
                            ".",
                            _attributes({**base, "exon_number": exon_number[exon]}),
                        ]
                    )
                )
    return "\n".join(lines) + "\n"


def _splice_motifs(sequence: list[str]) -> None:
    for model in MODELS:
        for left, right in zip(model.exons, model.exons[1:]):
            if model.strand == "+":
                sequence[left[1] : left[1] + 2] = "GT"
                sequence[right[0] - 3 : right[0] - 1] = "AG"
            else:
                sequence[left[1] : left[1] + 2] = "CT"
                sequence[right[0] - 3 : right[0] - 1] = "AC"


def _mutate(sequence: str, rng: random.Random, rate: float) -> str:
    bases = "ACGT"
    output = list(sequence)
    for index, base in enumerate(output):
        if rng.random() < rate:
            output[index] = rng.choice(bases.replace(base, ""))
    return "".join(output)


def _write_fastq(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw_handle:
        with gzip.GzipFile(filename="", fileobj=raw_handle, mode="wb", mtime=0) as handle:
            for read_id, sequence in records:
                text = f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                handle.write(text.encode("ascii"))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def generate(output_dir: Path, *, seed: int, reads_per_transcript: int, substitution_rate: float) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(seed)
    sequence = [rng.choice("ACGT") for _ in range(20000)]
    _splice_motifs(sequence)
    reference_sequence = "".join(sequence)
    fasta = output_dir / "reference.fa"
    fasta.write_text(
        ">chrSmoke\n"
        + "\n".join(
            reference_sequence[index : index + 80]
            for index in range(0, len(reference_sequence), 80)
        )
        + "\n",
        encoding="ascii",
    )
    (output_dir / "truth.gtf").write_text(_gtf(MODELS), encoding="utf-8")
    (output_dir / "annotation-v1.gtf").write_text(
        _gtf(tuple(model for model in MODELS if model.annotation_release <= 1)),
        encoding="utf-8",
    )
    (output_dir / "annotation-v2.gtf").write_text(_gtf(MODELS), encoding="utf-8")

    truth_rows = []
    for model in MODELS:
        introns = [
            [left[1] + 1, right[0] - 1]
            for left, right in zip(model.exons, model.exons[1:])
        ]
        if model.strand == "-":
            introns.reverse()
        truth_rows.append(
            {
                "truth_id": model.transcript_id,
                "gene_id": model.gene_id,
                "contig": "chrSmoke",
                "strand": model.strand,
                "exons": json.dumps(model.exons, separators=(",", ":")),
                "introns": json.dumps(introns, separators=(",", ":")),
                "annotation_release": model.annotation_release,
            }
        )
    with (output_dir / "truth.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(truth_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(truth_rows)

    read_truth_rows = []
    count_rows = []
    by_id = {model.transcript_id: model for model in MODELS}
    for sample, configured_counts in SAMPLE_COUNTS.items():
        records = []
        for truth_id, configured_count in sorted(configured_counts.items()):
            count = max(2, round(configured_count * reads_per_transcript / 12))
            model = by_id[truth_id]
            transcript = "".join(
                reference_sequence[start - 1 : end] for start, end in model.exons
            )
            if model.strand == "-":
                transcript = _reverse_complement(transcript)
            count_rows.append({"sample": sample, "truth_id": truth_id, "count": count})
            for number in range(1, count + 1):
                read_id = f"{sample}|{truth_id}|r{number:04d}"
                observed = _mutate(transcript, rng, substitution_rate)
                records.append((read_id, observed))
                read_truth_rows.append(
                    {"read_id": read_id, "sample": sample, "truth_id": truth_id}
                )
        _write_fastq(output_dir / f"{sample}.fastq.gz", records)

    for name, rows, fields in (
        ("read-truth.tsv", read_truth_rows, ["read_id", "sample", "truth_id"]),
        ("truth-counts.tsv", count_rows, ["sample", "truth_id", "count"]),
    ):
        with (output_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    generated = sorted(path for path in output_dir.iterdir() if path.is_file())
    (output_dir / "SHA256SUMS").write_text(
        "".join(f"{_sha256(path)}  {path.name}\n" for path in generated),
        encoding="utf-8",
    )
    metadata = {
        "generator": "generate_smoke_data.py",
        "seed": seed,
        "reads_per_transcript_baseline": reads_per_transcript,
        "substitution_rate": substitution_rate,
        "transcripts": len(MODELS),
        "samples": len(SAMPLE_COUNTS),
        "purpose": "integration smoke test; not biological publication evidence",
    }
    (output_dir / "generation.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--reads-per-transcript", type=int, default=12)
    parser.add_argument("--substitution-rate", type=float, default=0.003)
    args = parser.parse_args()
    generate(
        args.output_dir,
        seed=args.seed,
        reads_per_transcript=args.reads_per_transcript,
        substitution_rate=args.substitution_rate,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
