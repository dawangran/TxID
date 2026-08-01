#!/usr/bin/env python3
"""Prepare one truth-defined input bundle for external comparator workflows."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

try:
    from txid.identity import identify
    from txid.models import TranscriptModel
    from txid.parser import parse_annotation
    from txid.reference import read_fasta
    from txid.simulation import run_simulation
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from txid.identity import identify
    from txid.models import TranscriptModel
    from txid.parser import parse_annotation
    from txid.reference import read_fasta
    from txid.simulation import run_simulation


SEED = 20260730


def _attributes(values: dict[str, str]) -> str:
    return " ".join(
        f'{key} "{value.replace(chr(34), chr(92) + chr(34))}";'
        for key, value in values.items()
    )


def _truth_id(model: TranscriptModel) -> str | None:
    value = model.attribute_dict().get("truth_id")
    if value:
        return value
    values = {
        dict(exon.attributes).get("truth_id")
        for exon in model.exons
        if dict(exon.attributes).get("truth_id")
    }
    if len(values) > 1:
        raise ValueError(f"conflicting truth_id values on {model.original_transcript_id}")
    return next(iter(values), None)


def _write_full_gtf(path: Path, models: list[TranscriptModel]) -> None:
    genes: dict[str, list[TranscriptModel]] = defaultdict(list)
    for model in models:
        genes[model.original_gene_id or f"unresolved-{model.original_transcript_id}"].append(model)
    lines: list[str] = []
    for gene_id in sorted(genes):
        transcripts = genes[gene_id]
        contigs = {model.contig for model in transcripts}
        strands = {model.strand for model in transcripts}
        if len(contigs) != 1 or len(strands) != 1:
            raise ValueError(f"gene {gene_id!r} spans multiple contigs or strands")
        gene_start = min(model.start for model in transcripts)
        gene_end = max(model.end for model in transcripts)
        contig = transcripts[0].contig
        strand = transcripts[0].strand
        gene_attrs = _attributes(
            {"gene_id": gene_id, "gene_name": gene_id, "gene_status": "KNOWN"}
        )
        lines.append(
            "\t".join(
                [contig, "simulation", "gene", str(gene_start), str(gene_end), ".", strand, ".", gene_attrs]
            )
        )
        for model in sorted(transcripts, key=lambda item: item.original_transcript_id):
            truth_id = _truth_id(model) or model.original_transcript_id
            tx_attrs = _attributes(
                {
                    "gene_id": gene_id,
                    "transcript_id": model.original_transcript_id,
                    "gene_name": gene_id,
                    "transcript_name": model.original_transcript_id,
                    "transcript_status": "KNOWN",
                    "truth_id": truth_id,
                }
            )
            lines.append(
                "\t".join(
                    [
                        model.contig,
                        model.source,
                        "transcript",
                        str(model.start),
                        str(model.end),
                        ".",
                        model.strand,
                        ".",
                        tx_attrs,
                    ]
                )
            )
            exon_order = model.exons if model.strand == "+" else tuple(reversed(model.exons))
            exon_number = {exon: number for number, exon in enumerate(exon_order, 1)}
            for exon in model.exons:
                exon_attrs = _attributes(
                    {
                        "gene_id": gene_id,
                        "transcript_id": model.original_transcript_id,
                        "exon_number": str(exon_number[exon]),
                        "truth_id": truth_id,
                    }
                )
                lines.append(
                    "\t".join(
                        [
                            model.contig,
                            model.source,
                            "exon",
                            str(exon.start),
                            str(exon.end),
                            ".",
                            model.strand,
                            ".",
                            exon_attrs,
                        ]
                    )
                )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _structure_key(model: TranscriptModel) -> str:
    return json.dumps(
        [model.contig, model.strand, [[exon.start, exon.end] for exon in model.exons]],
        separators=(",", ":"),
    )


def _intron_key(model: TranscriptModel) -> tuple[str, str, tuple[tuple[int, int], ...]]:
    return model.contig, model.strand, model.introns


def _sqanti_category(
    model: TranscriptModel,
    reference_by_gene: dict[str, list[TranscriptModel]],
) -> tuple[str, str, str]:
    gene_id = model.original_gene_id or f"unresolved-{model.original_transcript_id}"
    references = reference_by_gene.get(gene_id, [])
    exact = next((item for item in references if _structure_key(item) == _structure_key(model)), None)
    splice = next((item for item in references if _intron_key(item) == _intron_key(model)), None)
    if exact is not None or (model.exon_count > 1 and splice is not None):
        matched = exact or splice
        assert matched is not None
        return "full-splice_match", matched.original_transcript_id, str(matched.exon_count)
    if references:
        return "novel_not_in_catalog", "novel", "NA"
    truth_id = _truth_id(model) or ""
    truth_index = int(truth_id.split(":")[-1]) if truth_id.startswith("truth:") else 0
    return ("fusion" if truth_index % 10 == 0 else "intergenic"), "novel", "NA"


def _write_isoseql_inputs(
    root: Path,
    dataset: str,
    sample: str,
    tool: str,
    models: list[TranscriptModel],
    reference_by_gene: dict[str, list[TranscriptModel]],
) -> None:
    target = root / dataset
    target.mkdir(parents=True, exist_ok=True)
    class_header = [
        "isoform",
        "structural_category",
        "associated_gene",
        "associated_transcript",
        "ref_exons",
        "exons",
        "subcategory",
        "all_canonical",
        "FL",
        "length",
    ]
    class_rows: list[list[object]] = []
    gene_pred: list[str] = []
    for model in sorted(models, key=lambda item: item.original_transcript_id):
        gene_id = model.original_gene_id or f"unresolved-{model.original_transcript_id}"
        category, matched_tx, ref_exons = _sqanti_category(model, reference_by_gene)
        class_rows.append(
            [
                model.original_transcript_id,
                category,
                gene_id,
                matched_tx,
                ref_exons,
                model.exon_count,
                "mono-exon" if model.exon_count == 1 else "multi-exon",
                "canonical",
                1,
                sum(exon.end - exon.start + 1 for exon in model.exons),
            ]
        )
        exon_starts = ",".join(str(exon.start - 1) for exon in model.exons) + ","
        exon_ends = ",".join(str(exon.end) for exon in model.exons) + ","
        gene_pred.append(
            "\t".join(
                [
                    model.original_transcript_id,
                    model.contig,
                    model.strand,
                    str(model.start - 1),
                    str(model.end),
                    str(model.start - 1),
                    str(model.end),
                    str(model.exon_count),
                    exon_starts,
                    exon_ends,
                ]
            )
        )
    with (target / "classification.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(class_header)
        writer.writerows(class_rows)
    (target / "models.genePred").write_text("\n".join(gene_pred) + "\n", encoding="utf-8")
    with (target / "sample.config").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["sample_name", "tissue", "disease", "age", "sex"])
        writer.writerow([sample, "synthetic", "none", 0, "NA"])
    with (target / "experiment.config").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(
            [
                "RIN",
                "date",
                "platform",
                "method",
                "vMap",
                "vReference",
                "vAnnot",
                "vLima",
                "vCCS",
                "vIsoseq3",
                "vCupcake",
                "vSQANTI",
                "exp_name",
            ]
        )
        writer.writerow(
            [10, "2026-07-31", "synthetic", tool, "exact-SAM", "TxID-simulation-v1", "ref-v1", "NA", "NA", "NA", "NA", "compatible-fixture", dataset]
        )


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _write_sam(path: Path, models: list[TranscriptModel], sequence: str) -> None:
    lines = ["@HD\tVN:1.6\tSO:coordinate", f"@SQ\tSN:chrSim\tLN:{len(sequence)}"]
    for model in sorted(models, key=lambda item: (item.start, item.end, item.original_transcript_id)):
        cigar_parts: list[str] = []
        for index, exon in enumerate(model.exons):
            if index:
                previous = model.exons[index - 1]
                cigar_parts.append(f"{exon.start - previous.end - 1}N")
            cigar_parts.append(f"{exon.end - exon.start + 1}M")
        aligned = "".join(sequence[exon.start - 1 : exon.end] for exon in model.exons)
        if model.strand == "-":
            aligned = _reverse_complement(aligned)
        length = len(aligned)
        fields = [
            model.original_transcript_id,
            "0" if model.strand == "+" else "16",
            model.contig,
            str(model.start),
            "60",
            "".join(cigar_parts),
            "*",
            "0",
            "0",
            aligned,
            "I" * length,
            f"MD:Z:{length}",
            "NM:i:0",
            f"AS:i:{length}",
            "fA:f:0.0",
        ]
        lines.append("\t".join(fields))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare(output_dir: Path, work_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty output directory: {output_dir}")
    if work_dir.exists() and any(work_dir.iterdir()):
        raise FileExistsError(f"refusing to overwrite non-empty work directory: {work_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    run_simulation(
        output_dir / "txid-baseline",
        seed=SEED,
        transcript_count=120,
        sample_count=6,
        tool_count=4,
        artifact_rate=0.05,
        work_dir=work_dir,
    )
    fasta = work_dir / "simulation.fa"
    reference_models = parse_annotation(work_dir / "reference-v1.gtf", fmt="gtf")
    reference_gtf = output_dir / "reference.gtf"
    _write_full_gtf(reference_gtf, reference_models)
    (output_dir / "simulation.fa").write_bytes(fasta.read_bytes())
    sequence = "".join(
        line.strip().upper()
        for line in fasta.read_text(encoding="utf-8").splitlines()
        if not line.startswith(">")
    )
    reference = read_fasta(fasta, assembly_name="TxID-simulation-v1")
    reference_by_gene: dict[str, list[TranscriptModel]] = defaultdict(list)
    for model in reference_models:
        if model.original_gene_id:
            reference_by_gene[model.original_gene_id].append(model)

    raw_paths = sorted(
        [
            path
            for path in [*work_dir.glob("sample-*.gtf"), *work_dir.glob("sample-*.gff3")]
            if not path.name.endswith(".txid.gtf")
        ],
        key=lambda path: path.name,
    )
    query_dir = output_dir / "queries"
    sam_dir = output_dir / "sam"
    isoseql_dir = output_dir / "isoseql"
    manifest_rows: list[dict[str, object]] = []
    truth_rows: list[dict[str, object]] = []
    talon_config: list[str] = []
    for path in raw_paths:
        sample, tool, _suffix = path.name.split(".", 2)
        dataset = f"{sample}__{tool}"
        models = parse_annotation(path)
        normalized_gtf = query_dir / f"{dataset}.gtf"
        _write_full_gtf(normalized_gtf, models)
        sam_path = sam_dir / f"{dataset}.sam"
        _write_sam(sam_path, models, sequence)
        _write_isoseql_inputs(isoseql_dir, dataset, sample, tool, models, reference_by_gene)
        talon_config.append(f"{dataset},synthetic,{tool},{sam_path.resolve()}")
        manifest_rows.append(
            {
                "dataset": dataset,
                "sample": sample,
                "tool": tool,
                "query_gtf": str(normalized_gtf.resolve()),
                "sam": str(sam_path.resolve()),
                "isoseql_dir": str((isoseql_dir / dataset).resolve()),
                "transcripts": len(models),
            }
        )
        for model in models:
            truth_id = _truth_id(model)
            if not truth_id:
                raise ValueError(f"missing truth_id on {model.original_transcript_id}")
            identity = identify(model, reference.fingerprint)
            truth_rows.append(
                {
                    "dataset": dataset,
                    "sample": sample,
                    "tool": tool,
                    "original_transcript_id": model.original_transcript_id,
                    "truth_id": truth_id,
                    "structural_truth": _structure_key(model),
                    "splice_truth": json.dumps(_intron_key(model), separators=(",", ":")),
                    "txid_form": identity.form.public_id,
                }
            )

    with (output_dir / "inputs-manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    with (output_dir / "truth-manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(truth_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(truth_rows, key=lambda row: (str(row["dataset"]), str(row["original_transcript_id"]))))
    (output_dir / "talon-config.csv").write_text("\n".join(talon_config) + "\n", encoding="utf-8")
    shuffled = talon_config.copy()
    random.Random(SEED + 1).shuffle(shuffled)
    (output_dir / "talon-config-shuffled.csv").write_text("\n".join(shuffled) + "\n", encoding="utf-8")

    checksummed = sorted(
        path for path in output_dir.rglob("*") if path.is_file() and path.name != "SHA256SUMS"
    )
    (output_dir / "SHA256SUMS").write_text(
        "".join(f"{_sha256(path)}  {path.relative_to(output_dir)}\n" for path in checksummed),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.output_dir, args.work_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
