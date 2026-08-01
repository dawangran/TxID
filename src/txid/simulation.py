"""Reproducible truth-defined simulation and interoperability evaluation."""

from __future__ import annotations

import csv
import json
import math
import platform
import random
import sqlite3
import tempfile
import time
import tracemalloc
from collections import Counter, defaultdict
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

from . import __version__
from .identity import identify
from .models import Exon, TranscriptModel
from .registry import Registry
from .workflow import add_annotation_context, import_file, initialize_registry
from .writers import _atomic_text, write_json


@dataclass(frozen=True, slots=True)
class TruthTranscript:
    truth_id: str
    category: str
    gene_id: str | None
    model: TranscriptModel


@dataclass(frozen=True, slots=True)
class ObservationTruth:
    sample: str
    tool: str
    original_transcript_id: str
    truth_id: str
    structural_truth: str
    expected_classification: str
    count: int
    perturbed: bool


def _model(
    index: int,
    exons: Iterable[tuple[int, int]],
    *,
    strand: str,
    transcript_id: str,
    gene_id: str | None,
) -> TranscriptModel:
    return TranscriptModel(
        "chrSim",
        strand,
        tuple(Exon(start, end, source="simulation") for start, end in exons),
        transcript_id,
        gene_id,
        source="simulation",
    )


def generate_truth(transcript_count: int) -> tuple[list[TruthTranscript], list[TranscriptModel], list[TranscriptModel]]:
    if transcript_count < 10:
        raise ValueError("simulation requires at least 10 truth transcripts")
    truth: list[TruthTranscript] = []
    reference_v1: list[TranscriptModel] = []
    reference_v2: list[TranscriptModel] = []
    for index in range(transcript_count):
        base = 1000 + index * 3000
        strand = "+" if index % 2 == 0 else "-"
        truth_id = f"truth:{index + 1:05d}"
        group = index % 10
        if group <= 4:
            category = "known"
            gene_id = f"G{index + 1:05d}"
            model = _model(
                index,
                [(base, base + 99), (base + 300, base + 399)],
                strand=strand,
                transcript_id=truth_id,
                gene_id=gene_id,
            )
            reference = replace(
                model,
                original_transcript_id=f"REF{index + 1:05d}",
                source="reference",
            )
            reference_v1.append(reference)
            reference_v2.append(reference)
        elif group in {5, 6}:
            category = "novel_in_known_gene"
            gene_id = f"G{index + 1:05d}"
            model = _model(
                index,
                [(base + 10, base + 99), (base + 300, base + 410)],
                strand=strand,
                transcript_id=truth_id,
                gene_id=gene_id,
            )
            surrogate = _model(
                index,
                [(base, base + 99), (base + 300, base + 399)],
                strand=strand,
                transcript_id=f"REF{index + 1:05d}",
                gene_id=gene_id,
            )
            reference_v1.append(surrogate)
            reference_v2.extend(
                [
                    surrogate,
                    replace(
                        model,
                        original_transcript_id=f"PROMOTED{index + 1:05d}",
                        source="reference-v2",
                    ),
                ]
            )
        elif group == 7:
            category = "new_locus"
            gene_id = None
            exons = (
                [(base + 50, base + 180)]
                if index % 20 == 7
                else [(base, base + 99), (base + 300, base + 399)]
            )
            model = _model(
                index,
                exons,
                strand=strand,
                transcript_id=truth_id,
                gene_id=None,
            )
        elif group == 8:
            category = "known"
            gene_id = f"G{index + 1:05d}"
            model = _model(
                index,
                [(base + 50, base + 180)],
                strand=strand,
                transcript_id=truth_id,
                gene_id=gene_id,
            )
            reference = replace(
                model,
                original_transcript_id=f"REF{index + 1:05d}",
                source="reference",
            )
            reference_v1.append(reference)
            reference_v2.append(reference)
        else:
            category = "ambiguous_gene"
            gene_id = None
            model = _model(
                index,
                [(base + 150, base + 210), (base + 400, base + 460)],
                strand=strand,
                transcript_id=truth_id,
                gene_id=None,
            )
            first = _model(
                index,
                [(base, base + 80), (base + 260, base + 320)],
                strand=strand,
                transcript_id=f"REF{index + 1:05d}A",
                gene_id=f"G{index + 1:05d}A",
            )
            second = _model(
                index,
                [(base + 250, base + 330), (base + 550, base + 620)],
                strand=strand,
                transcript_id=f"REF{index + 1:05d}B",
                gene_id=f"G{index + 1:05d}B",
            )
            reference_v1.extend([first, second])
            reference_v2.extend([first, second])
        truth.append(TruthTranscript(truth_id, category, gene_id, model))
    return truth, reference_v1, reference_v2


def _gtf(models: Iterable[TranscriptModel], rng: random.Random | None = None) -> str:
    lines: list[str] = []
    for model in models:
        exon_lines = []
        for exon in model.exons:
            attributes = [
                f'gene_id "{model.original_gene_id or "unresolved"}";',
                f'transcript_id "{model.original_transcript_id}";',
                f'truth_id "{model.attribute_dict().get("truth_id", model.original_transcript_id)}";',
            ]
            exon_lines.append(
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
                        " ".join(attributes),
                    ]
                )
            )
        if rng is not None:
            rng.shuffle(exon_lines)
        lines.extend(exon_lines)
    if rng is not None:
        rng.shuffle(lines)
    return "\n".join(lines) + "\n"


def _gff3(models: Iterable[TranscriptModel], rng: random.Random) -> str:
    blocks: list[list[str]] = []
    for model in models:
        gene = model.original_gene_id or f"unresolved-{model.original_transcript_id}"
        block = [
            "\t".join(
                [
                    model.contig,
                    model.source,
                    "mRNA",
                    str(model.start),
                    str(model.end),
                    ".",
                    model.strand,
                    ".",
                    f"ID={model.original_transcript_id};Parent={gene};truth_id={model.attribute_dict().get('truth_id', model.original_transcript_id)}",
                ]
            )
        ]
        exon_lines = [
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
                    f"ID={model.original_transcript_id}.e{number};Parent={model.original_transcript_id}",
                ]
            )
            for number, exon in enumerate(model.exons, 1)
        ]
        rng.shuffle(exon_lines)
        block.extend(exon_lines)
        blocks.append(block)
    rng.shuffle(blocks)
    return "##gff-version 3\n" + "\n".join(line for block in blocks for line in block) + "\n"


def _perturb(model: TranscriptModel) -> TranscriptModel:
    if model.exon_count == 1:
        first = model.exons[0]
        exons = (replace(first, start=first.start + 1),)
    else:
        first, *rest = model.exons
        exons = (replace(first, end=first.end + 1), *rest)
    return replace(model, exons=tuple(exons))


def _structural_key(model: TranscriptModel) -> str:
    return json.dumps(
        [model.contig, model.strand, [[e.start, e.end] for e in model.exons]],
        separators=(",", ":"),
    )


def _write_annotation(path: Path, models: Iterable[TranscriptModel]) -> None:
    path.write_text(_gtf(models), encoding="utf-8", newline="\n")


def _pairwise_error_rates(expected: list[str], predicted: list[str]) -> dict[str, float | int]:
    if len(expected) != len(predicted):
        raise ValueError("expected and predicted arrays differ in length")
    expected_groups: dict[str, Counter[str]] = defaultdict(Counter)
    predicted_groups: dict[str, Counter[str]] = defaultdict(Counter)
    for truth, prediction in zip(expected, predicted):
        expected_groups[truth][prediction] += 1
        predicted_groups[prediction][truth] += 1
    choose2 = lambda value: value * (value - 1) // 2
    same_truth_pairs = sum(choose2(sum(group.values())) for group in expected_groups.values())
    correctly_joined = sum(
        choose2(count) for group in expected_groups.values() for count in group.values()
    )
    same_prediction_pairs = sum(choose2(sum(group.values())) for group in predicted_groups.values())
    correctly_unmerged = sum(
        choose2(count) for group in predicted_groups.values() for count in group.values()
    )
    split_pairs = same_truth_pairs - correctly_joined
    merge_pairs = same_prediction_pairs - correctly_unmerged
    return {
        "false_split_pairs": split_pairs,
        "eligible_same_truth_pairs": same_truth_pairs,
        "false_split_rate": split_pairs / same_truth_pairs if same_truth_pairs else 0.0,
        "false_merge_pairs": merge_pairs,
        "eligible_same_prediction_pairs": same_prediction_pairs,
        "false_merge_rate": merge_pairs / same_prediction_pairs if same_prediction_pairs else 0.0,
    }


def _write_matrix(path: Path, samples: list[str], columns: list[str], values: dict[tuple[str, str], int]) -> None:
    def emit(handle):
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample", *columns])
        for sample in samples:
            writer.writerow([sample, *(values.get((sample, column), 0) for column in columns)])

    _atomic_text(path, emit)


def _run_in_workspace(
    workspace: Path,
    output_dir: Path,
    *,
    seed: int,
    transcript_count: int,
    sample_count: int,
    tool_count: int,
    artifact_rate: float,
) -> dict[str, object]:
    rng = random.Random(seed)
    truth, reference_v1, reference_v2 = generate_truth(transcript_count)
    contig_length = transcript_count * 3000 + 5000
    fasta = workspace / "simulation.fa"
    reference1_path = workspace / "reference-v1.gtf"
    reference2_path = workspace / "reference-v2.gtf"
    fasta.write_text(">chrSim\n" + ("ACGT" * math.ceil(contig_length / 4))[:contig_length] + "\n", encoding="utf-8")
    _write_annotation(reference1_path, reference_v1)
    _write_annotation(reference2_path, reference_v2)
    database = workspace / "simulation.sqlite"
    initialize_registry(
        database,
        fasta=fasta,
        assembly_name="TxID-simulation-v1",
        annotation=reference1_path,
        annotation_name="ref-v1",
    )
    add_annotation_context(
        database,
        annotation=reference2_path,
        annotation_name="ref-v2",
    )

    tools = ["IsoQuant-like", "TALON-like", "FLAIR-like", "StringTie-like"][:tool_count]
    samples = [f"sample-{index + 1:02d}" for index in range(sample_count)]
    observation_truth: dict[tuple[str, str, str], ObservationTruth] = {}
    mapping_rows: list[dict[str, str]] = []
    first_input: tuple[Path, str, str, dict[tuple[str, str, str], ObservationTruth]] | None = None
    for sample_index, sample in enumerate(samples):
        for tool_index, tool in enumerate(tools):
            observed_models: list[TranscriptModel] = []
            local_truth: dict[tuple[str, str, str], ObservationTruth] = {}
            for truth_item in truth:
                if rng.random() > 0.88:
                    continue
                perturbed = rng.random() < artifact_rate
                model = _perturb(truth_item.model) if perturbed else truth_item.model
                original_id = f"{tool.replace('-like', '')}.{sample}.{truth_item.truth_id.split(':')[1]}"
                model = replace(
                    model,
                    original_transcript_id=original_id,
                    original_gene_id=truth_item.gene_id or f"upstream-{truth_item.truth_id}",
                    source=tool,
                    attributes=(("truth_id", truth_item.truth_id),),
                )
                observed_models.append(model)
                if truth_item.category == "known" and perturbed:
                    expected_classification = "novel_in_known_gene"
                else:
                    expected_classification = truth_item.category
                record = ObservationTruth(
                    sample,
                    tool,
                    original_id,
                    truth_item.truth_id,
                    _structural_key(model),
                    expected_classification,
                    rng.randint(1, 100),
                    perturbed,
                )
                local_truth[(sample, tool, original_id)] = record
                observation_truth[(sample, tool, original_id)] = record
            suffix = "gtf" if tool_index % 2 == 0 else "gff3"
            input_path = workspace / f"{sample}.{tool}.{suffix}"
            content = _gtf(observed_models, rng) if suffix == "gtf" else _gff3(observed_models, rng)
            input_path.write_text(content, encoding="utf-8", newline="\n")
            mapping_path = workspace / f"{sample}.{tool}.mapping.tsv"
            import_file(
                database,
                input_path=input_path,
                sample=sample,
                tool=tool,
                annotation_name="ref-v1",
                output_gtf=workspace / f"{sample}.{tool}.txid.gtf",
                mapping_path=mapping_path,
                annotation_format=suffix,
            )
            with mapping_path.open("r", encoding="utf-8", newline="") as handle:
                mapping_rows.extend(csv.DictReader(handle, delimiter="\t"))
            if first_input is None:
                first_input = (input_path, suffix, tool, local_truth)

    if first_input is None:
        raise RuntimeError("simulation generated no caller observations")
    with Registry(database) as registry:
        fuzzy_map = registry.run_fuzzy(
            splice_tolerance=2,
            end_tolerance=15,
            software_version=__version__,
        )

    # Re-import one identical caller result under annotation v2 to test status
    # changes without altering exact structural identifiers.
    release_input, release_format, release_tool, _ = first_input
    release_mapping = workspace / "release-v2.mapping.tsv"
    import_file(
        database,
        input_path=release_input,
        sample=samples[0],
        tool=release_tool,
        annotation_name="ref-v2",
        output_gtf=workspace / "release-v2.txid.gtf",
        mapping_path=release_mapping,
        annotation_format=release_format,
    )
    with release_mapping.open("r", encoding="utf-8", newline="") as handle:
        v2_rows = {row["original_transcript_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    v1_first_rows = {
        row["original_transcript_id"]: row
        for row in mapping_rows
        if row["sample"] == samples[0] and row["tool"] == release_tool
    }
    shared_release_ids = sorted(set(v1_first_rows) & set(v2_rows))
    release_txid_invariant = all(
        v1_first_rows[item]["txid_form"] == v2_rows[item]["txid_form"]
        for item in shared_release_ids
    )
    status_changes = sum(
        v1_first_rows[item]["classification"] != v2_rows[item]["classification"]
        for item in shared_release_ids
    )

    structural_truth: list[str] = []
    underlying_truth: list[str] = []
    txid_predictions: list[str] = []
    raw_predictions: list[str] = []
    classification_correct = 0
    txid_matrix: dict[tuple[str, str], int] = defaultdict(int)
    raw_columns: set[str] = set()
    truth_to_txids: dict[str, set[str]] = defaultdict(set)
    for row in mapping_rows:
        key = (row["sample"], row["tool"], row["original_transcript_id"])
        record = observation_truth[key]
        structural_truth.append(record.structural_truth)
        underlying_truth.append(record.truth_id)
        txid_predictions.append(row["txid_form"])
        raw_predictions.append(row["original_transcript_id"])
        classification_correct += row["classification"] == record.expected_classification
        txid_matrix[(record.sample, row["txid_form"])] += record.count
        raw_columns.add(row["original_transcript_id"])
        truth_to_txids[record.truth_id].add(row["txid_form"])

    txid_columns = sorted(set(txid_predictions))
    _write_matrix(output_dir / "sample-by-txid.tsv", samples, txid_columns, txid_matrix)
    structural_metrics = _pairwise_error_rates(structural_truth, txid_predictions)
    raw_metrics = _pairwise_error_rates(structural_truth, raw_predictions)
    biological_metrics = _pairwise_error_rates(underlying_truth, txid_predictions)
    fragmented_truth = sum(len(values) > 1 for values in truth_to_txids.values())

    scaling = []
    identity_models = [item.model for item in truth]
    assembly = "sha256:" + "1" * 64
    for multiplier in [1, 2, 4]:
        start = time.perf_counter()
        items = identity_models * multiplier
        for model in items:
            identify(model, assembly)
        elapsed = time.perf_counter() - start
        scaling.append(
            {
                "items": len(items),
                "seconds": elapsed,
                "items_per_second": len(items) / elapsed if elapsed else None,
                "synthetic_resampled_scaling": True,
            }
        )

    summary: dict[str, object] = {
        "benchmark": "TxID truth-defined synthetic interoperability evaluation",
        "software_version": __version__,
        "environment": {
            "python": platform.python_version(),
            "sqlite": sqlite3.sqlite_version,
            "platform": platform.system(),
            "machine": platform.machine(),
        },
        "seed": seed,
        "parameters": {
            "truth_transcripts": transcript_count,
            "samples": sample_count,
            "simulated_callers": tools,
            "artifact_rate": artifact_rate,
            "fuzzy_splice_tolerance": 2,
            "fuzzy_end_tolerance": 15,
        },
        "observations": len(mapping_rows),
        "perturbed_observations": sum(item.perturbed for item in observation_truth.values()),
        "classification_accuracy": classification_correct / len(mapping_rows),
        "exact_structural_identity": structural_metrics,
        "raw_upstream_identifier_baseline": raw_metrics,
        "underlying_biological_model_fragmentation": {
            **biological_metrics,
            "fragmented_truth_transcripts": fragmented_truth,
            "fraction_fragmented": fragmented_truth / len(truth_to_txids),
            "interpretation": "Includes intentionally injected caller boundary errors; exact TxID correctly keeps distinct emitted structures distinct.",
        },
        "annotation_release_check": {
            "shared_observations": len(shared_release_ids),
            "exact_txid_invariant": release_txid_invariant,
            "classification_status_changes": status_changes,
        },
        "sample_by_transcript_matrix": {
            "raw_upstream_columns": len(raw_columns),
            "txid_exact_columns": len(txid_columns),
            "column_reduction": len(raw_columns) - len(txid_columns),
            "output": "sample-by-txid.tsv",
        },
        "fuzzy_clusters": len(set(fuzzy_map.values())),
        "scaling": scaling,
        "comparator_status": {
            "gffcompare": "outside this simulation run; see external-comparison-summary.json",
            "isoSeQL": "outside this simulation run; see external-comparison-summary.json",
            "TALON_shared_database": "outside this simulation run; see external-comparison-summary.json",
        },
    }
    write_json(output_dir / "simulation-summary.json", summary)
    return summary


def run_simulation(
    output_dir: str | Path,
    *,
    seed: int = 20260730,
    transcript_count: int = 120,
    sample_count: int = 6,
    tool_count: int = 4,
    artifact_rate: float = 0.05,
    work_dir: str | Path | None = None,
) -> dict[str, object]:
    if sample_count < 1:
        raise ValueError("sample_count must be positive")
    if not 1 <= tool_count <= 4:
        raise ValueError("tool_count must be between 1 and 4")
    if not 0 <= artifact_rate <= 1:
        raise ValueError("artifact_rate must be between 0 and 1")
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    tracemalloc.start()
    started = time.perf_counter()
    if work_dir is None:
        with tempfile.TemporaryDirectory(prefix="txid-simulation-") as temporary:
            summary = _run_in_workspace(
                Path(temporary),
                destination,
                seed=seed,
                transcript_count=transcript_count,
                sample_count=sample_count,
                tool_count=tool_count,
                artifact_rate=artifact_rate,
            )
    else:
        workspace = Path(work_dir)
        workspace.mkdir(parents=True, exist_ok=True)
        summary = _run_in_workspace(
            workspace,
            destination,
            seed=seed,
            transcript_count=transcript_count,
            sample_count=sample_count,
            tool_count=tool_count,
            artifact_rate=artifact_rate,
        )
    elapsed = time.perf_counter() - started
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    summary["end_to_end_resources"] = {
        "seconds": elapsed,
        "peak_memory_mib": peak / (1024 * 1024),
    }
    write_json(destination / "simulation-summary.json", summary)
    return summary
