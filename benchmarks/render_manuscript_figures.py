#!/usr/bin/env python3
"""Render publication figures from checked TxID benchmark summaries.

All numbered figures are generated with Matplotlib.  The script reads only the
checked JSON summaries and writes editable SVG, 300-dpi PNG, and a tabular source
data file.  It does not recalculate benchmark metrics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Iterable

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.ticker import FuncFormatter, PercentFormatter


# Okabe-Ito colour-blind-safe palette. Grey is used for context and blue for
# TxID; comparator colours are stable across figures.
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#D55E00"
PURPLE = "#CC79A7"
SKY = "#56B4E9"
BLACK = "#222222"
MID_GREY = "#6B6B6B"
LIGHT_GREY = "#D9D9D9"
VERY_LIGHT_GREY = "#F2F2F2"

METHOD_COLOURS = {
    "TxID": BLUE,
    "TxID exact": BLUE,
    "gffcompare": ORANGE,
    "gffcompare exact": ORANGE,
    "gffcompare splice": "#B56F00",
    "isoSeQL": GREEN,
    "isoSeQL exact ends": GREEN,
    "isoSeQL junction": "#007C5A",
    "TALON": PURPLE,
}


def configure_style() -> None:
    """Apply a restrained, journal-oriented plotting style."""

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.0,
            "axes.titlesize": 9.0,
            "axes.labelsize": 8.0,
            "axes.linewidth": 0.7,
            "axes.edgecolor": BLACK,
            "axes.titleweight": "normal",
            "xtick.labelsize": 7.2,
            "ytick.labelsize": 7.2,
            "xtick.major.width": 0.7,
            "ytick.major.width": 0.7,
            "xtick.major.size": 3.0,
            "ytick.major.size": 3.0,
            "legend.fontsize": 7.2,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "svg.fonttype": "none",
            "svg.hashsalt": "txid-manuscript-figures-v1",
            "pdf.fonttype": 42,
        }
    )


def panel_label(ax: Axes, label: str) -> None:
    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        ha="left",
        va="bottom",
        clip_on=False,
    )


def clean_axis(ax: Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_axisbelow(True)


def save_figure(fig: Figure, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        format="svg",
        metadata={"Date": None, "Creator": "TxID Matplotlib figure renderer"},
    )
    fig.savefig(
        output_path.with_suffix(".png"),
        format="png",
        dpi=300,
        metadata={"Software": "TxID Matplotlib figure renderer"},
    )
    plt.close(fig)


def draw_transcript(
    ax: Axes,
    exons: Iterable[tuple[float, float]],
    y: float,
    *,
    colour: str = BLUE,
    height: float = 0.16,
    linewidth: float = 1.2,
) -> None:
    exon_list = list(exons)
    for left, right in zip(exon_list, exon_list[1:]):
        ax.plot([left[1], right[0]], [y, y], color=BLACK, linewidth=linewidth, zorder=1)
    for start, end in exon_list:
        ax.add_patch(
            Rectangle(
                (start, y - height / 2),
                end - start,
                height,
                facecolor=colour,
                edgecolor=BLACK,
                linewidth=0.6,
                zorder=2,
            )
        )


def simple_box(
    ax: Axes,
    x: float,
    y: float,
    width: float,
    height: float,
    text: str,
    *,
    facecolor: str = "white",
    edgecolor: str = MID_GREY,
    fontsize: float = 7.2,
    weight: str = "normal",
) -> None:
    ax.add_patch(
        Rectangle(
            (x, y),
            width,
            height,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=0.8,
        )
    )
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        linespacing=1.25,
    )


def arrow(ax: Axes, start: tuple[float, float], end: tuple[float, float]) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=8,
            linewidth=0.8,
            color=MID_GREY,
            shrinkA=2,
            shrinkB=2,
        )
    )


def horizontal_bars(
    ax: Axes,
    labels: list[str],
    values: list[float],
    colours: list[str],
    *,
    value_labels: list[str] | None = None,
    xlim: tuple[float, float] | None = None,
) -> None:
    y = list(range(len(labels)))
    ax.barh(y, values, color=colours, edgecolor="none", height=0.58, zorder=2)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    if xlim:
        ax.set_xlim(*xlim)
    labels_to_draw = value_labels or [f"{value:g}" for value in values]
    xmax = ax.get_xlim()[1]
    for yi, value, value_label in zip(y, values, labels_to_draw):
        offset = 0.018 * xmax
        ax.text(value + offset, yi, value_label, va="center", ha="left", fontsize=7.2)
    clean_axis(ax)


def render_workflow(output_path: Path) -> None:
    fig = plt.figure(figsize=(7.2, 3.20))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.04, 1.12, 1.05], wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.set_title("Inputs", loc="left", pad=8)
    ax.text(0.02, 0.84, "sample A · local ID A_17", fontsize=7.2)
    draw_transcript(ax, [(0.05, 0.22), (0.40, 0.56), (0.76, 0.96)], 0.70, colour=SKY)
    ax.text(0.02, 0.47, "sample B · local ID B_204", fontsize=7.2)
    draw_transcript(ax, [(0.05, 0.22), (0.40, 0.56), (0.76, 0.96)], 0.33, colour=SKY)
    ax.text(0.02, 0.08, "Same assembly and exon coordinates", fontsize=7, color=MID_GREY)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Canonical identity", loc="left", pad=8)
    fields = ["assembly fingerprint", "primary contig + strand", "ordered intron chain", "TSS and TES"]
    for i, field in enumerate(fields):
        y = 0.82 - i * 0.145
        ax.text(0.04, y, field, va="center", fontsize=7.2)
        ax.plot([0.68, 0.92], [y, y], color=LIGHT_GREY, linewidth=2.2)
    arrow(ax, (0.50, 0.30), (0.50, 0.18))
    simple_box(ax, 0.08, 0.02, 0.84, 0.14, "txid:TF1.<digest>", facecolor=VERY_LIGHT_GREY, edgecolor=BLUE, fontsize=7.2, weight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax = fig.add_subplot(gs[0, 2])
    panel_label(ax, "C")
    ax.set_title("Outputs", loc="left", pad=8)
    outputs = [
        ("rewritten GTF/GFF3", "known reference IDs retained"),
        ("mapping table", "local ID ↔ exact TxID"),
        ("cohort catalogue", "stable matrix join key"),
        ("SQLite registry", "context and provenance"),
    ]
    for i, (name, note) in enumerate(outputs):
        y = 0.78 - i * 0.205
        simple_box(ax, 0.04, y, 0.92, 0.14, name, facecolor="white", edgecolor=LIGHT_GREY, fontsize=7.1, weight="bold")
        ax.text(0.50, y - 0.025, note, ha="center", va="top", fontsize=6.5, color=MID_GREY)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.text(
        0.5,
        0.005,
        "Annotation release, classification, sample and upstream tool are metadata—not identity inputs.",
        ha="center",
        va="bottom",
        fontsize=7.0,
        color=MID_GREY,
    )
    save_figure(fig, output_path)


def render_benchmark(data: dict[str, Any], output_path: Path) -> None:
    exact = data["exact_structural_identity"]
    caller = data["raw_upstream_identifier_baseline"]
    biological = data["underlying_biological_model_fragmentation"]
    matrix = data["sample_by_transcript_matrix"]
    annotation = data["annotation_release_check"]

    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.05))
    fig.subplots_adjust(left=0.15, right=0.98, bottom=0.24, top=0.82, wspace=0.95)

    ax = axes[0]
    panel_label(ax, "A")
    ax.set_title("Same-target pairs split", loc="left")
    labels = ["TxID\n(emitted form)", "Original\ncaller ID", "TxID vs\nsource model*"]
    values = [100 * exact["false_split_rate"], 100 * caller["false_split_rate"], 100 * biological["false_split_rate"]]
    horizontal_bars(ax, labels, values, [BLUE, MID_GREY, RED], value_labels=[f"{v:.2f}%" for v in values], xlim=(0, 112))
    ax.set_xlabel("False-split rate (%)")

    ax = axes[1]
    panel_label(ax, "B")
    ax.set_title("Matrix columns", loc="left")
    values = [matrix["raw_upstream_columns"], matrix["txid_exact_columns"]]
    horizontal_bars(ax, ["Original caller IDs", "Exact TxID"], values, [MID_GREY, BLUE], value_labels=[f"{v:,}" for v in values], xlim=(0, 2800))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))
    ax.set_xlabel("Distinct columns")

    ax = axes[2]
    panel_label(ax, "C")
    ax.set_title("Annotation update", loc="left")
    shared = annotation["shared_observations"]
    values = [0, annotation["classification_status_changes"] / shared * 100]
    horizontal_bars(ax, ["Exact TxID", "Classification"], values, [BLUE, ORANGE], value_labels=["0 / 109", f"22 / {shared}"], xlim=(0, 26))
    ax.set_xlabel("Shared observations (%)")

    fig.text(0.02, 0.035, "*Includes deliberately injected 1-bp boundary errors; distinct emitted structures are expected to remain distinct.", fontsize=6.5, color=MID_GREY)
    save_figure(fig, output_path)


def render_identity_semantics(output_path: Path) -> None:
    fig = plt.figure(figsize=(7.2, 3.85))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.34)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.set_title("Exact identifier families", loc="left", pad=8)
    rows = [
        (0.82, [(0.07, 0.22), (0.42, 0.56), (0.76, 0.94)], "reference form", "SC1=a · TF1=x"),
        (0.61, [(0.02, 0.22), (0.42, 0.56), (0.76, 0.99)], "different ends", "SC1=a · TF1=y"),
        (0.40, [(0.07, 0.22), (0.43, 0.56), (0.76, 0.94)], "1-bp splice shift", "SC1=b · TF1=z"),
        (0.19, [(0.31, 0.70)], "single exon", "SE1=s"),
    ]
    for y, exons, label, ids in rows:
        draw_transcript(ax, exons, y, colour=BLUE if label != "1-bp splice shift" else RED, height=0.11)
        ax.text(-0.02, y, label, ha="right", va="center", fontsize=7.2)
        ax.text(1.01, y, ids, ha="left", va="center", fontsize=7.1, family="monospace")
    ax.plot([0.425, 0.425], [0.34, 0.46], color=RED, linewidth=0.8, linestyle="--")
    ax.text(0.43, 0.335, "+1 bp", ha="center", va="top", fontsize=6.5, color=RED)
    ax.set_xlim(-0.27, 1.33)
    ax.set_ylim(0.05, 0.96)
    ax.axis("off")

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Annotation is versioned metadata", loc="left", pad=8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    simple_box(ax, 0.08, 0.77, 0.84, 0.12, "exact structure on GRCh38", facecolor=VERY_LIGHT_GREY, edgecolor=BLUE, weight="bold")
    arrow(ax, (0.50, 0.76), (0.50, 0.64))
    ax.text(0.50, 0.58, "txid:TF1.<digest>", ha="center", va="center", fontsize=7.6, family="monospace", fontweight="bold", color=BLUE)
    ax.plot([0.50, 0.25], [0.52, 0.40], color=MID_GREY, linewidth=0.8)
    ax.plot([0.50, 0.75], [0.52, 0.40], color=MID_GREY, linewidth=0.8)
    simple_box(ax, 0.02, 0.16, 0.44, 0.23, "GENCODE v29\nnovel_in_known_gene\nreference ID:\nnone", facecolor="white", edgecolor=LIGHT_GREY, fontsize=5.5)
    simple_box(ax, 0.54, 0.16, 0.44, 0.23, "GENCODE v49\nknown\nreference ID:\nretained", facecolor="white", edgecolor=LIGHT_GREY, fontsize=5.5)
    ax.text(0.50, 0.07, "The exact TxID is unchanged.", ha="center", fontsize=7.1, color=MID_GREY)
    ax.axis("off")
    save_figure(fig, output_path)


def render_registry_and_verification(output_path: Path, test_count: int) -> None:
    fig = plt.figure(figsize=(7.2, 3.85))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1], wspace=0.30)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.set_title("Registry data model", loc="left", pad=8)
    boxes = {
        "reference": (0.02, 0.79, 0.41, 0.12, "reference context"),
        "annotation": (0.57, 0.79, 0.41, 0.12, "annotation context"),
        "structure": (0.24, 0.51, 0.52, 0.15, "exact structural object\ncanonical JSON + full digest"),
        "observation": (0.02, 0.23, 0.41, 0.13, "sample observation\nupstream provenance"),
        "interpretation": (0.57, 0.23, 0.41, 0.13, "reference interpretation\nclass + retained aliases"),
        "optional": (0.29, 0.02, 0.42, 0.11, "FC groups / GL accessions"),
    }
    for key, (x, y, w, h, label) in boxes.items():
        fill = VERY_LIGHT_GREY if key == "structure" else "white"
        edge = BLUE if key == "structure" else LIGHT_GREY
        simple_box(ax, x, y, w, h, label, facecolor=fill, edgecolor=edge, fontsize=6.2, weight="bold" if key == "structure" else "normal")
    arrow(ax, (0.23, 0.78), (0.40, 0.67))
    arrow(ax, (0.77, 0.78), (0.60, 0.67))
    arrow(ax, (0.39, 0.50), (0.23, 0.37))
    arrow(ax, (0.61, 0.50), (0.77, 0.37))
    arrow(ax, (0.37, 0.22), (0.46, 0.14))
    arrow(ax, (0.63, 0.22), (0.54, 0.14))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Import and verification", loc="left", pad=8)
    steps = [
        "Parse complete input",
        "Validate assembly and annotation context",
        "Derive and independently\nrecompute digests",
        "Commit registry transaction",
        "Sort and serialize outputs deterministically",
    ]
    for i, step in enumerate(steps):
        y = 0.87 - i * 0.135
        ax.text(0.09, y, f"{i + 1}", ha="center", va="center", fontsize=7, color="white", bbox={"boxstyle": "square,pad=0.22", "facecolor": BLUE, "edgecolor": BLUE})
        ax.text(0.18, y, step, ha="left", va="center", fontsize=6.7, linespacing=1.2)
        if i < len(steps) - 1:
            ax.plot([0.09, 0.09], [y - 0.035, y - 0.10], color=LIGHT_GREY, linewidth=0.8)
    ax.plot([0.03, 0.97], [0.26, 0.26], color=LIGHT_GREY, linewidth=0.8)
    ax.text(0.03, 0.21, f"Automated verification ({test_count} tests)", fontsize=7.2, fontweight="bold", ha="left")
    categories = ["parsing and coordinates", "identity and annotation", "transactions and retries", "invariance and conformance"]
    for i, category in enumerate(categories):
        ax.text(0.03, 0.15 - i * 0.045, f"— {category}", fontsize=6.3, ha="left", va="center")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    save_figure(fig, output_path)


def render_external_comparison(data: dict[str, Any], output_path: Path) -> None:
    benchmark = data["benchmark"]
    results = data["results"]
    methods = ["TxID", "gffcompare", "isoSeQL", "TALON"]
    colours = [METHOD_COLOURS[name] for name in methods]

    fig = plt.figure(figsize=(7.2, 5.05))
    gs = fig.add_gridspec(2, 2, hspace=0.58, wspace=0.55)

    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.set_title("Exact-form matrix columns", loc="left")
    values = [int(results[m]["exact_form"]["predicted_groups"]) for m in methods]
    horizontal_bars(ax, methods, values, colours, value_labels=[str(v) for v in values], xlim=(0, 215))
    ax.set_xlabel("Distinct columns")

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Exact-form false merges", loc="left")
    values = [100 * float(results[m]["exact_form"]["false_merge_rate"]) for m in methods]
    horizontal_bars(ax, methods, values, colours, value_labels=[f"{v:.3f}%" for v in values], xlim=(0, 1.75))
    ax.set_xlabel("False-merge rate (%)")

    cgs = gs[1, 0].subgridspec(2, 1, hspace=0.85)
    comparator_methods = ["gffcompare", "isoSeQL", "TALON"]
    comparator_colours = [METHOD_COLOURS[name] for name in comparator_methods]
    ax = fig.add_subplot(cgs[0, 0])
    panel_label(ax, "C")
    ax.set_title("Measured workflow stages", loc="left")
    values = [float(results[m]["runtime"]["elapsed_seconds"]) for m in comparator_methods]
    horizontal_bars(ax, comparator_methods, values, comparator_colours, value_labels=[f"{v:.2f} s" for v in values], xlim=(0, 3.0))
    ax.set_xlabel("Wall time (s)")
    ax = fig.add_subplot(cgs[1, 0])
    values = [float(results[m]["runtime"]["max_rss_kib"]) / 1024 for m in comparator_methods]
    horizontal_bars(ax, comparator_methods, values, comparator_colours, value_labels=[f"{v:.1f} MiB" for v in values], xlim=(0, 101))
    ax.set_xlabel("Peak RSS (MiB)")

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "D")
    ax.set_title("Literal labels changed after permutation", loc="left")
    observations = int(benchmark["observations"])
    values = [100 * int(results[m]["order_invariance"]["literal_label_changed_observations"]) / observations for m in methods]
    horizontal_bars(ax, methods, values, colours, value_labels=[f"{v:.1f}%" for v in values], xlim=(0, 103))
    ax.set_xlabel("Observations (%)")
    fig.text(0.02, 0.01, "Measured stages differ by workflow; SQANTI3 is excluded from isoSeQL timing. Structural partitions changed for 0 observations in all methods.", fontsize=6.4, color=MID_GREY)
    save_figure(fig, output_path)


def _ci_error(metric: dict[str, Any], key: str) -> tuple[float, float, float]:
    value = 100 * float(metric[f"{key}_rate"])
    lower = 100 * float(metric[f"{key}_ci95_lower"])
    upper = 100 * float(metric[f"{key}_ci95_upper"])
    return value, value - lower, upper - value


def render_encode_gtf_full(data: dict[str, Any], output_path: Path) -> None:
    evaluations = data["evaluations"]
    order = data["order_comparisons"]
    rows = [
        ("TxID exact", evaluations["txid-sorted"]),
        ("isoSeQL exact ends", evaluations["isoseql-exact-sorted"]),
        ("gffcompare exact", evaluations["gffcompare-sorted"]),
        ("gffcompare splice", evaluations["gffcompare-splice-sorted"]),
        ("isoSeQL junction", evaluations["isoseql-junction-sorted"]),
    ]
    labels = [name for name, _ in rows]
    colours = [METHOD_COLOURS[name] for name in labels]

    fig = plt.figure(figsize=(7.2, 5.55))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.1, 1.0], hspace=0.54, wspace=0.52)

    ags = gs[0, 0].subgridspec(1, 2, wspace=0.42)
    ax_merge = fig.add_subplot(ags[0, 0])
    ax_merge.text(-0.20, 1.06, "A", transform=ax_merge.transAxes, fontsize=10, fontweight="bold", ha="left", va="bottom", clip_on=False)
    ax_merge.set_title("False merge", loc="left")
    y = list(range(len(rows)))
    merge = [_ci_error(metric, "false_merge") for _, metric in rows]
    for yi, (value, lo, hi), colour in zip(y, merge, colours):
        ax_merge.errorbar(value, yi, xerr=[[lo], [hi]], fmt="o", markersize=4, color=colour, ecolor=colour, capsize=2, linewidth=0.9)
        ax_merge.text(value + 1.8, yi, f"{value:.2f}%" if value else "0", va="center", fontsize=6.2)
    ax_merge.set_yticks(y, labels)
    ax_merge.invert_yaxis()
    ax_merge.set_xlim(-2, 79)
    ax_merge.set_xlabel("Rate (%)")
    clean_axis(ax_merge)

    ax_split = fig.add_subplot(ags[0, 1])
    ax_split.set_title("False split", loc="left")
    split = [_ci_error(metric, "false_split") for _, metric in rows]
    for yi, (value, lo, hi), colour in zip(y, split, colours):
        ax_split.errorbar(value, yi, xerr=[[lo], [hi]], fmt="o", markersize=4, color=colour, ecolor=colour, capsize=2, linewidth=0.9)
        ax_split.text(value + 0.08, yi, f"{value:.3f}%" if value else "0", va="center", fontsize=6.2)
    ax_split.set_yticks(y, [""] * len(labels))
    ax_split.invert_yaxis()
    ax_split.set_xlim(-0.12, 3.25)
    ax_split.set_xlabel("Rate (%)")
    clean_axis(ax_split)

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Pair relations changed after input permutation", loc="left")
    order_rows = [
        ("TxID exact", int(order["txid"]["exact_form_partition"]["pair_relation_changes"]), BLUE),
        ("isoSeQL exact ends", int(order["isoseql_exact_ends"]["pair_relation_changes"]), GREEN),
        ("isoSeQL junction", int(order["isoseql_common_junction"]["pair_relation_changes"]), "#007C5A"),
        ("gffcompare", int(order["gffcompare_exact"]["pair_relation_changes"]), ORANGE),
        ("TxID GL locus", int(order["txid"]["gene_locus_partition"]["pair_relation_changes"]), MID_GREY),
    ]
    values = [v for _, v, _ in order_rows]
    horizontal_bars(ax, [x[0] for x in order_rows], values, [x[2] for x in order_rows], value_labels=[f"{v:,}" for v in values], xlim=(0, 465000))
    ax.set_xlabel("Changed observation-pair relations")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x / 1000):,}k" if x else "0"))

    ax = fig.add_subplot(gs[1, 0])
    panel_label(ax, "C")
    ax.set_title("Six-sample matrix dimensions", loc="left")
    values = [int(metric["matrix_columns"]) for _, metric in rows]
    sparsity = [100 * float(metric["matrix_sparsity"]) for _, metric in rows]
    horizontal_bars(ax, labels, values, colours, value_labels=[f"{v:,} · {s:.1f}% sparse" for v, s in zip(values, sparsity)], xlim=(0, 555000))
    ax.set_xlabel("Distinct transcript groups")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x / 1000):,}k" if x else "0"))

    ax = fig.add_subplot(gs[1, 1])
    panel_label(ax, "D")
    ax.set_title("Data and evaluation scope", loc="left")
    ax.axis("off")
    scope_rows = [
        ("Input", "6 ENCODE WTC11 TALON GTFs"),
        ("Models", "538,804 fixed GTF models"),
        ("Reference", "GRCh38; GENCODE v29 and v49"),
        ("Bootstrap", "2,000 replicates; structural units"),
        ("Emitted", "538,804 (TxID/isoSeQL); 538,802 (gffcompare)"),
        ("Claim", "released-model identity only"),
    ]
    table = ax.table(cellText=scope_rows, colWidths=[0.26, 0.74], cellLoc="left", loc="upper left", bbox=[0, 0.03, 1, 0.90])
    table.auto_set_font_size(False)
    table.set_fontsize(6.1)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor(LIGHT_GREY)
        cell.set_linewidth(0.6)
        cell.set_facecolor(VERY_LIGHT_GREY if col == 0 else "white")
        if col == 0:
            cell.get_text().set_fontweight("bold")

    fig.text(0.02, 0.005, "Points show pairwise rates; error bars are 95% bootstrap intervals. Targets differ where explicitly labelled exact, splice or junction.", fontsize=6.6, color=MID_GREY)
    save_figure(fig, output_path)


def render_encode_gtf_pilot(data: dict[str, Any], output_path: Path) -> None:
    """Render the retained, non-numbered pilot result in the same visual style."""

    evaluations = data["evaluations"]
    order = data["order_comparisons"]
    rows = [
        ("TxID exact", evaluations["txid-sorted"]),
        ("isoSeQL exact ends", evaluations["isoseql-exact-sorted"]),
        ("gffcompare exact", evaluations["gffcompare-sorted"]),
        ("isoSeQL junction", evaluations["isoseql-junction-sorted"]),
    ]
    labels = [name for name, _ in rows]
    colours = [METHOD_COLOURS[name] for name in labels]
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 4.7))
    fig.subplots_adjust(left=0.20, right=0.97, bottom=0.12, top=0.92, hspace=0.62, wspace=0.78)

    for ax, key, title, xmax, panel in (
        (axes[0, 0], "false_merge", "False merge", 6.2, "A"),
        (axes[0, 1], "false_split", "False split", 0.14, "B"),
    ):
        panel_label(ax, panel)
        ax.set_title(title, loc="left")
        y = list(range(len(rows)))
        values = [_ci_error(metric, key) for _, metric in rows]
        for yi, (value, lo, hi), colour in zip(y, values, colours):
            ax.errorbar(value, yi, xerr=[[lo], [hi]], fmt="o", markersize=4, color=colour, ecolor=colour, capsize=2, linewidth=0.9)
            ax.text(value + xmax * 0.025, yi, f"{value:.3f}%" if value else "0", va="center", fontsize=6.4)
        ax.set_yticks(y, labels if panel == "A" else [""] * len(labels))
        ax.invert_yaxis()
        ax.set_xlim(-0.02 * xmax, xmax)
        ax.set_xlabel("Rate (%)")
        clean_axis(ax)

    ax = axes[1, 0]
    panel_label(ax, "C")
    ax.set_title("Pilot matrix dimensions", loc="left")
    values = [int(metric["matrix_columns"]) for _, metric in rows]
    sparsity = [100 * float(metric["matrix_sparsity"]) for _, metric in rows]
    horizontal_bars(ax, labels, values, colours, value_labels=[f"{v / 1000:.1f}k · {s:.1f}%" for v, s in zip(values, sparsity)], xlim=(0, 31000))
    ax.set_xlabel("Distinct transcript groups")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x / 1000):,}k" if x else "0"))

    ax = axes[1, 1]
    panel_label(ax, "D")
    ax.set_title("Permutation effects", loc="left")
    order_rows = [
        ("TxID exact", int(order["txid"]["exact_form_partition"]["pair_relation_changes"]), BLUE),
        ("isoSeQL exact ends", int(order["isoseql_exact_ends"]["pair_relation_changes"]), GREEN),
        ("isoSeQL junction", int(order["isoseql_common_junction"]["pair_relation_changes"]), "#007C5A"),
        ("gffcompare", int(order["gffcompare_exact"]["pair_relation_changes"]), ORANGE),
        ("TxID GL locus", int(order["txid"]["gene_locus_partition"]["pair_relation_changes"]), MID_GREY),
    ]
    values = [value for _, value, _ in order_rows]
    horizontal_bars(ax, [name for name, _, _ in order_rows], values, [colour for _, _, colour in order_rows], value_labels=[f"{value:,}" for value in values], xlim=(0, 12800))
    ax.set_xlabel("Changed observation-pair relations")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x / 1000):,}k" if x else "0"))
    save_figure(fig, output_path)


def render_scaling(data: dict[str, Any], output_path: Path) -> None:
    scaling = data["scaling"]
    resources = data["end_to_end_resources"]
    items = [int(row["items"]) for row in scaling]
    throughput = [float(row["items_per_second"]) for row in scaling]

    fig = plt.figure(figsize=(7.2, 2.65))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1.0], wspace=0.40)
    ax = fig.add_subplot(gs[0, 0])
    panel_label(ax, "A")
    ax.set_title("Identity-only resampling", loc="left")
    ax.plot(items, throughput, color=BLUE, marker="o", markersize=4, linewidth=1.2)
    ax.set_xticks(items)
    ax.set_xlabel("Resampled transcript models")
    ax.set_ylabel("Models processed per second")
    ax.set_ylim(10500, 12000)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{int(x):,}"))
    clean_axis(ax)

    ax = fig.add_subplot(gs[0, 1])
    panel_label(ax, "B")
    ax.set_title("Full synthetic run", loc="left")
    ax.axis("off")
    resource_rows = [
        ("Wall time", f'{float(resources["seconds"]):.2f} s'),
        ("Peak traced memory", f'{float(resources["peak_memory_mib"]):.2f} MiB'),
        ("Host", "development host"),
    ]
    for i, (label, value) in enumerate(resource_rows):
        y = 0.62 - i * 0.19
        ax.text(0.02, y, label, fontsize=7.1, ha="left", va="center")
        ax.text(0.98, y, value, fontsize=7.1, ha="right", va="center", fontweight="bold")
        ax.plot([0.02, 0.98], [y - 0.075, y - 0.075], color=LIGHT_GREY, linewidth=0.6)
    ax.text(0.02, 0.02, "Resampled rows are scaling checks,\nnot independent biological samples.", fontsize=6.6, color=MID_GREY, va="top")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    save_figure(fig, output_path)


def write_source_data(
    data: dict[str, Any],
    external: dict[str, Any],
    encode_benchmark: dict[str, Any],
    output_path: Path,
) -> None:
    exact = data["exact_structural_identity"]
    raw = data["raw_upstream_identifier_baseline"]
    biological = data["underlying_biological_model_fragmentation"]
    matrix = data["sample_by_transcript_matrix"]
    annotation = data["annotation_release_check"]
    resources = data["end_to_end_resources"]
    rows: list[tuple[str, str, str, object, str]] = [
        ("Figure 2", "A", "TxID exact false-split rate", exact["false_split_rate"], "proportion"),
        ("Figure 2", "A", "Caller-ID false-split rate", raw["false_split_rate"], "proportion"),
        ("Figure 2", "A", "Underlying-truth fragmentation rate", biological["false_split_rate"], "proportion"),
        ("Figure 2", "B", "Original caller-ID columns", matrix["raw_upstream_columns"], "columns"),
        ("Figure 2", "B", "Exact TxID columns", matrix["txid_exact_columns"], "columns"),
        ("Figure 2", "C", "Annotation-release shared observations", annotation["shared_observations"], "observations"),
        ("Figure 2", "C", "Annotation classifications changed", annotation["classification_status_changes"], "observations"),
        ("Figure S1", "B", "End-to-end wall time", resources["seconds"], "seconds"),
        ("Figure S1", "B", "Peak Python-traced memory", resources["peak_memory_mib"], "MiB"),
    ]
    for row in data["scaling"]:
        rows.append(("Figure S1", "A", f'Throughput at {row["items"]} models', row["items_per_second"], "models/second"))

    external_results = external["results"]
    observations = int(external["benchmark"]["observations"])
    for method in ("TxID", "gffcompare", "isoSeQL", "TALON"):
        item = external_results[method]
        rows.extend(
            [
                ("Figure 5", "A", f"{method} exact-form columns", item["exact_form"]["predicted_groups"], "columns"),
                ("Figure 5", "B", f"{method} false-merge rate", item["exact_form"]["false_merge_rate"], "proportion"),
                ("Figure 5", "D", f"{method} literal-label changes", item["order_invariance"]["literal_label_changed_observations"], "observations"),
                ("Figure 5", "D", f"{method} literal-label change rate", int(item["order_invariance"]["literal_label_changed_observations"]) / observations, "proportion"),
            ]
        )
        if method != "TxID":
            rows.extend(
                [
                    ("Figure 5", "C", f"{method} measured-stage wall time", item["runtime"]["elapsed_seconds"], "seconds"),
                    ("Figure 5", "C", f"{method} measured-stage peak RSS", float(item["runtime"]["max_rss_kib"]) / 1024, "MiB"),
                ]
            )

    evaluations = encode_benchmark["evaluations"]
    encode_rows = [
        ("TxID exact", evaluations["txid-sorted"]),
        ("isoSeQL exact ends", evaluations["isoseql-exact-sorted"]),
        ("gffcompare exact", evaluations["gffcompare-sorted"]),
        ("gffcompare splice", evaluations["gffcompare-splice-sorted"]),
        ("isoSeQL junction", evaluations["isoseql-junction-sorted"]),
    ]
    for method, item in encode_rows:
        for key, label in (("false_merge", "false-merge"), ("false_split", "false-split")):
            rows.extend(
                [
                    ("Figure 6", "A", f"{method} {label} rate", item[f"{key}_rate"], "proportion"),
                    ("Figure 6", "A", f"{method} {label} CI95 lower", item[f"{key}_ci95_lower"], "proportion"),
                    ("Figure 6", "A", f"{method} {label} CI95 upper", item[f"{key}_ci95_upper"], "proportion"),
                ]
            )
        rows.extend(
            [
                ("Figure 6", "C", f"{method} matrix columns", item["matrix_columns"], "columns"),
                ("Figure 6", "C", f"{method} matrix sparsity", item["matrix_sparsity"], "proportion"),
            ]
        )
    order = encode_benchmark["order_comparisons"]
    rows.extend(
        [
            ("Figure 6", "B", "TxID exact pair-relation changes", order["txid"]["exact_form_partition"]["pair_relation_changes"], "pairs"),
            ("Figure 6", "B", "TxID GL pair-relation changes", order["txid"]["gene_locus_partition"]["pair_relation_changes"], "pairs"),
            ("Figure 6", "B", "gffcompare pair-relation changes", order["gffcompare_exact"]["pair_relation_changes"], "pairs"),
            ("Figure 6", "B", "isoSeQL exact pair-relation changes", order["isoseql_exact_ends"]["pair_relation_changes"], "pairs"),
            ("Figure 6", "B", "isoSeQL junction pair-relation changes", order["isoseql_common_junction"]["pair_relation_changes"], "pairs"),
            ("Figure 6", "D", "gffcompare omitted observations", evaluations["gffcompare-sorted"]["omitted_observations"], "observations"),
        ]
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        "figure\tpanel\tmetric\tvalue\tunit\n"
        + "".join("\t".join(str(value) for value in row) + "\n" for row in rows),
        encoding="utf-8",
    )


def render_all(
    summary_path: Path,
    external_summary_path: Path,
    encode_benchmark_path: Path,
    encode_pilot_path: Path | None,
    output_dir: Path,
    test_count: int,
) -> None:
    configure_style()
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    external = json.loads(external_summary_path.read_text(encoding="utf-8"))
    encode_benchmark = json.loads(encode_benchmark_path.read_text(encoding="utf-8"))
    render_workflow(output_dir / "figure1-txid-workflow.svg")
    render_benchmark(data, output_dir / "figure2-simulation.svg")
    render_identity_semantics(output_dir / "figure3-identity-semantics.svg")
    render_registry_and_verification(output_dir / "figure4-registry-and-verification.svg", test_count)
    render_external_comparison(external, output_dir / "figure5-external-comparison.svg")
    render_encode_gtf_full(encode_benchmark, output_dir / "figure6-encode-gtf-full.svg")
    if encode_pilot_path is not None and encode_pilot_path.is_file():
        pilot = json.loads(encode_pilot_path.read_text(encoding="utf-8"))
        render_encode_gtf_pilot(pilot, output_dir / "figure6-encode-gtf-pilot.svg")
    render_scaling(data, output_dir / "figure-s1-scaling.svg")
    write_source_data(data, external, encode_benchmark, output_dir / "source-data.tsv")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=Path("benchmarks/results/simulation-summary.json"))
    parser.add_argument("--external-summary", type=Path, default=Path("benchmarks/results/external-comparison-summary.json"))
    parser.add_argument(
        "--encode-gtf",
        "--encode-pilot",
        dest="encode_gtf",
        type=Path,
        default=Path("benchmarks/results/encode-gtf-full/identity-partitions.json"),
    )
    parser.add_argument(
        "--encode-pilot-summary",
        type=Path,
        default=Path("benchmarks/results/encode-gtf-pilot/identity-partitions.json"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("manuscript/figures"))
    parser.add_argument("--test-count", type=int, default=56)
    args = parser.parse_args()
    render_all(
        args.summary,
        args.external_summary,
        args.encode_gtf,
        args.encode_pilot_summary,
        args.output_dir,
        args.test_count,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
