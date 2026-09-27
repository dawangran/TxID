"""Three-panel TxID workflow, annotation transitions and cohort dimensions."""
from __future__ import annotations

import json
from collections import Counter
from pathlib import Path as FilePath

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.patches import FancyBboxPatch, PathPatch, Rectangle

WIDTH_MM, HEIGHT_MM = 178, 122
INK = "#20272d"
BLUE = "#346B8A"
TEAL = "#549A94"
AMBER = "#CE9B63"
PLUM = "#9385AC"
GRAY = "#84939E"
RULE = "#CBD4D9"
COLORS = {"novel_in_known_gene": TEAL, "new_locus": AMBER,
          "ambiguous_gene": PLUM, "known": BLUE}
LEFT_ORDER = ["novel_in_known_gene", "new_locus", "ambiguous_gene"]
RIGHT_ORDER = ["ambiguous_gene", "novel_in_known_gene", "known", "new_locus"]
NAMES = {"novel_in_known_gene": "Novel in gene", "new_locus": "New locus",
         "ambiguous_gene": "Ambiguous", "known": "Known"}


def render(data, output):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 6.5,
        "text.color": INK, "axes.edgecolor": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK, "axes.linewidth": .5,
        "svg.fonttype": "none", "svg.hashsalt": "txid-editorial-v4",
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    fig = plt.figure(figsize=(WIDTH_MM / 25.4, HEIGHT_MM / 25.4), dpi=150)
    page = fig.add_axes([0, 0, 1, 1], xlim=(0, WIDTH_MM), ylim=(0, HEIGHT_MM))
    page.set_axis_off()
    texts = []

    def text(x, y, value, size=6.5, ha="left", va="center", gid=None, **kw):
        artist = page.text(x, y, value, fontsize=size, ha=ha, va=va,
                           color=INK, linespacing=1.1, **kw)
        if gid: artist.set_gid(gid)
        texts.append(artist)
        return artist

    def line(x1, y1, x2, y2, color=RULE, lw=.5, **kw):
        page.plot([x1, x2], [y1, y2], color=color, lw=lw, **kw)

    def box(x, y, w, h, face="white", edge=RULE, radius=.8):
        page.add_patch(FancyBboxPatch((x, y), w, h,
                      boxstyle=f"round,pad=0,rounding_size={radius}",
                      facecolor=face, edgecolor=edge, linewidth=.6))

    def arrow(x1, y1, x2, y2, dashed=False):
        page.annotate("", xy=(x2, y2), xytext=(x1, y1),
                      arrowprops={"arrowstyle": "-|>", "mutation_scale": 6,
                                  "lw": .6, "color": GRAY,
                                  "linestyle": (0, (2, 2)) if dashed else "solid"})

    def panel(letter, x, y, title):
        text(x, y, letter, 8, weight="bold", gid="panel-" + letter)
        text(x + 5, y, title, 7, weight="bold")

    # a: reference sequence identity and annotation context use separate paths.
    panel("a", 4, 117, "Reference-aware structural identity")
    text(8, 107, "Transcript models", 7)
    text(8, 102.5, "GTF / GFF3", 6)
    xpos = lambda position: 17 + (position - 90) * 30 / 330
    for key, name, y in [("A", "M1", 95), ("B", "M2", 88), ("C", "M3", 81)]:
        exons = data["schematic"][key]
        text(6, y, name, 6)
        line(xpos(exons[0][1]), y, xpos(exons[1][0]), y, GRAY, .6)
        for start, end in exons:
            page.add_patch(Rectangle((xpos(start), y - 1), xpos(end) - xpos(start), 2,
                                     facecolor="#AEBAC3", edgecolor="none"))
        if key == "B":
            for endpoint in [exons[0][0], exons[-1][1]]:
                line(xpos(endpoint), y - 1.3, xpos(endpoint), y + 1.3, BLUE, 1)
        if key == "C":
            line(xpos(exons[0][1]), y - 1.3, xpos(exons[0][1]), y + 1.3, AMBER, 1)
            text(xpos(exons[0][1]), 76.5, "200 → 201", 6, ha="center")
    arrow(49, 88, 57, 88)
    box(59, 79, 44, 23, face="#F3F7F9", edge="#B8CCD8")
    text(81, 98, "Canonical structure", 7, ha="center", weight="bold")
    text(81, 92.5, "Reference · contig · strand", 6, ha="center")
    text(81, 87.5, "SC: ordered introns", 6.5, ha="center")
    text(81, 82.8, "TF: introns + exact ends", 6.5, ha="center")
    box(66, 106, 30, 6, face="white")
    text(81, 109, "Reference FASTA", 6.5, ha="center")
    arrow(81, 106, 81, 102)
    arrow(104, 88, 113, 88)
    text(109, 94, "SHA-256", 5.8, ha="center")
    text(127.5, 108.5, "Exact IDs", 7, ha="center")
    text(121, 102.5, "SC", 6.5, ha="center")
    text(135, 102.5, "TF", 6.5, ha="center")
    for y, sc, tf, scfill, tffill in [
        (95, "SCa", "TFa", "#DFEAF0", "#DFEAF0"),
        (88, "SCa", "TFb", "#DFEAF0", "#E3EFEC"),
        (81, "SCb", "TFc", "#F4E8D9", "#F4E8D9"),
    ]:
        for x, label, fill in [(121, sc, scfill), (135, tf, tffill)]:
            box(x - 5, y - 2.2, 10, 4.4, face=fill, edge=fill, radius=.7)
            text(x, y, label, 6.3, ha="center")
    arrow(142, 88, 150, 88)
    box(152, 79, 22, 23, face="#F3F7F9", edge="#B8CCD8")
    text(163, 97, "Registry", 7, ha="center", weight="bold")
    text(163, 91, "Mappings", 6.5, ha="center")
    text(163, 86.5, "Catalog", 6.5, ha="center")
    text(163, 82, "Classifications", 5.8, ha="center")
    box(141, 67, 33, 6, face="#FAF6F0", edge="#DBC9B0")
    text(157.5, 70, "Annotation versions", 6.5, ha="center")
    arrow(163, 73, 163, 79, dashed=True)

    # b: observed transitions only; ribbon thickness is proportional to count.
    panel("b", 4, 58, "Annotation changes, identities persist")
    totals = data["annotation"]
    total = totals["shared_observations"]
    changed = totals["classification_changes"]
    text(5, 52.5, f"Changed observations only: {changed:,} / {total:,} ({100 * changed / total:.3f}%)", 6.1,
         gid="b-classification-summary")
    transitions = totals["classification_transitions"]
    left_total, right_total = Counter(), Counter()
    flows = []
    for key, count in transitions.items():
        old, new = key.split(" -> ")
        left_total[old] += count; right_total[new] += count
        flows.append({"from": old, "to": new, "count": count})
    assert sum(left_total.values()) == sum(right_total.values()) == changed
    scale = 25 / changed

    def node_positions(order, counts, top, gap):
        positions = {}
        for category in order:
            height = counts[category] * scale
            positions[category] = (top - height, top)
            top -= height + gap
        return positions

    left_pos = node_positions(LEFT_ORDER, left_total, 43, 3)
    right_pos = node_positions(RIGHT_ORDER, right_total, 44.5, 3)
    left_cursor = {key: high for key, (low, high) in left_pos.items()}
    right_cursor = {key: high for key, (low, high) in right_pos.items()}
    flows.sort(key=lambda item: (LEFT_ORDER.index(item["from"]), RIGHT_ORDER.index(item["to"])))
    for item in flows:
        old, new, count = item["from"], item["to"], item["count"]
        height = count * scale
        ltop, rtop = left_cursor[old], right_cursor[new]
        lbottom, rbottom = ltop - height, rtop - height
        left_cursor[old], right_cursor[new] = lbottom, rbottom
        vertices = [(26.5, ltop), (47, ltop), (61, rtop), (81, rtop),
                    (81, rbottom), (61, rbottom), (47, lbottom), (26.5, lbottom), (26.5, ltop)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
                 Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
        ribbon = PathPatch(Path(vertices, codes), facecolor=COLORS[old],
                           alpha=.52, edgecolor="none", zorder=1)
        ribbon.set_gid(f'flow-{old}-to-{new}')
        page.add_patch(ribbon)
        item.update({"height_mm": height, "left_interval": [lbottom, ltop],
                     "right_interval": [rbottom, rtop]})
        if count in (14517, 10732):
            text(53.9375, (ltop + lbottom + rtop + rbottom) / 4,
                 f"{count:,}", 6, ha="center", zorder=4,
                 gid=f"b-flow-count-{old}-to-{new}")
    text(25.6, 47.5, "v29", 6.5, ha="center", weight="bold")
    text(81.8, 47.5, "v49", 6.5, ha="center", weight="bold")
    for side, positions, counts, x, label_x, ha in [
        ("left", left_pos, left_total, 25, 23.5, "right"),
        ("right", right_pos, right_total, 81, 84, "left"),
    ]:
        for category, (low, high) in positions.items():
            page.add_patch(Rectangle((x, low), 1.5, high - low,
                                     facecolor=COLORS[category], edgecolor="none", zorder=3))
            label = f'{NAMES[category]}\n{counts[category]:,}' if category in LEFT_ORDER else f'{NAMES[category]} {counts[category]:,}'
            if side == "right" and category == "new_locus": label = f'New locus {counts[category]:,}'
            text(label_x, (low + high) / 2, label, 6, ha=ha,
                 gid=f'b-{side}-{category}')
    text(5, 4.5, f'Exact form IDs changed: {totals["txid_form_changes"]:,}', 6.5,
         gid="b-exact-id-changes")

    # c: group definitions differ; lengths encode counts, not performance ranks.
    panel("c", 115, 58, "Cohort matrix")
    text(125, 50, "6 sample rows", 6.5)
    ax = fig.add_axes([124 / WIDTH_MM, 17 / HEIGHT_MM, 49 / WIDTH_MM, 29 / HEIGHT_MM])
    ax.set_facecolor("none")
    matrix_by_layer = {row["Layer"]: row for row in data["matrix"]}
    layer_names = ["TxID exact form", "isoSeQL exact ends", "gffcompare tracking", "isoSeQL common junction"]
    matrix = [{"layer": layer, "columns": matrix_by_layer[layer]["Columns"], "rows": matrix_by_layer[layer]["Rows"]} for layer in layer_names]
    ax.bar(range(4), [row["columns"] for row in matrix], width=.55,
           color=[BLUE, "#A8B5BE", "#A8B5BE", "#A8B5BE"], zorder=2)
    ax.set_xlim(-.6, 3.6); ax.set_ylim(0, 550000)
    ax.set_yticks([0, 250000, 500000], ["0", "250k", "500k"])
    ax.set_xticks(range(4), ["TxID\nexact form", "isoSeQL\nexact ends", "gffcompare\ntracking", "isoSeQL\ncommon\njunction"])
    ax.tick_params(axis="both", labelsize=5.8, width=.5, length=2, pad=2)
    ax.tick_params(axis="x", length=0, pad=4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_ylabel("Matrix columns (n)", fontsize=6.5, labelpad=3)
    for i, row in enumerate(matrix):
        x = 124 + (i + .6) / 4.2 * 49
        y = 17 + row["columns"] / 550000 * 29 + 1.8
        text(x, y, f'{row["columns"]:,}', 5.8, ha="center", gid=f'c-columns-{i}')

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    visible = texts + ax.get_xticklabels() + ax.get_yticklabels() + [ax.yaxis.label]
    outside, overlaps = [], []
    for item in visible:
        bb = item.get_window_extent(renderer)
        if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fig.bbox.x1 or bb.y1 > fig.bbox.y1:
            outside.append(item.get_text())
    for i, first in enumerate(visible):
        a = first.get_window_extent(renderer)
        for second in visible[i + 1:]:
            b = second.get_window_extent(renderer)
            if min(a.x1, b.x1) - max(a.x0, b.x0) > .5 and min(a.y1, b.y1) - max(a.y0, b.y0) > .5:
                overlaps.append([first.get_text(), second.get_text()])
    if outside or overlaps: raise ValueError({"outside": outside, "overlaps": overlaps})
    output.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "svg", "pdf", "tif"]:
        options = {}
        if suffix == "svg": options["metadata"] = {"Date": None}
        if suffix == "pdf": options["metadata"] = {"CreationDate": None, "ModDate": None}
        if suffix == "png": options["metadata"] = {"Software": "TxID editorial figure renderer"}
        if suffix == "tif": options["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(output / f"figure1.{suffix}", dpi=1200 if suffix == "tif" else 600, **options)
    annotation = [{"metric": metric, "changed": totals[key], "total": total, "percent": 100 * totals[key] / total}
                  for metric, key in [("Exact form ID", "txid_form_changes"), ("Classification", "classification_changes")]]
    plotted = {"annotation": annotation, "flows": flows, "left_totals": dict(left_total),
               "right_totals": dict(right_total), "flow_scale_mm_per_observation": scale,
               "matrix": matrix, "order_comparison_location": "Supplementary Table S3"}
    (output / "figure1-plot-data.json").write_text(json.dumps(plotted, indent=2) + "\n")
    (output / "figure1-layout-check.json").write_text(json.dumps({
        "style": "workflow and two linked results", "width_mm": WIDTH_MM, "height_mm": HEIGHT_MM,
        "png_dpi": 600, "tiff_dpi": 1200, "panels": ["a", "b", "c"],
        "minimum_font_pt": min(t.get_fontsize() for t in visible),
        "maximum_body_font_pt": max(t.get_fontsize() for t in visible if not (t.get_gid() or "").startswith("panel-")),
        "panel_font_pt": 8, "text_elements_checked": len(visible),
        "out_of_bounds_text": outside, "overlapping_text": overlaps,
        "all_text_black": all(matplotlib.colors.to_hex(t.get_color()) == INK for t in visible),
        "background_gridlines": False, "quantitative_y_limits": [list(ax.get_ylim())],
    }, indent=2) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    package = FilePath(__file__).resolve().parents[1]
    render(json.loads((package / "figures/source-data.json").read_text()), package / "figures")
