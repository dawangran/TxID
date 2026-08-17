"""Dependency-free SVG plots for a TxID registry."""

from __future__ import annotations

import hashlib
import html
import json
from pathlib import Path
from typing import Any, Iterable

from .registry import Registry
from .writers import _atomic_text


def _escape(value: object) -> str:
    return html.escape(str(value), quote=True)


def write_registry_svg(path: str | Path, registry: Registry) -> None:
    summary = registry.summary()
    classifications = summary["classifications"]
    assert isinstance(classifications, dict)
    labels = ["known", "novel_in_known_gene", "ambiguous_gene", "new_locus"]
    counts = [int(classifications.get(label, 0)) for label in labels]
    maximum = max(counts, default=1) or 1
    width = 820
    height = 170 + 58 * len(labels)

    def emit(handle):
        handle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        handle.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">\n'
        )
        handle.write('<title id="title">TxID registry overview</title>\n')
        handle.write(
            f'<desc id="desc">Classification counts for {_escape(summary["assembly"])}</desc>\n'
        )
        handle.write('<rect width="100%" height="100%" fill="#f8fafc"/>\n')
        handle.write(
            '<text x="32" y="42" font-family="sans-serif" font-size="25" '
            'font-weight="700" fill="#172554">TxID registry overview</text>\n'
        )
        handle.write(
            f'<text x="32" y="72" font-family="sans-serif" font-size="14" fill="#475569">'
            f'{_escape(summary["assembly"])} · {summary["imports"]} imports · '
            f'{summary["exact_forms"]} exact forms</text>\n'
        )
        for index, (label, count) in enumerate(zip(labels, counts)):
            y = 112 + index * 58
            bar_width = 470 * count / maximum
            handle.write(
                f'<text x="32" y="{y + 21}" font-family="monospace" font-size="14" '
                f'fill="#334155">{_escape(label)}</text>\n'
            )
            handle.write(
                f'<rect x="230" y="{y}" width="470" height="28" rx="4" fill="#e2e8f0"/>\n'
            )
            handle.write(
                f'<rect x="230" y="{y}" width="{bar_width:.2f}" height="28" rx="4" fill="#2563eb"/>\n'
            )
            handle.write(
                f'<text x="715" y="{y + 20}" font-family="sans-serif" font-size="14" '
                f'fill="#172554">{count}</text>\n'
            )
        handle.write("</svg>\n")

    _atomic_text(path, emit)


def _exons(value: str) -> tuple[tuple[int, int], ...]:
    return tuple((int(start), int(end)) for start, end in json.loads(value))


def _wrap_items(prefix: str, values: Iterable[str], *, width: int = 58) -> list[str]:
    items = list(values)
    if not items:
        return []
    lines: list[str] = []
    current = prefix
    continuation = " " * len(prefix)
    for item in items:
        separator = "" if current in {prefix, continuation} else " · "
        if separator and len(current) + len(separator) + len(item) > width:
            lines.append(current)
            current = continuation + item
        else:
            current += separator + item
    lines.append(current)
    return lines


def _reference_tracks(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[object, ...], dict[str, Any]] = {}
    for row in rows:
        exons = _exons(str(row["exons_json"]))
        key = (
            row["reference_transcript_id"],
            row["form_id"],
            row["contig"],
            row["strand"],
            exons,
        )
        track = grouped.setdefault(
            key,
            {
                "reference_transcript_id": str(row["reference_transcript_id"]),
                "form_id": str(row["form_id"]),
                "contig": str(row["contig"]),
                "strand": str(row["strand"]),
                "exons": exons,
                "splice_chain_id": (
                    None
                    if row["splice_chain_id"] is None
                    else str(row["splice_chain_id"])
                ),
                "annotations": set(),
            },
        )
        track["annotations"].add(str(row["annotation_name"]))
    return sorted(
        grouped.values(),
        key=lambda row: (
            row["exons"][0][0],
            row["exons"][-1][1],
            row["exons"],
            row["form_id"],
            row["reference_transcript_id"],
        ),
    )


def _observed_tracks(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[str, dict[str, Any]] = {}
    for row in rows:
        form_id = str(row["form_id"])
        exons = _exons(str(row["exons_json"]))
        track = grouped.setdefault(
            form_id,
            {
                "form_id": form_id,
                "contig": str(row["contig"]),
                "strand": str(row["strand"]),
                "exons": exons,
                "splice_chain_ids": set(),
                "old_ids": set(),
                "output_ids": set(),
                "classifications": set(),
                "annotations": set(),
                "samples": set(),
                "tools": set(),
                "fuzzy": set(),
                "observations": 0,
            },
        )
        if track["exons"] != exons:
            raise ValueError(f"exact form {form_id} has inconsistent stored exon models")
        if row["splice_chain_id"] is not None:
            track["splice_chain_ids"].add(str(row["splice_chain_id"]))
        track["old_ids"].add(
            (str(row["sample"]), str(row["upstream_tool"]), str(row["original_transcript_id"]))
        )
        track["output_ids"].add(str(row["output_transcript_id"]))
        track["classifications"].add(str(row["classification"]))
        track["annotations"].add(str(row["annotation_name"]))
        track["samples"].add(str(row["sample"]))
        track["tools"].add(str(row["upstream_tool"]))
        if row["fuzzy_cluster_id"] is not None:
            track["fuzzy"].add(
                (
                    str(row["fuzzy_cluster_id"]),
                    str(row["bridge_status"]),
                    str(row["fuzzy_algorithm"]),
                    int(row["splice_tolerance"]),
                    int(row["end_tolerance"]),
                )
            )
        track["observations"] += 1
    return sorted(
        grouped.values(),
        key=lambda row: (
            row["exons"][0][0],
            row["exons"][-1][1],
            row["exons"],
            row["form_id"],
        ),
    )


def _cluster_color(cluster_id: str | None) -> str:
    if cluster_id is None:
        return "#2563eb"
    hue = int(hashlib.sha256(cluster_id.encode("utf-8")).hexdigest()[:8], 16) % 360
    return f"hsl({hue} 58% 42%)"


def _track_lines(track: dict[str, Any], *, reference: bool) -> tuple[list[str], list[str]]:
    if reference:
        left = [f'TxID (exact): {track["form_id"]}']
        left += _wrap_items("Reference ID: ", [track["reference_transcript_id"]])
        right = _wrap_items("Annotation: ", sorted(track["annotations"]))
        right += _wrap_items(
            "Splice-chain TxID: ",
            [] if track["splice_chain_id"] is None else [track["splice_chain_id"]],
        )
        return left, right

    old_ids = [
        f"{sample}/{tool}:{transcript_id}"
        for sample, tool, transcript_id in sorted(track["old_ids"])
    ]
    left = [f'TxID (new): {track["form_id"]}']
    left += _wrap_items("Upstream (old): ", old_ids)

    fuzzy = sorted(track["fuzzy"])
    if fuzzy:
        fuzzy_ids = [item[0] for item in fuzzy]
        right = _wrap_items("Fuzzy: ", fuzzy_ids)
        if any(item[1] == "ambiguous_bridge" for item in fuzzy):
            right.append("Warning: ambiguous bridge")
    else:
        right = ["Fuzzy: unassigned (exact only)"]
    right.append(
        f'{track["observations"]} observations · {len(track["samples"])} samples · '
        f'{len(track["tools"])} tools'
    )
    right += _wrap_items("Class: ", sorted(track["classifications"]))
    right += _wrap_items("Splice-chain TxID: ", sorted(track["splice_chain_ids"]))
    different_outputs = sorted(
        output_id for output_id in track["output_ids"] if output_id != track["form_id"]
    )
    right += _wrap_items("Output/reference ID: ", different_outputs)
    return left, right


def _draw_text_lines(
    handle: Any,
    lines: list[str],
    *,
    x: int,
    y: int,
    fill: str,
    font_size: int = 12,
) -> None:
    if not lines:
        return
    handle.write(
        f'<text x="{x}" y="{y}" font-family="ui-monospace,monospace" '
        f'font-size="{font_size}" fill="{fill}">\n'
    )
    for index, line in enumerate(lines):
        dy = 0 if index == 0 else 17
        handle.write(f'<tspan x="{x}" dy="{dy}">{_escape(line)}</tspan>\n')
    handle.write("</text>\n")


def _draw_model(
    handle: Any,
    exons: tuple[tuple[int, int], ...],
    strand: str,
    *,
    center_y: float,
    plot_x: int,
    plot_width: int,
    coordinate_min: int,
    coordinate_max: int,
    color: str,
) -> None:
    denominator = coordinate_max + 1 - coordinate_min

    def x(coordinate: int) -> float:
        return plot_x + plot_width * (coordinate - coordinate_min) / denominator

    start = exons[0][0]
    end = exons[-1][1]
    handle.write(
        f'<line x1="{x(start):.2f}" y1="{center_y:.2f}" x2="{x(end + 1):.2f}" '
        f'y2="{center_y:.2f}" stroke="{color}" stroke-width="2"/>\n'
    )
    for exon_start, exon_end in exons:
        exon_x = x(exon_start)
        exon_width = max(2.0, x(exon_end + 1) - exon_x)
        handle.write(
            f'<rect x="{exon_x:.2f}" y="{center_y - 8:.2f}" width="{exon_width:.2f}" '
            f'height="16" rx="1.5" fill="{color}"/>\n'
        )
    tip_x = x(end + 1) if strand == "+" else x(start)
    direction = 1 if strand == "+" else -1
    points = (
        f"{tip_x:.2f},{center_y:.2f} "
        f"{tip_x - direction * 8:.2f},{center_y - 5:.2f} "
        f"{tip_x - direction * 8:.2f},{center_y + 5:.2f}"
    )
    handle.write(f'<polygon points="{points}" fill="{color}"/>\n')


def write_gene_svg(path: str | Path, registry: Registry, gene: str) -> dict[str, Any]:
    """Write an exact-form gene structure plot with optional fuzzy membership."""

    data = registry.gene_plot_data(gene)
    references = _reference_tracks(data["references"])
    observed = _observed_tracks(data["observations"])
    all_tracks = [*references, *observed]
    all_exons = [exon for track in all_tracks for exon in track["exons"]]
    transcript_min = min(start for start, _ in all_exons)
    transcript_max = max(end for _, end in all_exons)
    span = transcript_max - transcript_min + 1
    padding = max(5, round(span * 0.03))
    coordinate_min = max(1, transcript_min - padding)
    coordinate_max = transcript_max + padding

    reference_layouts = [
        (*_track_lines(track, reference=True), track) for track in references
    ]
    observed_layouts = [
        (*_track_lines(track, reference=False), track) for track in observed
    ]

    def row_height(left: list[str], right: list[str]) -> int:
        return max(64, 24 + 17 * max(len(left), len(right)))

    header_height = 192
    section_height = 37
    reference_height = sum(row_height(left, right) for left, right, _ in reference_layouts)
    observed_height = sum(row_height(left, right) for left, right, _ in observed_layouts)
    height = (
        header_height
        + (section_height + reference_height if references else 0)
        + section_height
        + observed_height
        + 36
    )
    width = 1800
    plot_x = 480
    plot_width = 690
    right_x = 1205

    fuzzy_specs = sorted(
        {
            (item[2], item[3], item[4])
            for track in observed
            for item in track["fuzzy"]
        }
    )
    if fuzzy_specs:
        fuzzy_summary = "Fuzzy: " + " · ".join(
            f"{algorithm}, splice ±{splice} bp, ends ±{ends} bp"
            for algorithm, splice, ends in fuzzy_specs
        )
    else:
        fuzzy_summary = "Fuzzy: no cluster assignment; exact forms are shown independently"

    def emit(handle: Any) -> None:
        handle.write('<?xml version="1.0" encoding="UTF-8"?>\n')
        handle.write(
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
            f'viewBox="0 0 {width} {height}" role="img" aria-labelledby="title desc">\n'
        )
        handle.write(
            f'<title id="title">Transcript structures for gene {_escape(data["query"])}</title>\n'
        )
        handle.write(
            f'<desc id="desc">Reference and observed exact transcript forms on '
            f'{_escape(data["assembly"])} {_escape(data["contig"])}; fuzzy cluster colors '
            'are diagnostic groupings and do not replace exact TxIDs.</desc>\n'
        )
        handle.write('<rect width="100%" height="100%" fill="#f8fafc"/>\n')
        handle.write(
            f'<text x="32" y="42" font-family="sans-serif" font-size="25" '
            f'font-weight="700" fill="#0f172a">Gene transcript structure · '
            f'{_escape(data["query"])}</text>\n'
        )
        resolved = ", ".join(data["gene_ids"])
        handle.write(
            f'<text x="32" y="70" font-family="sans-serif" font-size="14" fill="#475569">'
            f'Resolved by {_escape(data["resolved_by"])} to {_escape(resolved)} · '
            f'{_escape(data["assembly"])} · {_escape(data["contig"])}:'
            f'{transcript_min:,}-{transcript_max:,} ({_escape(data["strand"])})</text>\n'
        )
        handle.write(
            f'<text x="32" y="94" font-family="sans-serif" font-size="14" fill="#475569">'
            f'{len(references)} reference tracks · {len(observed)} observed exact forms · '
            f'{len(data["observations"])} observations</text>\n'
        )
        handle.write(
            f'<text x="32" y="118" font-family="sans-serif" font-size="14" fill="#334155">'
            f'{_escape(fuzzy_summary)}</text>\n'
        )
        handle.write(
            '<rect x="32" y="137" width="18" height="12" rx="1.5" fill="#64748b"/>'
            '<text x="58" y="148" font-family="sans-serif" font-size="12" fill="#475569">'
            'reference exon</text>\n'
        )
        handle.write(
            '<rect x="178" y="137" width="18" height="12" rx="1.5" fill="#2563eb"/>'
            '<text x="204" y="148" font-family="sans-serif" font-size="12" fill="#475569">'
            'observed exon; equal color means equal fuzzy cluster</text>\n'
        )
        handle.write(
            '<text x="32" y="171" font-family="sans-serif" font-size="12" fill="#64748b">'
            'Each observed row remains a distinct exact form. Compare exon boundaries across '
            'same- and different-color rows to assess possible over- or under-grouping.</text>\n'
        )

        axis_y = 184
        handle.write(
            f'<line x1="{plot_x}" y1="{axis_y}" x2="{plot_x + plot_width}" y2="{axis_y}" '
            'stroke="#94a3b8" stroke-width="1"/>\n'
        )
        coordinate_span = coordinate_max - coordinate_min
        for index in range(6):
            fraction = index / 5
            tick_x = plot_x + plot_width * fraction
            coordinate = round(coordinate_min + coordinate_span * fraction)
            handle.write(
                f'<line x1="{tick_x:.2f}" y1="{axis_y - 4}" x2="{tick_x:.2f}" '
                f'y2="{axis_y + 4}" stroke="#94a3b8"/>\n'
            )
            handle.write(
                f'<text x="{tick_x:.2f}" y="{axis_y - 8}" text-anchor="middle" '
                f'font-family="sans-serif" font-size="11" fill="#64748b">{coordinate:,}</text>\n'
            )

        y = header_height

        def section(label: str, count: int) -> None:
            nonlocal y
            handle.write(
                f'<rect x="20" y="{y}" width="1760" height="29" rx="5" fill="#e2e8f0"/>\n'
            )
            handle.write(
                f'<text x="32" y="{y + 20}" font-family="sans-serif" font-size="14" '
                f'font-weight="700" fill="#334155">{_escape(label)} ({count})</text>\n'
            )
            y += section_height

        def row(
            left: list[str],
            right: list[str],
            track: dict[str, Any],
            *,
            reference: bool,
            index: int,
        ) -> None:
            nonlocal y
            current_height = row_height(left, right)
            background = "#ffffff" if index % 2 == 0 else "#f1f5f9"
            fuzzy = sorted(track.get("fuzzy", ()))
            cluster_id = None if not fuzzy else fuzzy[0][0]
            ambiguous = any(item[1] == "ambiguous_bridge" for item in fuzzy)
            color = "#64748b" if reference else _cluster_color(cluster_id)
            css_class = "reference-track" if reference else "observed-track"
            handle.write(
                f'<g class="{css_class}" data-form-id="{_escape(track["form_id"])}" '
                f'data-fuzzy-cluster="{_escape(cluster_id or "")}">\n'
            )
            handle.write(
                f'<rect x="20" y="{y}" width="1760" height="{current_height - 4}" rx="4" '
                f'fill="{background}"'
                + (' stroke="#b91c1c" stroke-width="1.5" stroke-dasharray="5 4"' if ambiguous else "")
                + '/>\n'
            )
            _draw_text_lines(handle, left, x=32, y=y + 20, fill="#0f172a")
            _draw_model(
                handle,
                track["exons"],
                track["strand"],
                center_y=y + 27,
                plot_x=plot_x,
                plot_width=plot_width,
                coordinate_min=coordinate_min,
                coordinate_max=coordinate_max,
                color=color,
            )
            _draw_text_lines(handle, right, x=right_x, y=y + 20, fill="#334155")
            handle.write("</g>\n")
            y += current_height

        if reference_layouts:
            section("Reference transcripts", len(reference_layouts))
            for index, (left, right, track) in enumerate(reference_layouts):
                row(left, right, track, reference=True, index=index)
        section("Observed exact forms", len(observed_layouts))
        for index, (left, right, track) in enumerate(observed_layouts):
            row(left, right, track, reference=False, index=index)
        handle.write("</svg>\n")

    _atomic_text(path, emit)
    return {
        "query": data["query"],
        "resolved_by": data["resolved_by"],
        "gene_ids": data["gene_ids"],
        "reference_tracks": len(references),
        "exact_forms": len(observed),
        "observations": len(data["observations"]),
    }
