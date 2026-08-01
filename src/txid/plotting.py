"""Dependency-free SVG overview plots for a TxID registry."""

from __future__ import annotations

import html
from pathlib import Path

from .registry import Registry
from .writers import _atomic_text


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
            f'<desc id="desc">Classification counts for {html.escape(str(summary["assembly"]))}</desc>\n'
        )
        handle.write('<rect width="100%" height="100%" fill="#f8fafc"/>\n')
        handle.write(
            f'<text x="32" y="42" font-family="sans-serif" font-size="25" font-weight="700" fill="#172554">TxID registry overview</text>\n'
        )
        handle.write(
            f'<text x="32" y="72" font-family="sans-serif" font-size="14" fill="#475569">'
            f'{html.escape(str(summary["assembly"]))} · {summary["imports"]} imports · '
            f'{summary["exact_forms"]} exact forms</text>\n'
        )
        for index, (label, count) in enumerate(zip(labels, counts)):
            y = 112 + index * 58
            bar_width = 470 * count / maximum
            handle.write(
                f'<text x="32" y="{y + 21}" font-family="monospace" font-size="14" fill="#334155">{html.escape(label)}</text>\n'
            )
            handle.write(
                f'<rect x="230" y="{y}" width="470" height="28" rx="4" fill="#e2e8f0"/>\n'
            )
            handle.write(
                f'<rect x="230" y="{y}" width="{bar_width:.2f}" height="28" rx="4" fill="#2563eb"/>\n'
            )
            handle.write(
                f'<text x="715" y="{y + 20}" font-family="sans-serif" font-size="14" fill="#172554">{count}</text>\n'
            )
        handle.write("</svg>\n")

    _atomic_text(path, emit)

