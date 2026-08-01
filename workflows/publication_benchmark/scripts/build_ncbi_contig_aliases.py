#!/usr/bin/env python3
"""Build a validated TxID contig-alias table from an NCBI assembly report."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path


def _file_record(path: Path) -> dict[str, object]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(block)
            digest.update(block)
    return {"path": str(path.resolve()), "bytes": size, "sha256": digest.hexdigest()}


def _read_fai(path: Path) -> dict[str, int]:
    contigs: dict[str, int] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: malformed FASTA index row")
            name = fields[0]
            length = int(fields[1])
            if not name or length <= 0:
                raise ValueError(f"{path}:{line_number}: invalid contig name or length")
            if name in contigs:
                raise ValueError(f"{path}:{line_number}: duplicate contig {name!r}")
            contigs[name] = length
    if not contigs:
        raise ValueError(f"{path}: FASTA index is empty")
    return contigs


def build(
    assembly_report: Path,
    fasta_fai: Path,
    output: Path,
) -> dict[str, object]:
    fasta_contigs = _read_fai(fasta_fai)
    aliases: dict[str, str] = {}
    matched_targets: set[str] = set()
    report_rows = 0
    metadata: dict[str, str] = {}
    with assembly_report.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                if ":" in line:
                    key, value = line[1:].split(":", 1)
                    metadata[key.strip()] = value.strip()
                continue
            fields = line.split("\t")
            if len(fields) < 10:
                raise ValueError(
                    f"{assembly_report}:{line_number}: expected at least 10 columns"
                )
            report_rows += 1
            sequence_name, genbank, refseq = fields[0], fields[4], fields[6]
            sequence_length = int(fields[8])
            ucsc_name = fields[9]
            if ucsc_name == "na" or ucsc_name not in fasta_contigs:
                continue
            if fasta_contigs[ucsc_name] != sequence_length:
                raise ValueError(
                    f"{assembly_report}:{line_number}: {ucsc_name!r} length "
                    f"{sequence_length} conflicts with FASTA length "
                    f"{fasta_contigs[ucsc_name]}"
                )
            matched_targets.add(ucsc_name)
            for alias in (sequence_name, genbank, refseq):
                if alias in {"", "na", ucsc_name}:
                    continue
                previous = aliases.get(alias)
                if previous is not None and previous != ucsc_name:
                    raise ValueError(
                        f"assembly report maps alias {alias!r} to both "
                        f"{previous!r} and {ucsc_name!r}"
                    )
                aliases[alias] = ucsc_name
    if not matched_targets:
        raise ValueError("assembly report has no UCSC targets present in the FASTA index")

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        handle.write("alias\tprimary\n")
        for alias, primary in sorted(aliases.items()):
            handle.write(f"{alias}\t{primary}\n")
    os.replace(temporary, output)
    return {
        "schema": "txid.ncbi-assembly-report-aliases.v1",
        "assembly_report_metadata": dict(sorted(metadata.items())),
        "report_rows": report_rows,
        "fasta_contigs": len(fasta_contigs),
        "report_targets_present_in_fasta": len(matched_targets),
        "fasta_contigs_without_report_target": sorted(
            set(fasta_contigs) - matched_targets
        ),
        "aliases": len(aliases),
        "validation": (
            "Every emitted alias is backed by an NCBI assembly-report row whose "
            "UCSC-style target exists in the FASTA index with the declared length."
        ),
        "inputs": [_file_record(assembly_report), _file_record(fasta_fai)],
        "output": _file_record(output),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--fasta-fai", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    args = parser.parse_args()
    started = time.time()
    record = build(args.assembly_report, args.fasta_fai, args.output)
    record["started_unix"] = started
    record["ended_unix"] = time.time()
    args.audit.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.audit.with_name(f".{args.audit.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, args.audit)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
