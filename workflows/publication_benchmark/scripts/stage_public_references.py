#!/usr/bin/env python3
"""Stage checked LRGASP/GENCODE references and build the v49+SIRV context."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
from typing import TextIO

from stage_encode_files import download_and_verify


SOURCE_COLUMNS = {
    "source_id",
    "filename",
    "target_path",
    "source_url",
    "bytes",
    "md5",
    "license",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        rows = list(reader)
        fields = list(reader.fieldnames)
    return fields, rows


def _write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def _open_maybe_gzip(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open(encoding="utf-8")


def build_combined_gtf(
    annotation: Path, sirv_truth: Path, output: Path
) -> dict[str, object]:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    line_counts = {}
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with __import__("io").TextIOWrapper(
                compressed, encoding="utf-8", newline="\n"
            ) as writer:
                for source in (annotation, sirv_truth):
                    count = 0
                    last_had_newline = True
                    with _open_maybe_gzip(source) as handle:
                        for line in handle:
                            if count and not last_had_newline:
                                writer.write("\n")
                            writer.write(line)
                            last_had_newline = line.endswith("\n")
                            count += 1
                    if count and not last_had_newline:
                        writer.write("\n")
                    line_counts[str(source.resolve())] = count
    os.replace(temporary, output)
    return {
        "path": str(output.resolve()),
        "bytes": output.stat().st_size,
        "sha256": _sha256(output),
        "inputs": line_counts,
        "compression": "gzip with mtime=0 and empty original filename",
    }


def stage(
    sources_path: Path,
    references_path: Path,
    output_manifest: Path,
    audit_path: Path,
    *,
    repo: Path,
    download: bool,
    curl_program: str,
) -> dict[str, object]:
    source_fields, sources = _read_tsv(sources_path)
    missing = sorted(SOURCE_COLUMNS - set(source_fields))
    if missing:
        raise ValueError(
            f"{sources_path}: missing columns: {', '.join(missing)}"
        )
    if output_manifest.resolve() == references_path.resolve():
        raise ValueError("refusing to overwrite the source reference manifest")
    records = []
    by_id = {}
    for row in sources:
        target = Path(row["target_path"])
        if not target.is_absolute():
            target = repo / target
        if download:
            file_record, command = download_and_verify(
                target,
                url=row["source_url"],
                size=int(row["bytes"]),
                md5=row["md5"],
                curl_program=curl_program,
            )
        else:
            if not target.is_file():
                raise ValueError(f"staged reference not found: {target}")
            from stage_encode_files import verify_file

            file_record = verify_file(
                target, size=int(row["bytes"]), md5=row["md5"]
            )
            command = None
        record = {
            **row,
            **file_record,
            "downloaded": command is not None,
            "download_argv": command,
        }
        records.append(record)
        by_id[row["source_id"]] = target

    combined_v49 = repo / (
        "benchmark-inputs/lrgasp/reference/gencode-v49-plus-sirv4.gtf.gz"
    )
    combined_record = build_combined_gtf(
        by_id["gencode-v49"],
        by_id["lrgasp-sirv4-truth"],
        combined_v49,
    )

    reference_fields, references = _read_tsv(references_path)
    fasta = by_id["lrgasp-grch38-sirv4"]
    annotations = {
        "gencode-v38": by_id["lrgasp-gencode-v38-sirv4"],
        "gencode-v49": combined_v49,
    }
    for row in references:
        annotation_id = row["annotation_id"]
        if annotation_id not in annotations:
            continue
        if Path(row["fasta"]) != fasta.relative_to(repo):
            raise ValueError(
                f"{references_path}: unexpected FASTA path for {annotation_id}"
            )
        if Path(row["annotation"]) != annotations[annotation_id].relative_to(repo):
            raise ValueError(
                f"{references_path}: unexpected annotation path for {annotation_id}"
            )
        row["fasta_sha256"] = _sha256(fasta)
        row["annotation_sha256"] = _sha256(annotations[annotation_id])
        row["enabled"] = "true"
    _write_tsv(output_manifest, reference_fields, references)
    audit = {
        "schema": "txid.public-reference-staging.v1",
        "source_registry": str(sources_path.resolve()),
        "source_reference_manifest": str(references_path.resolve()),
        "output_reference_manifest": str(output_manifest.resolve()),
        "sources": sorted(records, key=lambda value: value["source_id"]),
        "combined_v49_sirv": combined_record,
        "shared_assembly_fasta_sha256": _sha256(fasta),
    }
    _atomic_json(audit_path, audit)
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output-manifest", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--curl-program", default="curl")
    args = parser.parse_args()
    stage(
        args.sources,
        args.references,
        args.output_manifest,
        args.audit,
        repo=args.repo.resolve(),
        download=args.download,
        curl_program=args.curl_program,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
