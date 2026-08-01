#!/usr/bin/env python3
"""Stage ENCODE files with resumable downloads and audited checksums."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
import time
import urllib.request
from pathlib import Path
from typing import Iterable


ENCODE_API = "https://www.encodeproject.org/files/{accession}/?format=json"
ENCODE_DOWNLOAD = "https://www.encodeproject.org{href}"
TRUE = {"1", "true", "yes"}


def _digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def _read_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"{path}: missing header")
        rows = list(reader)
        fields = list(reader.fieldnames)
    required = {"public_accession", "input_path", "source_url", "sha256", "enabled"}
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError(f"{path}: missing columns: {', '.join(missing)}")
    return fields, rows


def _write_manifest(
    path: Path, fields: list[str], rows: Iterable[dict[str, str]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, path)


def fetch_metadata(accession: str, path: Path) -> dict[str, object]:
    request = urllib.request.Request(
        ENCODE_API.format(accession=accession),
        headers={"Accept": "application/json", "User-Agent": "TxID-benchmark/1"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        metadata = json.load(response)
    _atomic_json(path, metadata)
    return metadata


def load_metadata(
    accession: str, metadata_dir: Path, *, refresh: bool
) -> tuple[dict[str, object], Path]:
    path = metadata_dir / f"{accession}.json"
    if refresh or not path.is_file():
        metadata = fetch_metadata(accession, path)
    else:
        metadata = json.loads(path.read_text(encoding="utf-8"))
    return metadata, path


def validate_metadata(
    accession: str,
    metadata: dict[str, object],
    *,
    expected_file_format: str = "fastq",
) -> dict[str, object]:
    if metadata.get("accession") != accession:
        raise ValueError(
            f"{accession}: metadata accession is {metadata.get('accession')!r}"
        )
    if metadata.get("status") != "released":
        raise ValueError(f"{accession}: ENCODE status is not released")
    if metadata.get("file_format") != expected_file_format:
        raise ValueError(
            f"{accession}: ENCODE file format is "
            f"{metadata.get('file_format')!r}, expected {expected_file_format!r}"
        )
    href = metadata.get("href")
    size = metadata.get("file_size")
    md5sum = metadata.get("md5sum")
    if not isinstance(href, str) or not href.startswith("/files/"):
        raise ValueError(f"{accession}: missing or invalid download href")
    if not isinstance(size, int) or size <= 0:
        raise ValueError(f"{accession}: missing or invalid file size")
    if not isinstance(md5sum, str) or len(md5sum) != 32:
        raise ValueError(f"{accession}: missing or invalid MD5")
    return {
        "accession": accession,
        "status": metadata["status"],
        "file_format": metadata["file_format"],
        "file_size": size,
        "md5": md5sum.lower(),
        "download_url": ENCODE_DOWNLOAD.format(href=href),
        "dataset": metadata.get("dataset"),
        "biological_replicates": metadata.get("biological_replicates", []),
        "technical_replicates": metadata.get("technical_replicates", []),
    }


def verify_file(path: Path, *, size: int, md5: str) -> dict[str, object]:
    observed_size = path.stat().st_size
    if observed_size != size:
        raise ValueError(
            f"{path}: size mismatch: observed {observed_size}, expected {size}"
        )
    observed_md5 = _digest(path, "md5")
    if observed_md5.lower() != md5.lower():
        raise ValueError(
            f"{path}: MD5 mismatch: observed {observed_md5}, expected {md5}"
        )
    return {
        "path": str(path.resolve()),
        "bytes": observed_size,
        "md5": observed_md5,
        "sha256": _digest(path, "sha256"),
    }


def download_and_verify(
    path: Path,
    *,
    url: str,
    size: int,
    md5: str,
    curl_program: str,
) -> tuple[dict[str, object], list[str] | None]:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file():
        return verify_file(path, size=size, md5=md5), None
    partial = path.with_name(path.name + ".part")
    if partial.is_file() and partial.stat().st_size == size:
        record = verify_file(partial, size=size, md5=md5)
        os.replace(partial, path)
        record["path"] = str(path.resolve())
        return record, None
    command = [
        curl_program,
        "--fail",
        "--location",
        "--retry",
        "8",
        "--retry-delay",
        "5",
        "--continue-at",
        "-",
        "--output",
        str(partial),
        url,
    ]
    subprocess.run(command, check=True)
    record = verify_file(partial, size=size, md5=md5)
    os.replace(partial, path)
    record["path"] = str(path.resolve())
    return record, command


def stage(
    manifest: Path,
    metadata_dir: Path,
    audit_path: Path,
    *,
    repo: Path,
    accessions: set[str],
    output_manifest: Path | None,
    download: bool,
    refresh_metadata: bool,
    curl_program: str,
) -> dict[str, object]:
    fields, rows = _read_manifest(manifest)
    available = {
        row["public_accession"]
        for row in rows
        if row["public_accession"].startswith("ENCFF")
    }
    selected = accessions or available
    unknown = sorted(selected - available)
    if unknown:
        raise ValueError(f"accessions absent from manifest: {', '.join(unknown)}")
    if download and output_manifest is None:
        raise ValueError("--output-manifest is required with --download")
    if output_manifest is not None and output_manifest.resolve() == manifest.resolve():
        raise ValueError("refusing to overwrite the source manifest")

    staged = []
    started = time.time()
    for row in rows:
        accession = row["public_accession"]
        if accession not in selected:
            continue
        metadata, metadata_path = load_metadata(
            accession, metadata_dir, refresh=refresh_metadata
        )
        expected_file_format = row.get("expected_file_format", "fastq") or "fastq"
        official = validate_metadata(
            accession,
            metadata,
            expected_file_format=expected_file_format,
        )
        declared_bytes = row.get("expected_bytes", "")
        if declared_bytes not in {"", "NA"}:
            if int(declared_bytes) != int(official["file_size"]):
                raise ValueError(
                    f"{accession}: manifest expected_bytes={declared_bytes}, "
                    f"official metadata reports {official['file_size']}"
                )
        declared_md5 = row.get("expected_md5", "")
        if declared_md5 not in {"", "NA"}:
            if declared_md5.lower() != str(official["md5"]).lower():
                raise ValueError(
                    f"{accession}: manifest expected_md5={declared_md5}, "
                    f"official metadata reports {official['md5']}"
                )
        record: dict[str, object] = {
            **official,
            "dataset_id": row.get("dataset_id"),
            "metadata_path": str(metadata_path.resolve()),
            "metadata_sha256": _digest(metadata_path, "sha256"),
            "downloaded": False,
        }
        if download:
            target = Path(row["input_path"])
            if not target.is_absolute():
                target = repo / target
            file_record, command = download_and_verify(
                target,
                url=str(official["download_url"]),
                size=int(official["file_size"]),
                md5=str(official["md5"]),
                curl_program=curl_program,
            )
            record.update(file_record)
            record["downloaded"] = command is not None
            record["download_argv"] = command
            row["input_kind"] = "local"
            row["source_url"] = str(official["download_url"])
            row["sha256"] = str(file_record["sha256"])
            row["enabled"] = "true"
        staged.append(record)

    if len(staged) != len(selected):
        raise ValueError("not every selected accession was staged")
    if output_manifest is not None:
        _write_manifest(output_manifest, fields, rows)
    audit = {
        "schema": "txid.encode-staging.v2",
        "source_manifest": str(manifest.resolve()),
        "output_manifest": (
            str(output_manifest.resolve()) if output_manifest is not None else None
        ),
        "metadata_only": not download,
        "started_unix": started,
        "ended_unix": time.time(),
        "files": sorted(staged, key=lambda item: str(item["accession"])),
    }
    _atomic_json(audit_path, audit)
    return audit


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--metadata-dir", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--accession", action="append", default=[])
    parser.add_argument("--output-manifest", type=Path)
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--refresh-metadata", action="store_true")
    parser.add_argument("--curl-program", default="curl")
    args = parser.parse_args()
    stage(
        args.manifest,
        args.metadata_dir,
        args.audit,
        repo=args.repo.resolve(),
        accessions=set(args.accession),
        output_manifest=args.output_manifest,
        download=args.download,
        refresh_metadata=args.refresh_metadata,
        curl_program=args.curl_program,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
