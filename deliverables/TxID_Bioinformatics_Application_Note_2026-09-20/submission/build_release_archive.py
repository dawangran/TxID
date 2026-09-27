#!/usr/bin/env python3
"""Build a small, deterministic local source snapshot; never publish or upload."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import tempfile
import zipfile

VERSION = "0.1.3"
STEM = f"TxID-{VERSION}-submission-source"
PACKAGE = "deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/submission"
MAX_FILE_BYTES = 2 * 1024 * 1024
MAX_TOTAL_BYTES = 20 * 1024 * 1024
FIXED_TIME = (1980, 1, 1, 0, 0, 0)

# Every selected file is text. Large data/results and authoring history are absent.
REQUIRED_FILES = (
    ".gitignore", "LICENSE", "README.md", "CITATION.cff", "pyproject.toml",
    "Dockerfile", "conda-recipe/meta.yaml", ".github/workflows/tests.yml",
    "docs/cli.md", "docs/tutorial.md", "docs/reproducibility.md",
    "docs/benchmark.md", "docs/independent-identity.md",
    "docs/txid-batch-wdl-guide.zh-CN.md", "docs/txid-multi-add-wdl-guide.zh-CN.md",
    "docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md",
    "tests/__init__.py", "tests/helpers.py", "tests/conformance/exact-v1.json",
    "benchmarks/run_simulation.py", "benchmarks/prepare_external_comparison.py",
    "benchmarks/evaluate_external_comparison.py", "benchmarks/run_gffcompare_comparison.py",
    "benchmarks/run_isoseql_comparison.py", "benchmarks/independent_identity.py",
    "benchmarks/run_real_interop_callers.py", "benchmarks/evaluate_real_interop.py",
    f"{PACKAGE}/build_release_archive.py", f"{PACKAGE}/test_release_archive.py",
    f"{PACKAGE}/verify_release_archive.py",
)
TEST_MODULES = (
    "test_caller_annotation_policy", "test_classify", "test_cli", "test_conformance", "test_determinism",
    "test_external_comparison", "test_fuzzy", "test_gtf_contig_subset",
    "test_identity", "test_independent_identity", "test_isoseql_compat",
    "test_packaging", "test_parser", "test_plotting", "test_publication_benchmark",
    "test_real_interop", "test_simulation", "test_wdl_contract", "test_workflow",
)
PATTERNS = (
    "src/txid/*.py", "schemas/*.json", "schemas/*.sql",
    "workflows/*.wdl", "workflows/*.json",
    "workflows/publication_benchmark/README.md", "workflows/publication_benchmark/Snakefile",
    "workflows/publication_benchmark/scripts/*.py",
    "workflows/publication_benchmark/config/*.yaml",
    "workflows/publication_benchmark/config/*.tsv",
    "workflows/publication_benchmark/envs/*.yaml",
    "tests/fixtures/*.json", "tests/fixtures/*.tsv", "tests/fixtures/*.gtf",
    "tests/fixtures/*.gff3", "tests/fixtures/*.fa", "tests/fixtures/*.fasta",
)
SENSITIVE_PARTS = {".git", ".ssh", ".aws", ".codex", ".agents", "credentials", "secrets"}
BIO_SUFFIXES = {".bam", ".cram", ".sam", ".fastq", ".fq", ".fa", ".fasta", ".fna", ".gtf", ".gff", ".gff3"}


def selected_paths(root: Path) -> list[str]:
    selected = set(REQUIRED_FILES)
    selected.update(f"tests/{name}.py" for name in TEST_MODULES)
    for pattern in PATTERNS:
        selected.update(path.relative_to(root).as_posix() for path in root.glob(pattern))
    missing = [name for name in selected if not (root / name).is_file()]
    if missing:
        raise ValueError("Required release files missing: " + ", ".join(sorted(missing)))
    return sorted(selected)


def validate_path(root: Path, relative: str, allowed: set[str]) -> Path:
    pure = PurePosixPath(relative)
    if pure.is_absolute() or ".." in pure.parts or "\\" in relative or pure.as_posix() != relative:
        raise ValueError(f"Unsafe archive path: {relative}")
    if relative not in allowed:
        raise ValueError(f"Path not in release allowlist: {relative}")
    if any(part.lower() in SENSITIVE_PARTS or part.lower().startswith(".env") for part in pure.parts):
        raise ValueError(f"Sensitive archive path: {relative}")
    if Path(relative).suffix.lower() in {".pem", ".key", ".p12", ".sqlite", ".db", ".gz", ".zip"}:
        raise ValueError(f"Excluded archive file type: {relative}")
    if Path(relative).suffix.lower() in BIO_SUFFIXES and not relative.startswith("tests/fixtures/"):
        raise ValueError(f"Biological data outside explicit tiny fixtures: {relative}")
    path = root
    for part in pure.parts:
        path = path / part
        if path.is_symlink():
            raise ValueError(f"Symlink excluded: {relative}")
    if not path.is_file() or not path.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Not a regular in-tree file: {relative}")
    return path


def collect_files(root: Path, paths: list[str], allowed: set[str]) -> dict[str, bytes]:
    contents: dict[str, bytes] = {}
    total = 0
    for relative in sorted(set(paths)):
        path = validate_path(root, relative, allowed)
        if path.stat().st_size > MAX_FILE_BYTES:
            raise ValueError(f"File exceeds 2 MiB: {relative}")
        with path.open("rb") as handle:
            data = handle.read(MAX_FILE_BYTES + 1)
        if len(data) > MAX_FILE_BYTES:
            raise ValueError(f"File exceeds 2 MiB: {relative}")
        if re.search(rb"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----", data):
            raise ValueError(f"Private-key material excluded: {relative}")
        data.decode("utf-8")
        total += len(data)
        if total > MAX_TOTAL_BYTES:
            raise ValueError("Release inputs exceed 20 MiB")
        contents[relative] = data
    return contents


def render_archive(contents: dict[str, bytes]) -> tuple[bytes, dict]:
    records = [
        {"path": name, "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        for name, data in sorted(contents.items())
    ]
    manifest = {
        "format": "txid-local-source-snapshot-1", "package_version": VERSION,
        "working_tree_snapshot": True, "public_deposit_verified": False,
        "software_archive_doi": None, "manuscript_data_archive_doi": None,
        "included_test_modules": [f"tests.{name}" for name in TEST_MODULES],
        "excluded_scope": ["large primary data", "registries", "historical benchmark results",
                           "manuscript and figures", "site and manuscript-figure tests"],
        "files": records,
    }
    payload = dict(contents)
    payload["RELEASE-MANIFEST.json"] = (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode()
    payload["SHA256SUMS"] = "".join(
        f"{hashlib.sha256(data).hexdigest()}  {name}\n" for name, data in sorted(payload.items())
    ).encode()
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_STORED) as archive:
        for name, data in sorted(payload.items()):
            info = zipfile.ZipInfo(f"{STEM}/{name}", FIXED_TIME)
            info.create_system = 3
            info.external_attr = (0o100644 << 16)
            archive.writestr(info, data)
    return output.getvalue(), manifest


def write_atomic(path: Path, data: bytes) -> None:
    with tempfile.NamedTemporaryFile(dir=path.parent, prefix=f".{path.name}.", delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(data)
    try:
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def build(root: Path, outdir: Path) -> dict:
    paths = selected_paths(root)
    contents = collect_files(root, paths, set(paths))
    project = contents["pyproject.toml"].decode()
    if not re.search(rf'^version\s*=\s*"{re.escape(VERSION)}"\s*$', project, re.MULTILINE):
        raise ValueError("Archive version does not match pyproject.toml")
    archive, manifest = render_archive(contents)
    digest = hashlib.sha256(archive).hexdigest()
    manifest.update({"archive_name": f"{STEM}.zip", "archive_sha256": digest,
                     "archive_bytes": len(archive)})
    outdir.mkdir(parents=True, exist_ok=True)
    write_atomic(outdir / f"{STEM}.zip", archive)
    write_atomic(outdir / f"{STEM}.manifest.json", (json.dumps(manifest, indent=2, sort_keys=True) + "\n").encode())
    write_atomic(outdir / f"{STEM}.zip.sha256", f"{digest}  {STEM}.zip\n".encode())
    return {"archive": str(outdir / f"{STEM}.zip"), "sha256": digest,
            "files": len(contents), "archive_bytes": len(archive), "public_deposit_verified": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[3])
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    try:
        report = build(args.repo.resolve(), args.output_dir.resolve())
    except (OSError, ValueError) as error:
        parser.exit(1, f"release archive: {error}\n")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
