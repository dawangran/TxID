#!/usr/bin/env python3
"""Build one TxID registry in a declared order and retain incremental snapshots."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import random
import sqlite3
import subprocess
import sys
import time
from contextlib import closing
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


def _parse_inputs(values: list[str]) -> list[tuple[str, str, Path]]:
    rows = []
    for value in values:
        fields = value.split("=", 2)
        if len(fields) != 3 or not all(fields):
            raise ValueError(f"invalid --input value {value!r}; expected caller=dataset=GTF")
        rows.append((fields[0], fields[1], Path(fields[2])))
    return rows


def _ordered(
    rows: list[tuple[str, str, Path]], order: str
) -> tuple[list[tuple[str, str, Path]], int | None]:
    rows.sort(key=lambda row: (row[1], row[0], str(row[2])))
    if order == "sorted":
        return rows, None
    prefix = "shuffled-"
    if not order.startswith(prefix):
        raise ValueError(f"unsupported order: {order}")
    seed = int(order[len(prefix) :])
    random.Random(seed).shuffle(rows)
    return rows, seed


def _compress(
    path: Path, *, log_prefix: Path
) -> tuple[Path, dict[str, object]]:
    compressed = Path(f"{path}.gz")
    if compressed.exists():
        raise FileExistsError(f"refusing to overwrite compressed artifact: {compressed}")
    record = run_command(
        ["gzip", "-n", "-6", path],
        log_prefix=log_prefix,
    )
    if not compressed.is_file():
        raise FileNotFoundError(f"gzip did not emit expected artifact: {compressed}")
    return compressed, record


def _read_mapping_rows(path: Path) -> list[dict[str, str]]:
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open("r", encoding="utf-8", newline="")
    with handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _registered_imports(database: Path) -> int:
    uri = f"file:{database.resolve()}?mode=ro"
    with closing(sqlite3.connect(uri, uri=True)) as connection:
        row = connection.execute("SELECT COUNT(*) FROM import_manifest").fetchone()
    if row is None:
        raise RuntimeError("TxID registry did not return an import count")
    return int(row[0])


def _write_catalog_from_observations(
    path: Path, rows: list[dict[str, str]]
) -> dict[str, object]:
    """Write the CLI-equivalent catalog from checked component mappings."""

    started = time.perf_counter()
    grouped: dict[
        tuple[str, str, str],
        tuple[int, set[str], set[str], set[str]],
    ] = {}
    for row in rows:
        key = (row["txid_form"], row["txid_sc"], row["txid_gene_id"])
        if key not in grouped:
            grouped[key] = (0, set(), set(), set())
        count, samples, tools, classifications = grouped[key]
        samples.add(row["sample"])
        tools.add(row["tool"])
        classifications.add(row["classification"])
        grouped[key] = (count + 1, samples, tools, classifications)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    if temporary.exists():
        raise FileExistsError(f"refusing to overwrite temporary catalog: {temporary}")
    columns = [
        "form_id",
        "splice_chain_id",
        "gene_id",
        "classifications",
        "observation_count",
        "sample_count",
        "tool_count",
    ]
    try:
        with temporary.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=columns,
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for (form_id, splice_chain_id, gene_id), values in sorted(
                grouped.items(),
                key=lambda item: (item[0][2], item[0][0]),
            ):
                count, samples, tools, classifications = values
                writer.writerow(
                    {
                        "form_id": form_id,
                        "splice_chain_id": splice_chain_id,
                        "gene_id": gene_id,
                        "classifications": ",".join(sorted(classifications)),
                        "observation_count": count,
                        "sample_count": len(samples),
                        "tool_count": len(tools),
                    }
                )
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise
    return {
        "method": "component_mapping_aggregation_v1",
        "observations": len(rows),
        "catalog_rows": len(grouped),
        "output": str(path.resolve()),
        "wall_seconds": time.perf_counter() - started,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--contig-aliases", type=Path)
    parser.add_argument("--annotation", type=Path, required=True)
    parser.add_argument("--annotation-name", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--input", action="append", required=True)
    parser.add_argument("--order", required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--snapshot-dir", type=Path, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument(
        "--snapshot-policy",
        choices=("all", "none"),
        default="all",
        help=(
            "Write a catalog after every import for the incremental run, or "
            "skip intermediate catalogs for order/release comparison runs."
        ),
    )
    parser.add_argument(
        "--compress-artifacts",
        action="store_true",
        help=(
            "After validation, deterministically gzip rewritten GTFs, component "
            "mapping tables, and the SQLite registry."
        ),
    )
    parser.add_argument(
        "--resume-completed",
        type=int,
        default=0,
        help=(
            "Resume an interrupted run after this many completed ordered imports. "
            "The existing registry, component artifacts, and --resume-audit must "
            "agree; zero starts a fresh registry."
        ),
    )
    parser.add_argument(
        "--resume-audit",
        type=Path,
        help="Checked recovery/validation audit required with --resume-completed.",
    )
    args = parser.parse_args()
    args.database.parent.mkdir(parents=True, exist_ok=True)
    args.snapshot_dir.mkdir(parents=True, exist_ok=True)
    args.artifact_dir.mkdir(parents=True, exist_ok=True)
    rows, seed = _ordered(_parse_inputs(args.input), args.order)
    if not 0 <= args.resume_completed <= len(rows):
        raise ValueError(
            f"--resume-completed must be between 0 and {len(rows)}, "
            f"not {args.resume_completed}"
        )
    resume_database_record = None
    if args.resume_completed:
        if not args.database.is_file():
            raise FileNotFoundError(
                f"resume registry does not exist: {args.database}"
            )
        journal = Path(f"{args.database}-journal")
        if journal.exists():
            raise RuntimeError(
                f"refusing to resume with a live rollback journal: {journal}"
            )
        if args.resume_audit is None or not args.resume_audit.is_file():
            raise FileNotFoundError(
                "--resume-audit must name a checked recovery audit when resuming"
            )
        registered_imports = _registered_imports(args.database)
        if registered_imports != args.resume_completed:
            raise RuntimeError(
                "resume checkpoint disagrees with registry: "
                f"--resume-completed={args.resume_completed}, "
                f"registered imports={registered_imports}"
            )
        resume_database_record = file_records([args.database])[0]
    else:
        if args.resume_audit is not None:
            raise ValueError("--resume-audit requires --resume-completed")
        if args.database.exists():
            raise FileExistsError(f"refusing to overwrite registry: {args.database}")
    environment = {
        "PYTHONPATH": str((args.repo / "src").resolve()),
        "MPLCONFIGDIR": str((args.artifact_dir / "matplotlib-cache").resolve()),
    }
    python = sys.executable
    init_command = [
        python,
        "-m",
        "txid",
        "init",
        "--db",
        args.database,
        "--fasta",
        args.reference,
        "--assembly",
        args.assembly,
        "--annotation",
        args.annotation,
        "--annotation-name",
        args.annotation_name,
    ]
    if args.contig_aliases is not None:
        init_command.extend(["--contig-aliases", args.contig_aliases])
    commands = []
    if not args.resume_completed:
        commands.append(
            run_command(
                init_command,
                log_prefix=args.artifact_dir / "00-init",
                env=environment,
                cwd=args.repo,
            )
        )
    observation_rows = []
    snapshot_records = []
    component_records = []
    catalog_builds = []
    for index, (caller, dataset, gtf) in enumerate(rows, 1):
        stem = f"{index:03d}-{dataset}-{caller}"
        rewritten = args.artifact_dir / f"{stem}.txid.gtf"
        mapping = args.artifact_dir / f"{stem}.mapping.tsv"
        if index <= args.resume_completed:
            if args.compress_artifacts:
                rewritten = Path(f"{rewritten}.gz")
                mapping = Path(f"{mapping}.gz")
            missing = [path for path in (rewritten, mapping) if not path.is_file()]
            if missing:
                raise FileNotFoundError(
                    "resume checkpoint is missing component artifacts: "
                    + ", ".join(str(path) for path in missing)
                )
            observation_rows.extend(_read_mapping_rows(mapping))
            component_records.extend(file_records([rewritten, mapping]))
            if args.snapshot_policy == "all":
                snapshot = args.snapshot_dir / f"{index:03d}.catalog.tsv"
                if not snapshot.is_file():
                    raise FileNotFoundError(
                        f"resume checkpoint is missing snapshot: {snapshot}"
                    )
                snapshot_records.append(file_records([snapshot])[0])
            continue
        commands.append(
            run_command(
                [
                    python,
                    "-m",
                    "txid",
                    "add",
                    "--db",
                    args.database,
                    "--input",
                    gtf,
                    "--sample",
                    dataset,
                    "--tool",
                    caller,
                    "--annotation-name",
                    args.annotation_name,
                    "--output-gtf",
                    rewritten,
                    "--mapping",
                    mapping,
                ],
                log_prefix=args.artifact_dir / stem,
                env=environment,
                cwd=args.repo,
            )
        )
        observation_rows.extend(_read_mapping_rows(mapping))
        if args.compress_artifacts:
            rewritten, compression = _compress(
                rewritten,
                log_prefix=args.artifact_dir / f"{stem}.rewritten-gzip",
            )
            commands.append(compression)
            mapping, compression = _compress(
                mapping,
                log_prefix=args.artifact_dir / f"{stem}.mapping-gzip",
            )
            commands.append(compression)
        component_records.extend(file_records([rewritten, mapping]))
        if args.snapshot_policy == "all":
            snapshot = args.snapshot_dir / f"{index:03d}.catalog.tsv"
            catalog_builds.append(
                _write_catalog_from_observations(snapshot, observation_rows)
            )
            snapshot_records.append(file_records([snapshot])[0])
    catalog_builds.append(
        _write_catalog_from_observations(args.catalog, observation_rows)
    )
    commands.append(
        run_command(
            [python, "-m", "txid", "validate", "--db", args.database],
            log_prefix=args.artifact_dir / "final-validate",
            env=environment,
            cwd=args.repo,
        )
    )
    fields = [
        "sample",
        "tool",
        "original_transcript_id",
        "original_gene_id",
        "txid_transcript_id",
        "txid_gene_id",
        "txid_sc",
        "txid_form",
        "classification",
        "annotation_name",
        "gene_candidates",
        "fuzzy_cluster",
        "fuzzy_bridge_status",
    ]
    observation_rows.sort(
        key=lambda row: (row["sample"], row["tool"], row["original_transcript_id"])
    )
    with args.observations.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(observation_rows)
    version = subprocess.run(
        [python, "-m", "txid", "--version"],
        cwd=args.repo,
        env=dict(os.environ, **environment),
        text=True,
        capture_output=True,
        check=False,
    )
    database_artifact = args.database
    if args.compress_artifacts:
        database_artifact, compression = _compress(
            args.database,
            log_prefix=args.artifact_dir / "final-registry-gzip",
        )
        commands.append(compression)
    write_json(
        args.run_json,
        {
            "stage": "txid_registry",
            "version": (version.stdout or version.stderr).strip(),
            "order": args.order,
            "order_seed": seed,
            "snapshot_policy": args.snapshot_policy,
            "resume": (
                None
                if not args.resume_completed
                else {
                    "completed_imports": args.resume_completed,
                    "audit": file_records([args.resume_audit])[0],
                    "database_before_resume": resume_database_record,
                    "timing_scope": (
                        "commands contains only resumed stages; pre-interruption "
                        "runtime is not reconstructed"
                    ),
                }
            ),
            "dataset_tool_order": [
                {"dataset": dataset, "tool": caller, "input": str(gtf.resolve())}
                for caller, dataset, gtf in rows
            ],
            "commands": commands,
            "catalog_builds": catalog_builds,
            "snapshots": snapshot_records,
            "component_artifacts": component_records,
            "artifact_compression": (
                "gzip -n -6" if args.compress_artifacts else None
            ),
            "inputs": file_records(
                [
                    args.reference,
                    args.annotation,
                    *(
                        []
                        if args.contig_aliases is None
                        else [args.contig_aliases]
                    ),
                    *[gtf for _, _, gtf in rows],
                ]
            ),
            "outputs": file_records(
                [database_artifact, args.catalog, args.observations]
            ),
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
