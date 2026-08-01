#!/usr/bin/env python3
"""Run TALON's supported shared-database workflow over identical alignments."""

from __future__ import annotations

import argparse
import csv
import random
import shutil
import subprocess
from pathlib import Path

from benchmark_lib import file_records, run_command, write_json


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-gtf", type=Path, required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--annotation-name", required=True)
    parser.add_argument("--input", action="append", required=True, help="dataset=BAM")
    parser.add_argument("--output-gtf", type=Path, required=True)
    parser.add_argument("--abundance", type=Path, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--run-json", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument(
        "--order",
        default="sorted",
        help="sorted, given, or shuffled-SEED",
    )
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)
    args.output_gtf.parent.mkdir(parents=True, exist_ok=True)
    pairs = []
    for value in args.input:
        dataset, separator, raw_bam = value.partition("=")
        if not separator or not dataset:
            raise ValueError(f"invalid --input value: {value!r}")
        pairs.append((dataset, Path(raw_bam)))
    order_seed = None
    if args.order == "given":
        pass
    else:
        pairs.sort()
        if args.order.startswith("shuffled-"):
            order_seed = int(args.order[len("shuffled-") :])
            random.Random(order_seed).shuffle(pairs)
        elif args.order != "sorted":
            raise ValueError(f"unsupported order: {args.order}")

    commands = []
    config = args.work_dir / "talon-config.csv"
    with config.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        for dataset, bam in pairs:
            writer.writerow([dataset, "TxID benchmark", "long-read", str(bam.resolve())])

    db_prefix = args.work_dir / "talon"
    commands.append(
        run_command(
            [
                "talon_initialize_database",
                f"--f={args.reference_gtf}",
                f"--g={args.assembly}",
                f"--a={args.annotation_name}",
                "--idprefix=TXSMK",
                f"--o={db_prefix}",
            ],
            log_prefix=args.work_dir / "initialize",
            env={"MPLCONFIGDIR": str((args.work_dir / "mpl").resolve())},
        )
    )
    database = db_prefix.with_suffix(".db")
    commands.append(
        run_command(
            [
                "talon",
                "--f",
                str(config),
                "--db",
                str(database),
                "--build",
                args.assembly,
                "--threads",
                str(args.threads),
                "--cov",
                "0.9",
                "--identity",
                "0.8",
                "--nsg",
                "--tmpDir",
                str(args.work_dir / "tmp"),
                "--verbosity",
                "1",
                "--o",
                str(args.work_dir / "talon-run"),
            ],
            log_prefix=args.work_dir / "annotate",
            env={"MPLCONFIGDIR": str((args.work_dir / "mpl").resolve())},
        )
    )
    gtf_prefix = args.work_dir / "observed"
    commands.append(
        run_command(
            [
                "talon_create_GTF",
                "--db",
                str(database),
                "-b",
                args.assembly,
                "-a",
                args.annotation_name,
                "--observed",
                "--o",
                str(gtf_prefix),
            ],
            log_prefix=args.work_dir / "create-gtf",
            env={"MPLCONFIGDIR": str((args.work_dir / "mpl").resolve())},
        )
    )
    abundance_prefix = args.work_dir / "abundance"
    commands.append(
        run_command(
            [
                "talon_abundance",
                "--db",
                str(database),
                "-b",
                args.assembly,
                "-a",
                args.annotation_name,
                "--o",
                str(abundance_prefix),
            ],
            log_prefix=args.work_dir / "abundance",
            env={"MPLCONFIGDIR": str((args.work_dir / "mpl").resolve())},
        )
    )
    # TALON 6.0.1 names observed-only exports
    # ``<prefix>_talon_observedOnly.gtf``.  Older releases used a shorter
    # suffix, so accept both documented families but still require one result.
    observed_gtfs = sorted(args.work_dir.glob("observed*talon*.gtf"))
    abundance_files = sorted(args.work_dir.glob("abundance*talon_abundance.tsv"))
    if len(observed_gtfs) != 1 or len(abundance_files) != 1:
        raise FileNotFoundError(
            f"unexpected TALON outputs: GTF={observed_gtfs}, abundance={abundance_files}"
        )
    shutil.copyfile(observed_gtfs[0], args.output_gtf)
    shutil.copyfile(abundance_files[0], args.abundance)
    shutil.copyfile(database, args.database)
    version_result = subprocess.run(
        ["talon", "--help"], text=True, capture_output=True, check=False
    )
    record = {
        "stage": "talon_shared_database",
        "tool": "TALON",
        "version_contract": "6.0.1",
        "order": args.order,
        "order_seed": order_seed,
        "dataset_order": [dataset for dataset, _ in pairs],
        "input_format": "coordinate-sorted BAM passed directly to TALON 6.0.1",
        "commands": commands,
        "help_signature": (version_result.stdout or version_result.stderr).splitlines()[0],
        "inputs": file_records(
            [args.reference_gtf, *[bam for _, bam in pairs]]
        ),
        "outputs": file_records([args.output_gtf, args.abundance, args.database]),
    }
    write_json(args.run_json, record)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
