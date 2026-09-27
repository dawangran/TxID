#!/usr/bin/env python3
"""Run the synthetic example using the installed TxID package."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="new output directory; defaults to a temporary directory")
    args = parser.parse_args()
    inputs = Path(__file__).resolve().parent
    if args.output_dir is None:
        output = Path(tempfile.mkdtemp(prefix="txid-example-"))
    else:
        output = args.output_dir.resolve()
        if output.exists():
            parser.error(f"output directory already exists: {output}")
        output.mkdir(parents=True)
    database = output / "cohort.sqlite"

    def run(*arguments: object) -> dict:
        command = [sys.executable, "-m", "txid", *map(str, arguments)]
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, text=True)
        return json.loads(result.stdout)

    run("init", "--db", database, "--fasta", inputs / "reference.fa",
        "--assembly", "synthetic-v1", "--annotation", inputs / "reference.gtf",
        "--annotation-name", "ref-v1")
    import_args = ("add", "--db", database, "--input", inputs / "caller.gtf",
                   "--sample", "demo", "--tool", "synthetic-fixture",
                   "--annotation-name", "ref-v1", "--output-gtf", output / "demo.txid.gtf",
                   "--mapping", output / "demo.mapping.tsv")
    imported = run(*import_args)
    catalog = run("export", "--db", database, "--catalog", output / "catalog.tsv")
    validation = run("validate", "--db", database)
    run("plot", "--db", database, "--gene", "g1", "--output", output / "g1.svg")
    retry = run(*import_args)

    with (output / "demo.mapping.tsv").open(encoding="utf-8", newline="") as handle:
        rows = {row["original_transcript_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
    known, end, splice = (rows[name] for name in ("alt_known", "end_variant", "one_bp_splice"))
    observed = {
        "transcripts": imported["transcripts"],
        "catalog_forms": catalog["rows"],
        "classifications": validation["summary"]["classifications"],
        "valid": validation["valid"],
        "reference_id_preserved": known["txid_transcript_id"] == "tx_ref1",
        "end_variant_same_splice_chain": known["txid_sc"] == end["txid_sc"],
        "end_variant_distinct_form": known["txid_form"] != end["txid_form"],
        "splice_change_distinct_chain": known["txid_sc"] != splice["txid_sc"],
        "retry_created": retry["created"],
    }
    expected = json.loads((inputs / "expected.json").read_text(encoding="utf-8"))
    (output / "summary.json").write_text(json.dumps(observed, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if observed != expected:
        print(json.dumps({"expected": expected, "observed": observed}, indent=2), file=sys.stderr)
        return 1
    print(json.dumps({"output_dir": str(output), "checks_passed": True, **observed}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
