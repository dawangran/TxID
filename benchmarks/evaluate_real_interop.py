#!/usr/bin/env python3
"""Audit two real caller outputs across fresh registries using a separate oracle.

Only new outputs are created. The experiment compares exact forms and multi-exon
splice chains; GL1/FC1 and abundance accuracy are outside its scope. The parser
and identity oracle are separately implemented in independent_identity.py.
"""

from __future__ import annotations

import argparse
from collections import Counter
from contextlib import closing
import csv
import json
import os
from pathlib import Path
import sqlite3
import subprocess
import sys
from typing import Callable, Iterable, Mapping


REPOSITORY = Path(__file__).resolve().parents[1]
if str(REPOSITORY) not in sys.path:
    sys.path.insert(0, str(REPOSITORY))
from benchmarks.independent_identity import IndependentTranscript, read_gtf_models, sha256_file


CALLERS = ("isoquant", "stringtie")


def _json(path: Path, value: object) -> None:
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, sort_keys=True, indent=2, ensure_ascii=False)
        handle.write("\n")


def _tsv(path: Path, columns: list[str], rows: Iterable[Mapping[str, object]]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _file(path: Path, role: str) -> dict[str, object]:
    return {"role": role, "path": str(path.resolve()), "bytes": path.stat().st_size,
            "sha256": sha256_file(path)}


class EvaluationError(ValueError):
    pass


class _Commands:
    def __init__(self, output: Path):
        self.directory = output / "commands"
        self.directory.mkdir()
        self.records: list[dict[str, object]] = []

    def run(self, name: str, *arguments: object, raw: bool = False):
        argv = [sys.executable, "-m", "txid", *map(str, arguments)]
        environment = dict(os.environ)
        environment["PYTHONPATH"] = str(REPOSITORY / "src")
        result = subprocess.run(argv, cwd=REPOSITORY, env=environment, text=True,
                                encoding="utf-8", capture_output=True, check=False)
        prefix = self.directory / f"{len(self.records):02d}-{name}"
        stdout, stderr = prefix.with_suffix(".stdout.log"), prefix.with_suffix(".stderr.log")
        stdout.write_text(result.stdout, encoding="utf-8")
        stderr.write_text(result.stderr, encoding="utf-8")
        record = {"argv": argv, "cwd": str(REPOSITORY), "returncode": result.returncode,
                  "stdout": _file(stdout, "stdout"), "stderr": _file(stderr, "stderr")}
        self.records.append(record)
        _json(prefix.with_suffix(".json"), record)
        if result.returncode:
            raise EvaluationError(f"TxID {name} failed; inspect {stderr}")
        return result.stdout.strip() if raw else json.loads(result.stdout)


def _connect(path: Path) -> sqlite3.Connection:
    connection = sqlite3.connect(path.resolve().as_uri() + "?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    return connection


def _reference_context(database: Path) -> tuple[str, dict[str, int], dict[str, str], str]:
    with closing(_connect(database)) as connection:
        context = connection.execute("SELECT fingerprint, fasta_checksum FROM reference_context").fetchone()
        lengths = dict(connection.execute("SELECT name, length FROM contig"))
        aliases = dict(connection.execute("SELECT alias, primary_name FROM contig_alias"))
    return context["fingerprint"], lengths, aliases, context["fasta_checksum"]


def _objects(models: Iterable[IndependentTranscript], fingerprint: str) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for model in models:
        for identity in model.identity(fingerprint).values():
            if identity is None:
                continue
            previous = result.setdefault(identity["public_id"], identity)
            if previous != identity:
                raise EvaluationError("independent oracle public-digest collision")
    return result


def _mapping(path: Path, caller: str, models: Mapping[str, IndependentTranscript],
             fingerprint: str, sample: str) -> dict[tuple[str, str], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    result = {}
    for row in rows:
        identifier = row["original_transcript_id"]
        key = caller, identifier
        if key in result or identifier not in models:
            raise EvaluationError(f"duplicate or unknown observation in {path}: {key}")
        expected = models[identifier].identity(fingerprint)
        expected_sc = "" if expected["splice_chain"] is None else expected["splice_chain"]["public_id"]
        if row["sample"] != sample or row["tool"] != caller:
            raise EvaluationError(f"mapping provenance mismatch for {key}")
        if row["txid_form"] != expected["form"]["public_id"] or row["txid_sc"] != expected_sc:
            raise EvaluationError(f"mapping identity differs from independent oracle for {key}")
        result[key] = row
    if {key[1] for key in result} != set(models):
        raise EvaluationError(f"mapping omits observations in {path}")
    return result


def _audit_database(database: Path, models: Mapping[tuple[str, str], IndependentTranscript],
                    reference_models: Mapping[str, IndependentTranscript], fingerprint: str,
                    sample: str) -> tuple[dict[str, object], dict[tuple[str, str], tuple[str, str]], dict[str, dict[str, str]]]:
    expected_objects = _objects([*reference_models.values(), *models.values()], fingerprint)
    with closing(_connect(database)) as connection:
        stored = {row["public_id"]: {name: row[name] for name in ("public_id", "canonical_json", "full_digest")}
                  for row in connection.execute("SELECT public_id, canonical_json, full_digest FROM structural_object")}
        rows = connection.execute(
            "SELECT m.sample, m.upstream_tool, o.original_transcript_id, o.contig, o.strand, "
            "o.exons_json, o.form_id, o.splice_chain_id FROM observation o "
            "JOIN import_manifest m ON m.id=o.import_id"
        ).fetchall()
    if stored != expected_objects:
        missing = len(set(expected_objects) - set(stored))
        extra = len(set(stored) - set(expected_objects))
        changed = sum(stored[key] != expected_objects[key] for key in set(stored) & set(expected_objects))
        raise EvaluationError(f"stored canonical objects/digests disagree with oracle: missing={missing}, extra={extra}, changed={changed}")
    observations = {}
    for row in rows:
        key = row["upstream_tool"], row["original_transcript_id"]
        if key not in models or key in observations:
            raise EvaluationError(f"unexpected or duplicate database observation {key}")
        model = models[key]
        identity = model.identity(fingerprint)
        expected_sc = "" if identity["splice_chain"] is None else identity["splice_chain"]["public_id"]
        expected_pair = identity["form"]["public_id"], expected_sc
        observed_pair = row["form_id"], row["splice_chain_id"] or ""
        if observed_pair != expected_pair:
            raise EvaluationError(f"database identity differs from oracle for {key}")
        coordinates = tuple(tuple(pair) for pair in json.loads(row["exons_json"]))
        if row["sample"] != sample or row["contig"] != model.contig or row["strand"] != model.strand or coordinates != model.exons:
            raise EvaluationError(f"database coordinates/provenance differ from independent GTF reader for {key}")
        observations[key] = observed_pair
    if set(observations) != set(models):
        raise EvaluationError("database observation set differs from independently parsed inputs")
    return ({"observations_checked": len(observations), "structural_objects_checked": len(stored),
             "canonical_or_full_digest_mismatches": 0, "observation_coordinate_or_identity_mismatches": 0},
            observations, stored)


def _form_target(model: IndependentTranscript) -> tuple:
    return model.contig, model.strand, model.exons


def _chain_target(model: IndependentTranscript) -> tuple:
    # Direct coordinate equality: no digest is used as a truth label.
    introns = tuple((left[1] + 1, right[0] - 1) for left, right in zip(model.exons, model.exons[1:]))
    return model.contig, model.strand, introns if model.strand == "+" else introns[::-1]


def _cross_caller_pairs(left: Mapping[str, IndependentTranscript], right: Mapping[str, IndependentTranscript],
                        predict_left: Callable[[str], str], predict_right: Callable[[str], str],
                        target: Callable[[IndependentTranscript], tuple]) -> dict[str, int | float | None]:
    truth_a = Counter(target(model) for model in left.values())
    truth_b = Counter(target(model) for model in right.values())
    predicted_a = Counter(predict_left(identifier) for identifier in left)
    predicted_b = Counter(predict_right(identifier) for identifier in right)
    joint_a = Counter((predict_left(identifier), target(model)) for identifier, model in left.items())
    joint_b = Counter((predict_right(identifier), target(model)) for identifier, model in right.items())
    true_pairs = sum(count * truth_b[key] for key, count in truth_a.items())
    joined_pairs = sum(count * predicted_b[key] for key, count in predicted_a.items())
    correct_pairs = sum(count * joint_b[key] for key, count in joint_a.items())
    merges, splits = joined_pairs - correct_pairs, true_pairs - correct_pairs
    return {"equal_coordinate_pairs": true_pairs, "joined_pairs": joined_pairs,
            "correctly_joined_pairs": correct_pairs, "false_merge_pairs": merges, "false_split_pairs": splits,
            "false_merge_denominator": joined_pairs, "false_split_denominator": true_pairs,
            "false_merge_rate": merges / joined_pairs if joined_pairs else None,
            "false_split_rate": splits / true_pairs if true_pairs else None}


def _caller_provenance(case: Path, reference: dict[str, object], annotation: dict[str, object],
                       callers: dict[str, dict[str, object]], evidence_kind: str) -> dict[str, object]:
    if evidence_kind == "synthetic-fixture":
        return {"evidence_kind": evidence_kind, "scope": "Synthetic test fixture; not a real caller experiment."}
    path = case / "inputs/selection.json"
    selection = json.loads(path.read_text())
    correction = None
    manifests = [_file(path, "selection_manifest")]
    if selection.get("caller_exit_codes") != {caller: 0 for caller in CALLERS}:
        correction_path = case / "execution-correction.json"
        correction = json.loads(correction_path.read_text())
        if (correction.get("final_caller_exit_codes") != {caller: 0 for caller in CALLERS}
                or correction.get("final_wrapper_exit_code") != 0
                or correction.get("initial_manifest") != "inputs/selection.json"
                or correction.get("retry_record") != "callers/isoquant.run.json"
                or any(correction.get(key) is not False for key in
                       ("changed_selection", "changed_reference_bytes", "changed_model_coordinates", "tuned_on_overlap"))):
            raise EvaluationError("failed initial caller run lacks a valid unchanged-input correction record")
        manifests.append(_file(correction_path, "execution_correction_manifest"))
    if selection["reference_sha256"] != reference["sha256"] or selection["annotation_sha256"] != annotation["sha256"]:
        raise EvaluationError("case reference or annotation checksum differs from caller selection manifest")
    records = {}
    for caller in CALLERS:
        record_path = case / "callers" / f"{caller}.run.json"
        record = json.loads(record_path.read_text())
        if record.get("caller") != caller or record.get("returncode") != 0:
            raise EvaluationError(f"invalid caller execution record for {caller}")
        if callers[caller]["sha256"] not in {entry["sha256"] for entry in record["outputs"]}:
            raise EvaluationError(f"caller output checksum mismatch for {caller}")
        hashes = {entry["sha256"] for entry in record["inputs"]}
        if not {reference["sha256"], annotation["sha256"], selection["chr22_bam_sha256"]}.issubset(hashes):
            raise EvaluationError(f"caller {caller} did not use the locked reference, annotation and common BAM")
        records[caller] = {"version": record["version"], "argv": record["argv"],
                           "compatibility_entry": record.get("compatibility_entry")}
        manifests.append(_file(record_path, f"{caller}_execution_manifest"))
    return {"evidence_kind": evidence_kind, "input_accession": selection["accession"],
            "read_selection": selection["selection"], "selected_reads": selection["reads"],
            "region": selection["region"], "primary_alignment_exclusion_mask": selection["primary_alignment_exclusion_mask"],
            "chr22_primary_alignments": selection["chr22_primary_alignments"],
            "shared_alignment_sha256": selection["chr22_bam_sha256"], "callers": records,
            "initial_caller_exit_codes": selection["caller_exit_codes"],
            "execution_correction": correction, "manifests": manifests}


def evaluate_case(case_dir: Path, reference: Path, *, output: Path | None = None,
                  sample: str = "ENCFF105WIJ_chr22", evidence_kind: str = "real-read-case") -> dict[str, object]:
    case_dir, reference = case_dir.resolve(), reference.resolve()
    output = (case_dir / "interop-evaluation") if output is None else output.resolve()
    if output.exists():
        raise EvaluationError(f"refusing to reuse an evaluation directory: {output}")
    if evidence_kind not in ("real-read-case", "synthetic-fixture"):
        raise EvaluationError("unknown evidence kind")
    annotation = case_dir / "inputs/gencode-v29.identity.chr22.gtf"
    input_paths = {caller: case_dir / "callers" / f"{caller}.gtf" for caller in CALLERS}
    reference_record, annotation_record = _file(reference, "reference_fasta"), _file(annotation, "reference_annotation")
    caller_records = {caller: _file(path, f"{caller}_gtf") for caller, path in input_paths.items()}
    provenance = _caller_provenance(case_dir, reference_record, annotation_record, caller_records, evidence_kind)
    sources = [_file(path, "source_code") for path in [
        Path(__file__), REPOSITORY / "benchmarks/independent_identity.py",
        *sorted((REPOSITORY / "src/txid").glob("*.py")), REPOSITORY / "schemas/registry-v1.sql",
        REPOSITORY / "docs/spec/identity-v1.md"]]
    output.mkdir(parents=True, exist_ok=False)
    commands = _Commands(output)
    _json(output / "input-provenance.json", {"inputs": [reference_record, annotation_record, *caller_records.values()],
                                             "caller_provenance": provenance, "sources": sources})
    try:
        version = commands.run("version", "--version", raw=True)
        models = None
        reference_models = None
        fingerprint = None
        context = None
        initial_mappings = {}
        final_states = []
        audits = []
        for first, second in (CALLERS, CALLERS[::-1]):
            directory = output / f"{first}-first"
            directory.mkdir()
            database = directory / "registry.sqlite"
            init = commands.run(f"{first}-init", "init", "--db", database, "--fasta", reference,
                                "--assembly", "GRCh38-selected-full-reference", "--annotation", annotation,
                                "--annotation-name", "GENCODE-v29-chr22")
            loaded = _reference_context(database)
            if loaded[3].removeprefix("sha256:") != reference_record["sha256"]:
                raise EvaluationError("registered FASTA checksum differs from exact supplied input")
            if context is None:
                context = loaded
                fingerprint = loaded[0]
                models = {caller: read_gtf_models(path, contig_aliases=loaded[2], contig_lengths=loaded[1])
                          for caller, path in input_paths.items()}
                reference_models = read_gtf_models(annotation, contig_aliases=loaded[2], contig_lengths=loaded[1])
            elif loaded != context:
                raise EvaluationError("the two fresh registries have different reference contexts")
            if init["assembly_fingerprint"] != fingerprint:
                raise EvaluationError("CLI and database reference fingerprints disagree")
            accumulated = {}
            previous_observations = None
            previous_objects = None
            for step, caller in enumerate((first, second), 1):
                prefix = directory / f"{step:02d}-{caller}"
                mapping = prefix.with_suffix(".mapping.tsv")
                result = commands.run(f"{first}-{step}-{caller}-add", "add", "--db", database,
                                      "--input", input_paths[caller], "--sample", sample, "--tool", caller,
                                      "--annotation-name", "GENCODE-v29-chr22", "--output-gtf", prefix.with_suffix(".txid.gtf"),
                                      "--mapping", mapping)
                if not result["created"]:
                    raise EvaluationError("a fresh experiment unexpectedly reused a previous import")
                current_mapping = _mapping(mapping, caller, models[caller], fingerprint, sample)
                if step == 1:
                    initial_mappings[caller] = current_mapping
                accumulated.update({(caller, identifier): model for identifier, model in models[caller].items()})
                audit, observations, objects = _audit_database(database, accumulated, reference_models, fingerprint, sample)
                validation = commands.run(f"{first}-{step}-validate", "validate", "--db", database)
                if not validation["valid"]:
                    raise EvaluationError("TxID validate returned an invalid registry")
                commands.run(f"{first}-{step}-export", "export", "--db", database,
                             "--catalog", prefix.with_suffix(".catalog.tsv"))
                changes = disappearances = object_changes = 0
                if previous_observations is not None:
                    disappearances = len(set(previous_observations) - set(observations))
                    changes = sum(observations.get(key) != pair for key, pair in previous_observations.items())
                    object_changes = sum(objects.get(key) != value for key, value in previous_objects.items())
                    if disappearances or changes or object_changes:
                        raise EvaluationError("incremental import changed an existing exact observation or canonical object")
                audit.update({"registry": first + "-first", "step": step, "imported_caller": caller,
                              "mapping_observations_checked": len(current_mapping),
                              "existing_exact_observation_changes": changes,
                              "existing_exact_observation_disappearances": disappearances,
                              "existing_canonical_or_digest_changes": object_changes})
                audits.append(audit)
                previous_observations, previous_objects = observations, objects
            final_states.append((observations, objects))
        if final_states[0] != final_states[1]:
            raise EvaluationError("reverse import order changed exact observations or stored canonical objects")

        # Join only the first-import mappings, before either registry saw the other caller.
        independently_assigned = {key: row for caller in CALLERS for key, row in initial_mappings[caller].items()}
        all_models = {(caller, identifier): model for caller in CALLERS for identifier, model in models[caller].items()}
        evidence = []
        for (caller, identifier), model in sorted(all_models.items()):
            identity = model.identity(fingerprint)
            form, chain = identity["form"], identity["splice_chain"]
            evidence.append({"caller": caller, "original_transcript_id": identifier, "original_gene_id": model.gene_id or "",
                             "contig": model.contig, "strand": model.strand, "exons": json.dumps(model.exons, separators=(",", ":")),
                             "form_id": form["public_id"], "splice_chain_id": chain["public_id"] if chain else "",
                             "form_full_digest": form["full_digest"], "form_canonical_json": form["canonical_json"],
                             "splice_full_digest": chain["full_digest"] if chain else "",
                             "splice_canonical_json": chain["canonical_json"] if chain else "",
                             "initial_mapping_matches_oracle": True, "both_final_registries_match_oracle": True})
        _tsv(output / "per-observation-evidence.tsv", list(evidence[0]), evidence)
        sets = {caller: {field: {row[field] for row in independently_assigned.values()
                               if row["tool"] == caller and row[field]} for field in ("txid_form", "txid_sc")}
                for caller in CALLERS}
        matrix_metrics = {}
        for field, name in (("txid_form", "form"), ("txid_sc", "chain")):
            columns = sorted(set.union(*(sets[caller][field] for caller in CALLERS)))
            _tsv(output / f"caller-by-{name}.tsv", ["caller", *columns],
                 ({"caller": caller, **{key: int(key in sets[caller][field]) for key in columns}} for caller in CALLERS))
            occupied = sum(len(sets[caller][field]) for caller in CALLERS)
            matrix_metrics[name] = {"rows": 2, "columns": len(columns), "occupied_cells": occupied,
                                    "shared_columns": len(set.intersection(*(sets[caller][field] for caller in CALLERS))),
                                    "sparsity": 1 - occupied / (2 * len(columns)) if columns else None}
        left, right = (models[caller] for caller in CALLERS)
        raw_pairs = _cross_caller_pairs(left, right, lambda name: name, lambda name: name, _form_target)
        exact_pairs = _cross_caller_pairs(left, right,
            lambda name: independently_assigned[(CALLERS[0], name)]["txid_form"],
            lambda name: independently_assigned[(CALLERS[1], name)]["txid_form"], _form_target)
        multi = {caller: {name: model for name, model in models[caller].items() if len(model.exons) > 1} for caller in CALLERS}
        chain_pairs = _cross_caller_pairs(multi[CALLERS[0]], multi[CALLERS[1]],
            lambda name: independently_assigned[(CALLERS[0], name)]["txid_sc"],
            lambda name: independently_assigned[(CALLERS[1], name)]["txid_sc"], _chain_target)
        shared_chains = sets[CALLERS[0]]["txid_sc"] & sets[CALLERS[1]]["txid_sc"]
        end_variant_chains = []
        for chain in sorted(shared_chains):
            forms = {row["txid_form"] for row in independently_assigned.values() if row["txid_sc"] == chain}
            if len(forms) > 1:
                end_variant_chains.append(chain)
        same_chain_form_pairs = _cross_caller_pairs(multi[CALLERS[0]], multi[CALLERS[1]],
            lambda name: independently_assigned[(CALLERS[0], name)]["txid_sc"],
            lambda name: independently_assigned[(CALLERS[1], name)]["txid_sc"], _form_target)
        for result in (exact_pairs, chain_pairs):
            if result["false_merge_pairs"] or result["false_split_pairs"]:
                raise EvaluationError("independently assigned exact keys do not reproduce coordinate equality")
        summary = {"schema": "txid.real-interop-evaluation.v1", "valid": True, "evidence_kind": evidence_kind,
                   "software_version": version, "python_version": sys.version.split()[0], "sample": sample,
                   "assembly_fingerprint": fingerprint, "caller_provenance": provenance,
                   "scope": "Two caller outputs from one biological input in the real-read case; separately initialized registries. Exact keys only; no GL1/FC1, discovery, abundance or biological-truth comparison. Reference fingerprint supplied by registered common FASTA; independent GTF parsing and coordinate/identity calculation.",
                   "independent_initial_registries": 2, "import_orders": [list(CALLERS), list(CALLERS[::-1])],
                   "caller_counts": {caller: {"observations": len(models[caller]), "multi_exon_observations": len(multi[caller]),
                                               "single_exon_observations": len(models[caller]) - len(multi[caller]),
                                               "distinct_forms": len(sets[caller]["txid_form"]),
                                               "distinct_multi_exon_chains": len(sets[caller]["txid_sc"])} for caller in CALLERS},
                   "original_id_intersection": len(set(left) & set(right)),
                   "coordinate_join_of_original_ids": raw_pairs, "exact_form_join": exact_pairs,
                   "multi_exon_chain_join": chain_pairs, "matrices": matrix_metrics,
                   "shared_chains_with_distinct_forms": len(end_variant_chains),
                   "same_chain_different_form_cross_caller_pairs": same_chain_form_pairs["false_merge_pairs"],
                   "end_variant_chain_ids": end_variant_chains,
                   "stage_audits": audits, "reverse_order_exact_observation_changes": 0,
                   "reverse_order_canonical_or_digest_changes": 0,
                   "inputs": [reference_record, annotation_record, *caller_records.values()], "sources": sources}
        metrics = []
        for section in ("exact_form_join", "multi_exon_chain_join", "coordinate_join_of_original_ids"):
            for key, value in summary[section].items():
                metrics.append({"metric": section + "." + key, "value": "NA" if value is None else value,
                                "unit": "fraction" if key.endswith("rate") else "cross-caller observation pairs"})
        for name, values in matrix_metrics.items():
            for key, value in values.items():
                unit = {"rows": "caller outputs", "columns": "identity groups", "shared_columns": "identity groups", "occupied_cells": "caller-group cells", "sparsity": "fraction"}[key]
                metrics.append({"metric": f"{name}_matrix.{key}", "value": "NA" if value is None else value, "unit": unit})
        for key, unit in (("original_id_intersection", "identifier strings"), ("shared_chains_with_distinct_forms", "splice chains"),
                          ("same_chain_different_form_cross_caller_pairs", "cross-caller observation pairs")):
            metrics.append({"metric": key, "value": summary[key], "unit": unit})
        _tsv(output / "metric.tsv", ["metric", "value", "unit"], metrics)
        _json(output / "commands.json", commands.records)
        summary["artifacts"] = [_file(output / name, "evidence_artifact") for name in
                                ("per-observation-evidence.tsv", "caller-by-form.tsv", "caller-by-chain.tsv", "metric.tsv", "commands.json", "input-provenance.json")]
        _json(output / "summary.json", summary)
        return summary
    except Exception as error:
        _json(output / "failure.json", {"valid": False, "error_type": type(error).__name__,
                                       "error": str(error), "completed_commands": len(commands.records),
                                       "note": "Partial artifacts retained for diagnosis; this is not a successful evaluation."})
        raise


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="new directory; default CASE/interop-evaluation")
    parser.add_argument("--sample", default="ENCFF105WIJ_chr22")
    parser.add_argument("--evidence-kind", choices=("real-read-case", "synthetic-fixture"), default="real-read-case")
    args = parser.parse_args(argv)
    try:
        result = evaluate_case(args.case_dir, args.reference, output=args.output,
                               sample=args.sample, evidence_kind=args.evidence_kind)
        print(json.dumps({"valid": result["valid"], "evidence_kind": result["evidence_kind"],
                          "caller_counts": result["caller_counts"], "matrices": result["matrices"]}, sort_keys=True))
        return 0
    except (OSError, ValueError, KeyError, TypeError, sqlite3.Error) as error:
        print(f"real-interop: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
