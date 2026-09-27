"""Small, explicitly synthetic tests of the real-caller interoperability audit."""

import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from benchmarks.evaluate_real_interop import (
    EvaluationError, REPOSITORY, _caller_provenance, _cross_caller_pairs, _form_target, evaluate_case,
)
from benchmarks.independent_identity import IndependentTranscript, sha256_file


def _gtf(models):
    lines = []
    for identifier, gene, strand, exons in models:
        attributes = f'gene_id "{gene}"; transcript_id "{identifier}";'
        lines.append(f"chr22\tfixture\ttranscript\t{exons[0][0]}\t{exons[-1][1]}\t.\t{strand}\t.\t{attributes}\n")
        for start, end in exons:
            lines.append(f"chr22\tfixture\texon\t{start}\t{end}\t.\t{strand}\t.\t{attributes}\n")
    return "".join(lines)


def _fixture(case):
    (case / "callers").mkdir(parents=True)
    (case / "inputs").mkdir()
    fasta = case / "reference.fa"
    fasta.write_text(">chr22\n" + "A" * 1200 + "\n")
    (case / "inputs/gencode-v29.identity.chr22.gtf").write_text(_gtf([
        ("R1", "GR", "+", [(100, 149), (200, 249)]),
    ]))
    (case / "callers/isoquant.gtf").write_text(_gtf([
        ("q_shared", "q1", "+", [(100, 149), (200, 249)]),
        ("q_end", "q2", "-", [(300, 349), (400, 449)]),
        ("q_unique", "q3", "+", [(700, 740)]),
        ("collision", "q4", "+", [(800, 849), (900, 949)]),
    ]))
    (case / "callers/stringtie.gtf").write_text(_gtf([
        ("s_shared", "s1", "+", [(100, 149), (200, 249)]),
        ("s_end", "s2", "-", [(290, 349), (400, 460)]),
        ("s_unique", "s3", "+", [(1000, 1040)]),
        ("collision", "s4", "+", [(50, 60)]),
    ]))
    return fasta


class RealInteropTests(unittest.TestCase):
    def test_two_fresh_registries_increment_reverse_order_and_matrix(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            fasta = _fixture(case)
            result = subprocess.run([
                sys.executable, str(REPOSITORY / "benchmarks/evaluate_real_interop.py"),
                "--case-dir", str(case), "--reference", str(fasta),
                "--evidence-kind", "synthetic-fixture", "--sample", "fixture_one_sample",
            ], capture_output=True, text=True, check=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(json.loads(result.stdout)["valid"])
            output = case / "interop-evaluation"
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual(summary["evidence_kind"], "synthetic-fixture")
            self.assertEqual(summary["matrices"]["form"]["columns"], 7)
            self.assertEqual(summary["matrices"]["form"]["shared_columns"], 1)
            self.assertEqual(summary["matrices"]["form"]["occupied_cells"], 8)
            self.assertEqual(summary["matrices"]["chain"]["columns"], 3)
            self.assertEqual(summary["matrices"]["chain"]["shared_columns"], 2)
            self.assertEqual(summary["shared_chains_with_distinct_forms"], 1)
            self.assertEqual(summary["same_chain_different_form_cross_caller_pairs"], 1)
            self.assertEqual(summary["original_id_intersection"], 1)
            raw = summary["coordinate_join_of_original_ids"]
            self.assertEqual(raw["false_merge_pairs"], 1)
            self.assertEqual(raw["false_split_pairs"], 1)
            self.assertEqual(raw["false_merge_rate"], 1)
            self.assertEqual(raw["false_split_rate"], 1)
            for key, denominator in (("exact_form_join", 1), ("multi_exon_chain_join", 2)):
                self.assertEqual(summary[key]["false_merge_pairs"], 0)
                self.assertEqual(summary[key]["false_split_pairs"], 0)
                self.assertEqual(summary[key]["equal_coordinate_pairs"], denominator)
            self.assertEqual(len(summary["stage_audits"]), 4)
            for audit in summary["stage_audits"]:
                self.assertEqual(audit["observations_checked"], audit["step"] * 4)
                for key in ("canonical_or_full_digest_mismatches", "observation_coordinate_or_identity_mismatches",
                            "existing_exact_observation_changes", "existing_canonical_or_digest_changes"):
                    self.assertEqual(audit[key], 0)
            self.assertEqual(summary["reverse_order_exact_observation_changes"], 0)
            self.assertEqual(summary["reverse_order_canonical_or_digest_changes"], 0)
            with (output / "per-observation-evidence.tsv").open() as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 8)
            self.assertEqual({row["form_id"].split(".")[0] for row in rows}, {"txid:TF1", "txid:SE1"})
            for row in rows:
                self.assertEqual(row["initial_mapping_matches_oracle"], "True")
                self.assertEqual(row["both_final_registries_match_oracle"], "True")
                self.assertEqual(len(row["form_full_digest"]), 64)
            for artifact in summary["artifacts"]:
                self.assertEqual(sha256_file(artifact["path"]), artifact["sha256"])
            for layer in ("form", "chain"):
                with (output / f"caller-by-{layer}.tsv").open() as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertEqual(len(rows), 2)
                self.assertTrue(all(value in ("0", "1") for row in rows for key, value in row.items() if key != "caller"))
            previous_hash = sha256_file(output / "summary.json")
            with self.assertRaisesRegex(EvaluationError, "refusing to reuse"):
                evaluate_case(case, fasta, evidence_kind="synthetic-fixture")
            self.assertEqual(sha256_file(output / "summary.json"), previous_hash)

    def test_real_evidence_requires_caller_manifests(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            fasta = _fixture(case)
            with self.assertRaises(FileNotFoundError):
                evaluate_case(case, fasta)
            self.assertFalse((case / "interop-evaluation").exists())

    def test_pair_counts_are_cross_caller_observation_pairs(self):
        def model(identifier, exons):
            return IndependentTranscript(identifier, None, "chr22", "+", tuple(exons))
        left = {"a": model("a", [(1, 10)]), "b": model("b", [(1, 10)])}
        right = {"c": model("c", [(1, 10)]), "d": model("d", [(1, 10)]), "e": model("e", [(20, 30)])}
        result = _cross_caller_pairs(left, right, lambda _: "joined", lambda _: "joined", _form_target)
        self.assertEqual(result["equal_coordinate_pairs"], 4)
        self.assertEqual(result["joined_pairs"], 6)
        self.assertEqual(result["correctly_joined_pairs"], 4)
        self.assertEqual(result["false_merge_pairs"], 2)
        self.assertEqual(result["false_merge_denominator"], 6)
        self.assertEqual(result["false_split_denominator"], 4)
        self.assertAlmostEqual(result["false_merge_rate"], 1 / 3)
        empty = _cross_caller_pairs({}, {}, str, str, _form_target)
        self.assertIsNone(empty["false_merge_rate"])
        self.assertIsNone(empty["false_split_rate"])

    def test_corrected_execution_keeps_failed_history_and_locks_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            case = Path(directory)
            _fixture(case)
            selection = {
                "caller_exit_codes": {"isoquant": 1, "stringtie": 0},
                "reference_sha256": "reference", "annotation_sha256": "annotation",
                "chr22_bam_sha256": "alignment", "accession": "synthetic",
                "selection": "synthetic fixture", "reads": 1, "region": "chr22",
                "primary_alignment_exclusion_mask": 2308, "chr22_primary_alignments": 1,
            }
            (case / "inputs/selection.json").write_text(json.dumps(selection))
            correction = {
                "final_caller_exit_codes": {"isoquant": 0, "stringtie": 0},
                "final_wrapper_exit_code": 0, "initial_manifest": "inputs/selection.json",
                "retry_record": "callers/isoquant.run.json", "changed_selection": False,
                "changed_reference_bytes": False, "changed_model_coordinates": False,
                "tuned_on_overlap": False,
            }
            correction_path = case / "execution-correction.json"
            correction_path.write_text(json.dumps(correction))
            callers = {name: {"sha256": name} for name in ("isoquant", "stringtie")}
            for name in callers:
                record = {"caller": name, "returncode": 0, "version": "fixture", "argv": [],
                          "inputs": [{"sha256": key} for key in ("reference", "annotation", "alignment")],
                          "outputs": [callers[name]]}
                (case / "callers" / f"{name}.run.json").write_text(json.dumps(record))
            result = _caller_provenance(case, {"sha256": "reference"}, {"sha256": "annotation"}, callers, "real-read-case")
            self.assertEqual(result["initial_caller_exit_codes"]["isoquant"], 1)
            self.assertEqual(result["execution_correction"]["final_caller_exit_codes"]["isoquant"], 0)
            self.assertEqual(len(result["manifests"]), 4)
            with self.assertRaisesRegex(EvaluationError, "checksum differs"):
                _caller_provenance(case, {"sha256": "wrong"}, {"sha256": "annotation"}, callers, "real-read-case")
            correction["changed_reference_bytes"] = True
            correction_path.write_text(json.dumps(correction))
            with self.assertRaisesRegex(EvaluationError, "correction record"):
                _caller_provenance(case, {"sha256": "reference"}, {"sha256": "annotation"}, callers, "real-read-case")


if __name__ == "__main__":
    unittest.main()
