from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from contextlib import closing
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflows" / "publication_benchmark" / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def load(name: str):
    path = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BENCHMARK_LIB = load("benchmark_lib")
GENERATOR = load("generate_smoke_data")
VALIDATOR = load("validate_manifests")
IDENTITY_EVALUATOR = load("evaluate_identity_partitions")
TALON_EVALUATOR = load("evaluate_talon_partitions")
ENCODE_STAGER = load("stage_encode_files")
SIRV_TRUTH = load("build_sirv_truth")
REFERENCE_STAGER = load("stage_public_references")
ENCODE_TALON_EXTRACTOR = load("extract_encode_talon_models")
GTF_IDENTITY_VIEW = load("make_gtf_identity_view")
NCBI_CONTIG_ALIASES = load("build_ncbi_contig_aliases")
TXID_REGISTRY_RUNNER = load("run_txid_registry")
TXID_INVARIANCE_EVALUATOR = load("evaluate_txid_invariance")
FIXED_MODEL_CLASSIFIER = load("classify_fixed_models")


class PublicationBenchmarkTests(unittest.TestCase):
    def test_txid_registry_resume_helpers_count_imports_and_read_gzip_mapping(self):
        import gzip
        import sqlite3

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = root / "registry.sqlite"
            with closing(sqlite3.connect(database)) as connection:
                connection.execute(
                    "CREATE TABLE import_manifest (id INTEGER PRIMARY KEY)"
                )
                connection.executemany(
                    "INSERT INTO import_manifest(id) VALUES (?)",
                    [(1,), (2,), (3,)],
                )
                connection.commit()
            self.assertEqual(TXID_REGISTRY_RUNNER._registered_imports(database), 3)

            payload = "sample\ttool\toriginal_transcript_id\ns\tt\ttx1\n"
            plain = root / "mapping.tsv"
            plain.write_text(payload, encoding="utf-8")
            compressed = root / "mapping.tsv.gz"
            with gzip.open(compressed, "wt", encoding="utf-8", newline="") as handle:
                handle.write(payload)
            expected = [
                {
                    "sample": "s",
                    "tool": "t",
                    "original_transcript_id": "tx1",
                }
            ]
            self.assertEqual(
                TXID_REGISTRY_RUNNER._read_mapping_rows(plain), expected
            )
            self.assertEqual(
                TXID_REGISTRY_RUNNER._read_mapping_rows(compressed), expected
            )
            self.assertEqual(TXID_INVARIANCE_EVALUATOR._read(plain), expected)
            self.assertEqual(TXID_INVARIANCE_EVALUATOR._read(compressed), expected)
            self.assertEqual(
                TXID_INVARIANCE_EVALUATOR._payload_bytes(plain),
                TXID_INVARIANCE_EVALUATOR._payload_bytes(compressed),
            )
            self.assertEqual(
                TXID_INVARIANCE_EVALUATOR._sha256(plain),
                TXID_INVARIANCE_EVALUATOR._sha256(compressed),
            )

            catalog = root / "catalog.tsv"
            rows = [
                {
                    "sample": "s2",
                    "tool": "t2",
                    "txid_form": "txid:TF1.b",
                    "txid_sc": "txid:SC1.b",
                    "txid_gene_id": "GENE2",
                    "classification": "novel_in_known_gene",
                },
                {
                    "sample": "s1",
                    "tool": "t1",
                    "txid_form": "txid:TF1.a",
                    "txid_sc": "txid:SC1.a",
                    "txid_gene_id": "GENE1",
                    "classification": "known",
                },
                {
                    "sample": "s2",
                    "tool": "t2",
                    "txid_form": "txid:TF1.a",
                    "txid_sc": "txid:SC1.a",
                    "txid_gene_id": "GENE1",
                    "classification": "novel_in_known_gene",
                },
            ]
            record = TXID_REGISTRY_RUNNER._write_catalog_from_observations(
                catalog, rows
            )
            self.assertEqual(record["catalog_rows"], 2)
            self.assertEqual(
                catalog.read_text(encoding="utf-8"),
                "form_id\tsplice_chain_id\tgene_id\tclassifications\t"
                "observation_count\tsample_count\ttool_count\n"
                "txid:TF1.a\ttxid:SC1.a\tGENE1\t"
                "known,novel_in_known_gene\t2\t2\t2\n"
                "txid:TF1.b\ttxid:SC1.b\tGENE2\t"
                "novel_in_known_gene\t1\t1\t1\n",
            )

            snapshots = root / "snapshots"
            snapshots.mkdir()
            first_snapshot = snapshots / "001.catalog.tsv"
            second_snapshot = snapshots / "002.catalog.tsv"
            header = (
                "form_id\tsplice_chain_id\tgene_id\tclassifications\t"
                "observation_count\tsample_count\ttool_count\n"
            )
            first_snapshot.write_text(
                header
                + "txid:TF1.a\ttxid:SC1.a\tGENE1\tknown\t1\t1\t1\n",
                encoding="utf-8",
            )
            second_snapshot.write_text(
                header
                + "txid:TF1.a\ttxid:SC1.a\tGENE1\tknown\t1\t1\t1\n"
                + "txid:TF1.a\ttxid:SC1.a\tGENE2\t"
                "novel_in_known_gene\t1\t1\t1\n",
                encoding="utf-8",
            )
            incremental = TXID_INVARIANCE_EVALUATOR._incremental_metrics(
                snapshots, second_snapshot
            )
            self.assertEqual(
                incremental["incremental_identifier_disappearances"], 0
            )
            self.assertEqual(
                incremental["incremental_structural_metadata_changes"], 0
            )
            self.assertEqual(
                incremental["incremental_classification_regressions"], 0
            )
            self.assertTrue(
                incremental["incremental_final_catalog_bytes_identical"]
            )

            from tests.helpers import make_registry

            classify_root = root / "classification-only"
            classify_root.mkdir()
            classify_database, classify_input = make_registry(classify_root)
            classify_output = classify_root / "classification.tsv.gz"
            classification_record = FIXED_MODEL_CLASSIFIER.classify_models(
                database=classify_database,
                annotation_name="ref-v1",
                input_path=classify_input,
                sample="sample",
                tool="caller",
                output=classify_output,
            )
            self.assertEqual(classification_record["transcripts"], 9)
            with closing(sqlite3.connect(classify_database)) as connection:
                self.assertEqual(
                    connection.execute(
                        "SELECT COUNT(*) FROM import_manifest"
                    ).fetchone()[0],
                    0,
                )
            with gzip.open(
                classify_output, "rt", encoding="utf-8", newline=""
            ) as handle:
                classified_rows = list(
                    csv.DictReader(handle, delimiter="\t")
                )
            self.assertEqual(len(classified_rows), 9)
            self.assertEqual(
                {row["sample"] for row in classified_rows}, {"sample"}
            )

    def test_encode_gtf_reference_contexts_pin_primary_assembly(self):
        path = (
            ROOT
            / "workflows"
            / "publication_benchmark"
            / "config"
            / "references.encode-gtf.public.tsv"
        )
        with path.open(newline="", encoding="utf-8") as handle:
            rows = {
                row["annotation_id"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(set(rows), {"gencode-v29", "gencode-v49"})
        self.assertEqual(
            rows["gencode-v49"]["annotation_md5"],
            "8486a6bdcd27a8a7a08232d01cc13b77",
        )
        self.assertTrue(
            rows["gencode-v49"]["annotation"].endswith(
                "gencode.v49.primary_assembly.annotation.gtf.gz"
            )
        )
        self.assertNotIn(".v49.annotation.gtf.gz", rows["gencode-v49"]["annotation"])

    def test_ncbi_assembly_report_aliases_require_matching_fasta_lengths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fai = root / "reference.fa.fai"
            fai.write_text("chr1\t10\t0\t10\t11\nchrEBV\t5\t11\t5\t6\n", encoding="utf-8")
            report = root / "assembly_report.txt"
            report.write_text(
                "# Assembly name: GRCh38\n"
                "1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\t"
                "NC_000001.11\tPrimary Assembly\t10\tchr1\n",
                encoding="utf-8",
            )
            output = root / "aliases.tsv"
            record = NCBI_CONTIG_ALIASES.build(report, fai, output)
            self.assertEqual(record["report_targets_present_in_fasta"], 1)
            self.assertEqual(record["fasta_contigs_without_report_target"], ["chrEBV"])
            self.assertEqual(
                output.read_text(encoding="utf-8"),
                "alias\tprimary\n"
                "1\tchr1\n"
                "CM000663.2\tchr1\n"
                "NC_000001.11\tchr1\n",
            )
            report.write_text(
                "1\tassembled-molecule\t1\tChromosome\tCM000663.2\t=\t"
                "NC_000001.11\tPrimary Assembly\t11\tchr1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "conflicts with FASTA length"):
                NCBI_CONTIG_ALIASES.build(report, fai, root / "invalid.tsv")

    def test_clustered_ratio_bootstrap_is_deterministic_and_contains_point(self):
        strata = [(0, 10), (2, 10), (5, 10), (0, 0)]
        first = BENCHMARK_LIB.clustered_ratio_bootstrap(
            strata, replicates=500, seed=20260731
        )
        second = BENCHMARK_LIB.clustered_ratio_bootstrap(
            strata, replicates=500, seed=20260731
        )
        self.assertEqual(first, second)
        point = 7 / 30
        self.assertLessEqual(first["lower"], point)
        self.assertGreaterEqual(first["upper"], point)

    def test_transcript_bootstrap_is_deterministic_and_contains_point(self):
        first = BENCHMARK_LIB.binary_precision_recall_bootstrap(
            [True, True, False, True],
            [True, False, True],
            replicates=500,
            seed=20260731,
        )
        second = BENCHMARK_LIB.binary_precision_recall_bootstrap(
            [True, True, False, True],
            [True, False, True],
            replicates=500,
            seed=20260731,
        )
        self.assertEqual(first, second)
        self.assertLessEqual(first["recall_ci95_lower"], 0.75)
        self.assertGreaterEqual(first["recall_ci95_upper"], 0.75)
        self.assertLessEqual(first["precision_ci95_lower"], 2 / 3)
        self.assertGreaterEqual(first["precision_ci95_upper"], 2 / 3)

    def test_git_source_fetch_can_use_a_checked_local_cache(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            subprocess.run(["git", "init", source], check=True, capture_output=True)
            (source / "tool.py").write_text("print('checked')\n", encoding="utf-8")
            subprocess.run(["git", "-C", source, "add", "tool.py"], check=True)
            subprocess.run(
                [
                    "git",
                    "-C",
                    source,
                    "-c",
                    "user.name=TxID test",
                    "-c",
                    "user.email=txid-test@localhost",
                    "commit",
                    "-m",
                    "fixture",
                ],
                check=True,
                capture_output=True,
            )
            commit = subprocess.run(
                ["git", "-C", source, "rev-parse", "HEAD"],
                check=True,
                capture_output=True,
                text=True,
            ).stdout.strip()
            target = root / "checkout"
            record = root / "fetch.json"
            subprocess.run(
                [
                    sys.executable,
                    SCRIPTS / "fetch_git_source.py",
                    "--repository",
                    "https://example.invalid/upstream.git",
                    "--cache",
                    source,
                    "--revision",
                    commit,
                    "--target",
                    target,
                    "--record",
                    record,
                ],
                check=True,
            )
            audit = json.loads(record.read_text(encoding="utf-8"))
            self.assertEqual(audit["repository"], "https://example.invalid/upstream.git")
            self.assertTrue(audit["used_local_cache"])
            self.assertEqual(audit["resolved_commit"], commit)
            self.assertEqual((target / "tool.py").read_text(), "print('checked')\n")

    def test_run_command_preserves_non_utf8_tool_diagnostics(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record = BENCHMARK_LIB.run_command(
                [
                    sys.executable,
                    "-c",
                    "import sys; sys.stderr.buffer.write(bytes([255, 10]))",
                ],
                log_prefix=root / "binary-diagnostic",
            )
            self.assertEqual(record["returncode"], 0)
            diagnostic = (root / "binary-diagnostic.stderr.log").read_text(
                encoding="utf-8"
            )
            self.assertEqual(diagnostic, "\ufffd\n")

    def test_git_source_failure_never_occupies_final_checkout(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            target = root / "checkout"
            result = subprocess.run(
                [
                    sys.executable,
                    SCRIPTS / "fetch_git_source.py",
                    "--repository",
                    str(root / "absent-repository"),
                    "--revision",
                    "main",
                    "--target",
                    target,
                    "--record",
                    root / "fetch.json",
                ],
                capture_output=True,
                text=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertFalse(target.exists())
            source = root / "absent-repository"
            subprocess.run(["git", "init", source], check=True, capture_output=True)
            (source / "tool.py").write_text("print('retry')\n", encoding="utf-8")
            subprocess.run(["git", "-C", source, "add", "tool.py"], check=True)
            subprocess.run(
                [
                    "git",
                    "-C",
                    source,
                    "-c",
                    "user.name=TxID test",
                    "-c",
                    "user.email=txid-test@localhost",
                    "commit",
                    "-m",
                    "retry fixture",
                ],
                check=True,
                capture_output=True,
            )
            subprocess.run(
                [
                    sys.executable,
                    SCRIPTS / "fetch_git_source.py",
                    "--repository",
                    str(source),
                    "--revision",
                    "HEAD",
                    "--target",
                    target,
                    "--record",
                    root / "retry.json",
                ],
                check=True,
            )
            self.assertEqual((target / "tool.py").read_text(), "print('retry')\n")

    def test_encode_stager_checks_metadata_and_never_overwrites_source_manifest(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = b"checked ENCODE fixture\n"
            accession = "ENCFFTEST001"
            target = root / "inputs" / f"{accession}.fastq.gz"
            target.parent.mkdir()
            target.write_bytes(data)
            metadata_dir = root / "metadata"
            metadata_dir.mkdir()
            metadata = {
                "accession": accession,
                "status": "released",
                "file_format": "fastq",
                "file_size": len(data),
                "md5sum": __import__("hashlib").md5(data).hexdigest(),
                "href": f"/files/{accession}/@@download/{accession}.fastq.gz",
                "dataset": "/experiments/ENCSRTEST/",
                "biological_replicates": [1],
                "technical_replicates": ["1_1"],
            }
            (metadata_dir / f"{accession}.json").write_text(
                json.dumps(metadata), encoding="utf-8"
            )
            manifest = root / "datasets.tsv"
            fields = [
                "dataset_id",
                "public_accession",
                "input_kind",
                "input_path",
                "source_url",
                "sha256",
                "enabled",
                "expected_bytes",
                "expected_md5",
            ]
            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "dataset_id": "public-r1",
                        "public_accession": accession,
                        "input_kind": "local_required",
                        "input_path": str(target),
                        "source_url": "landing-page",
                        "sha256": "NA",
                        "enabled": "false",
                        "expected_bytes": str(len(data)),
                        "expected_md5": __import__("hashlib").md5(data).hexdigest(),
                    }
                )
            original = manifest.read_bytes()
            output_manifest = root / "staged.tsv"
            audit = ENCODE_STAGER.stage(
                manifest,
                metadata_dir,
                root / "audit.json",
                repo=root,
                accessions={accession},
                output_manifest=output_manifest,
                download=True,
                refresh_metadata=False,
                curl_program="curl",
            )
            self.assertEqual(manifest.read_bytes(), original)
            with output_manifest.open(newline="", encoding="utf-8") as handle:
                staged = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(staged["enabled"], "true")
            self.assertEqual(staged["input_kind"], "local")
            self.assertEqual(
                staged["sha256"], __import__("hashlib").sha256(data).hexdigest()
            )
            self.assertFalse(audit["files"][0]["downloaded"])

            with manifest.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "dataset_id": "public-r1",
                        "public_accession": accession,
                        "input_kind": "local_required",
                        "input_path": str(target),
                        "source_url": "landing-page",
                        "sha256": "NA",
                        "enabled": "false",
                        "expected_bytes": str(len(data) + 1),
                        "expected_md5": metadata["md5sum"],
                    }
                )
            with self.assertRaisesRegex(ValueError, "expected_bytes"):
                ENCODE_STAGER.stage(
                    manifest,
                    metadata_dir,
                    root / "mismatch-audit.json",
                    repo=root,
                    accessions={accession},
                    output_manifest=None,
                    download=False,
                    refresh_metadata=False,
                    curl_program="curl",
                )

            gtf_accession = "ENCFFTESTGTF"
            gtf_metadata = dict(
                metadata,
                accession=gtf_accession,
                file_format="gtf",
                href=f"/files/{gtf_accession}/@@download/{gtf_accession}.gtf.gz",
            )
            validated = ENCODE_STAGER.validate_metadata(
                gtf_accession,
                gtf_metadata,
                expected_file_format="gtf",
            )
            self.assertEqual(validated["file_format"], "gtf")
            with self.assertRaisesRegex(ValueError, "expected 'fastq'"):
                ENCODE_STAGER.validate_metadata(gtf_accession, gtf_metadata)

    def test_smoke_generation_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "first"
            second = root / "second"
            GENERATOR.generate(
                first, seed=20260730, reads_per_transcript=4, substitution_rate=0.003
            )
            GENERATOR.generate(
                second, seed=20260730, reads_per_transcript=4, substitution_rate=0.003
            )
            first_files = sorted(path.name for path in first.iterdir() if path.is_file())
            self.assertEqual(first_files, sorted(path.name for path in second.iterdir()))
            for name in first_files:
                self.assertEqual((first / name).read_bytes(), (second / name).read_bytes())

    def test_encode_talon_extractor_follows_selected_transcript_membership(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "merged.gtf"
            source.write_text(
                'chr1\tHAVANA\ttranscript\t1\t30\t.\t+\t.\t'
                'gene_id "known"; transcript_id "known-tx";\n'
                'chr1\tHAVANA\texon\t1\t30\t.\t+\t.\t'
                'gene_id "known"; transcript_id "known-tx";\n'
                'chr1\tTALON\ttranscript\t101\t230\t.\t+\t.\t'
                'gene_id "known"; transcript_id "novel-tx"; talon_transcript "7";\n'
                'chr1\tHAVANA\texon\t101\t130\t.\t+\t.\t'
                'gene_id "known"; transcript_id "novel-tx"; talon_exon "10";\n'
                'chr1\tTALON\texon\t201\t230\t.\t+\t.\t'
                'gene_id "known"; transcript_id "novel-tx"; talon_exon "11";\n',
                encoding="utf-8",
            )
            first = root / "first.gtf.gz"
            second = root / "second.gtf.gz"
            record = ENCODE_TALON_EXTRACTOR.extract(source, first)
            ENCODE_TALON_EXTRACTOR.extract(source, second)
            self.assertEqual(record["selected_transcripts"], 1)
            self.assertEqual(record["selected_genes"], 1)
            self.assertEqual(record["selected_lines"], 3)
            self.assertEqual(record["selected_exon_lines"], 2)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            import gzip

            with gzip.open(first, "rt", encoding="utf-8") as handle:
                output = handle.read()
            self.assertNotIn("known-tx", output)
            self.assertEqual(output.count("novel-tx"), 3)
            self.assertIn("\tHAVANA\texon\t", output)

    def test_gtf_identity_view_audits_multivalued_nonidentity_attributes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.gtf"
            source.write_text(
                'chr1\tHAVANA\ttranscript\t101\t230\t.\t+\t.\t'
                'gene_id "g1"; transcript_id "tx1"; '
                'ont "PGO:1"; ont "PGO:2";\n'
                'chr1\tHAVANA\texon\t101\t130\t.\t+\t.\t'
                'gene_id "g1"; transcript_id "tx1"; '
                'tag "basic"; tag "MANE";\n'
                'chr1\tHAVANA\texon\t201\t230\t.\t+\t.\t'
                'gene_id "g1"; transcript_id "tx1";\n',
                encoding="utf-8",
            )
            first = root / "first.gtf.gz"
            second = root / "second.gtf.gz"
            record = GTF_IDENTITY_VIEW.project(source, first, repo=ROOT)
            GTF_IDENTITY_VIEW.project(source, second, repo=ROOT)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(record["duplicate_attribute_occurrences"], 2)
            self.assertEqual(
                record["duplicate_attribute_occurrences_by_key"],
                {"ont": 1, "tag": 1},
            )
            import gzip

            with gzip.open(first, "rt", encoding="utf-8") as handle:
                output = handle.read()
            self.assertNotIn("PGO:", output)
            self.assertNotIn('tag "', output)
            self.assertEqual(output.count('gene_id "g1";'), 3)
            self.assertEqual(output.count('transcript_id "tx1";'), 3)

            conflict = root / "conflict.gtf"
            conflict.write_text(
                'chr1\tx\texon\t1\t2\t.\t+\t.\t'
                'gene_id "g1"; gene_id "g2"; transcript_id "tx";\n',
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "conflicting duplicate gene_id"):
                GTF_IDENTITY_VIEW.project(
                    conflict, root / "conflict-view.gtf.gz", repo=ROOT
                )

    def test_manifest_validation_allows_declared_generated_outputs(self):
        with tempfile.TemporaryDirectory() as temporary:
            outdir = Path(temporary) / "work"
            result = VALIDATOR.validate(
                ROOT / "workflows/publication_benchmark/config/datasets.smoke.tsv",
                ROOT / "workflows/publication_benchmark/config/references.smoke.tsv",
                repo=ROOT,
                outdir=outdir,
                allow_deferred_generated=True,
            )
            self.assertEqual(result["status"], "valid")
            self.assertEqual(len(result["enabled_datasets"]), 2)
            self.assertEqual(len(result["enabled_reference_contexts"]), 2)

    def test_discovery_evaluator_recovers_truth_from_identical_gtf(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = root / "data"
            GENERATOR.generate(
                data, seed=20260730, reads_per_transcript=2, substitution_rate=0.0
            )
            summary_tsv = root / "summary.tsv"
            summary_json = root / "summary.json"
            subprocess.run(
                [
                    sys.executable,
                    SCRIPTS / "evaluate_discovery.py",
                    "--truth",
                    data / "truth.gtf",
                    "--input",
                    f"truth=smoke={data / 'truth.gtf'}",
                    "--output-tsv",
                    summary_tsv,
                    "--output-json",
                    summary_json,
                    "--repo",
                    ROOT,
                ],
                check=True,
            )
            with summary_tsv.open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["exact_precision"], "1.000000")
            self.assertEqual(row["exact_recall"], "1.000000")
            self.assertEqual(
                json.loads(summary_json.read_text())["truth_exact_forms"], 8
            )

    def test_identity_partition_metrics_count_false_merge_and_split_pairs(self):
        members = {
            ("caller", "sample", "a"): ("chr1", "+", ((1, 10), (21, 30))),
            ("caller", "sample", "b"): ("chr1", "+", ((1, 10), (21, 30))),
            ("caller", "sample", "c"): ("chr1", "+", ((2, 10), (21, 30))),
        }
        splices = {
            member: ("multi-exon", "chr1", "+", ((11, 20),))
            for member in members
        }
        metrics = IDENTITY_EVALUATOR._partition_metrics(
            {
                "merged": {
                    ("caller", "sample", "a"),
                    ("caller", "sample", "c"),
                },
                "split": {("caller", "sample", "b")},
            },
            members,
            splices,
        )
        self.assertEqual(metrics["clusters_mixing_exact_forms"], 1)
        self.assertEqual(metrics["clusters_merging_end_variants_only"], 1)
        self.assertEqual(metrics["exact_forms_split_across_clusters"], 1)
        self.assertEqual(metrics["observation_pairs_false_merged"], 1)
        self.assertEqual(metrics["observation_pairs_false_split"], 1)

    def test_identity_partition_metrics_report_explicit_tool_omissions(self):
        retained = ("caller", "sample", "retained")
        omitted = ("caller", "sample", "omitted")
        forms = {
            retained: ("chr1", "+", ((1, 10),)),
            omitted: ("chr1", "+", ((20, 30),)),
        }
        metrics = IDENTITY_EVALUATOR._partition_metrics(
            {"retained": {retained}},
            forms,
            forms,
            allow_missing=True,
            bootstrap_replicates=10,
        )
        self.assertEqual(metrics["observations_expected"], 2)
        self.assertEqual(metrics["observations_emitted"], 1)
        self.assertEqual(metrics["observations_missing"], 1)

    def test_order_comparison_separates_partition_and_identifier_changes(self):
        left = ("caller", "sample", "a")
        right = ("caller", "sample", "b")
        comparison = IDENTITY_EVALUATOR._order_comparison(
            {"1": {left, right}},
            {"9": {left, right}},
        )
        self.assertEqual(comparison["pair_relation_changes"], 0)
        self.assertEqual(comparison["cluster_identifier_changes"], 2)

    def test_order_comparison_counts_relation_changes_without_pair_expansion(self):
        first = ("caller", "sample", "a")
        second = ("caller", "sample", "b")
        third = ("caller", "sample", "c")
        comparison = IDENTITY_EVALUATOR._order_comparison(
            {"left": {first, second}, "single": {third}},
            {"single": {first}, "right": {second, third}},
        )
        self.assertEqual(comparison["pair_relation_changes"], 2)

    def test_txid_order_comparison_separates_exact_and_locus_accessions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            header = (
                "sample\ttool\toriginal_transcript_id\ttxid_form\ttxid_sc\t"
                "txid_transcript_id\ttxid_gene_id\tclassification\t"
                "gene_candidates\n"
            )
            sorted_path = root / "sorted.tsv"
            shuffled_path = root / "shuffled.tsv"
            sorted_path.write_text(
                header
                + "s\tt\ta\tform-a\tsc-a\tform-a\ttxid:GL1.000001\tnew_locus\t\n"
                + "s\tt\tb\tform-b\tsc-b\tform-b\ttxid:GL1.000001\tnew_locus\t\n",
                encoding="utf-8",
            )
            shuffled_path.write_text(
                header
                + "s\tt\ta\tform-a\tsc-a\tform-a\ttxid:GL1.000002\tnew_locus\t\n"
                + "s\tt\tb\tform-b\tsc-b\tform-b\ttxid:GL1.000003\tnew_locus\t\n",
                encoding="utf-8",
            )
            comparison = IDENTITY_EVALUATOR._txid_order_comparison(
                sorted_path, shuffled_path
            )
        self.assertEqual(comparison["field_changes"]["txid_form"], 0)
        self.assertEqual(comparison["field_changes"]["txid_gene_id"], 2)
        self.assertEqual(
            comparison["exact_form_partition"]["pair_relation_changes"], 0
        )
        self.assertEqual(
            comparison["gene_locus_partition"]["pair_relation_changes"], 1
        )

    def test_order_comparison_separates_omissions_from_shared_relations(self):
        shared = ("caller", "sample", "shared")
        sorted_only = ("caller", "sample", "sorted-only")
        shuffled_only = ("caller", "sample", "shuffled-only")
        comparison = IDENTITY_EVALUATOR._order_comparison(
            {"sorted": {shared, sorted_only}},
            {"shuffled": {shared, shuffled_only}},
        )
        self.assertEqual(comparison["members_shared"], 1)
        self.assertEqual(comparison["members_only_sorted"], 1)
        self.assertEqual(comparison["members_only_shuffled"], 1)
        self.assertEqual(comparison["pair_relation_changes"], 0)

    def test_partition_metrics_can_disclose_comparator_omissions(self):
        first = ("caller", "sample", "a")
        omitted = ("caller", "sample", "b")
        forms = {
            first: ("chr1", "+", ((1, 10),)),
            omitted: ("chr1", "+", ((20, 30),)),
        }
        metrics = IDENTITY_EVALUATOR._partition_metrics(
            {"reported": {first}},
            forms,
            forms,
            allow_missing=True,
        )
        self.assertEqual(metrics["input_observations"], 2)
        self.assertEqual(metrics["observations"], 1)
        self.assertEqual(metrics["omitted_observations"], 1)
        self.assertEqual(metrics["matrix_rows"], 1)
        self.assertEqual(metrics["matrix_columns"], 1)
        self.assertEqual(metrics["matrix_occupied_cells"], 1)

    def test_talon_partition_metrics_count_truth_pair_errors(self):
        first = ("sample", "read-a")
        second = ("sample", "read-b")
        third = ("sample", "read-c")
        metrics = TALON_EVALUATOR._partition_metrics(
            {"merged": {first, third}, "split": {second}},
            {first: "truth-a", second: "truth-a", third: "truth-b"},
        )
        self.assertEqual(metrics["clusters_mixing_truth_groups"], 1)
        self.assertEqual(metrics["truth_groups_split_across_clusters"], 1)
        self.assertEqual(metrics["observation_pairs_false_merged"], 1)
        self.assertEqual(metrics["observation_pairs_false_split"], 1)

    def test_talon_truth_subset_reports_unlabelled_and_missing_reads(self):
        labelled = ("sample", "labelled")
        missing = ("sample", "missing")
        unlabelled = ("sample", "human")
        restricted, ignored = TALON_EVALUATOR._truth_subset(
            {"1": {labelled, unlabelled}},
            {labelled, missing},
        )
        self.assertEqual(ignored, 1)
        metrics = TALON_EVALUATOR._partition_metrics(
            restricted,
            {labelled: "sirv-1", missing: "sirv-2"},
            allow_missing_truth=True,
        )
        self.assertEqual(metrics["observations"], 1)
        self.assertEqual(metrics["truth_observations_requested"], 2)
        self.assertEqual(metrics["truth_observations_missing"], 1)

    def test_sirv_transcriptome_and_unique_assignment_are_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            reference = root / "sirv.fa"
            reference.write_text(
                ">SIRV1\nAACCGGTTAACCGGTT\n>SIRV2\nAAAACCCCGGGGTTTT\n",
                encoding="utf-8",
            )
            truth_gtf = root / "sirv.gtf"
            truth_gtf.write_text(
                'SIRV1\tt\texon\t1\t4\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
                'SIRV1\tt\texon\t9\t12\t.\t+\t.\tgene_id "g1"; transcript_id "tx1";\n'
                'SIRV2\tt\texon\t1\t4\t.\t-\t.\tgene_id "g2"; transcript_id "tx2";\n'
                'SIRV2\tt\texon\t9\t12\t.\t-\t.\tgene_id "g2"; transcript_id "tx2";\n',
                encoding="utf-8",
            )
            transcript_fasta = root / "transcripts.fa"
            transcript_truth = root / "truth.tsv"
            transcripts = SIRV_TRUTH.build_transcriptome(
                reference, truth_gtf, transcript_fasta, transcript_truth
            )
            self.assertEqual([value.truth_id for value in transcripts], ["tx1", "tx2"])
            self.assertEqual(
                transcript_fasta.read_text(encoding="utf-8"),
                ">tx1\nAACCAACC\n>tx2\nCCCCTTTT\n",
            )
            paf = root / "reads.paf"
            paf.write_text(
                "unique\t100\t0\t100\t+\ttx1\t100\t0\t100\t98\t100\t60\tAS:i:200\n"
                "unique\t100\t0\t100\t+\ttx2\t100\t0\t100\t90\t100\t20\tAS:i:100\n"
                "ambiguous\t100\t0\t100\t+\ttx1\t100\t0\t100\t98\t100\t10\tAS:i:200\n"
                "ambiguous\t100\t0\t100\t+\ttx2\t100\t0\t100\t98\t100\t10\tAS:i:180\n",
                encoding="utf-8",
            )
            assignments, counts = SIRV_TRUTH.assign_candidates(
                SIRV_TRUTH.parse_paf(paf),
                min_identity=0.95,
                min_query_coverage=0.90,
                min_target_coverage=0.90,
                min_score_margin=50,
            )
            self.assertEqual(assignments, [("unique", "tx1")])
            self.assertEqual(counts["accepted"], 1)
            self.assertEqual(counts["ambiguous_score_margin"], 1)

    def test_public_v49_sirv_gtf_build_is_byte_deterministic(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            annotation = root / "v49.gtf.gz"
            sirv = root / "sirv.gtf.gz"
            import gzip

            with gzip.open(annotation, "wt", encoding="utf-8") as handle:
                handle.write("# v49\nchr1\tt\texon\t1\t2\t.\t+\t.\tgene_id \"g\";\n")
            with gzip.open(sirv, "wt", encoding="utf-8") as handle:
                handle.write("SIRV1\tt\texon\t1\t2\t.\t+\t.\tgene_id \"s\";\n")
            first = root / "first.gtf.gz"
            second = root / "second.gtf.gz"
            REFERENCE_STAGER.build_combined_gtf(annotation, sirv, first)
            REFERENCE_STAGER.build_combined_gtf(annotation, sirv, second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            with gzip.open(first, "rt", encoding="utf-8") as handle:
                self.assertEqual(
                    handle.read(),
                    "# v49\n"
                    'chr1\tt\texon\t1\t2\t.\t+\t.\tgene_id "g";\n'
                    'SIRV1\tt\texon\t1\t2\t.\t+\t.\tgene_id "s";\n',
                )


if __name__ == "__main__":
    unittest.main()
