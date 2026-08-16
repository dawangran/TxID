from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from txid.errors import ReferenceError, TxIDError
from txid.parser import parse_annotation
from txid.registry import Registry
from txid.workflow import add_annotation_context, batch_import, import_file, multi_import

from tests.helpers import REFERENCE_GTF, make_registry


class WorkflowTests(unittest.TestCase):
    def test_multi_import_matches_equivalent_manifest_bytes(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            manifest_root = root / "manifest-case"
            direct_root = root / "direct-case"
            manifest_root.mkdir()
            direct_root.mkdir()
            manifest_db, manifest_first = make_registry(manifest_root)
            direct_db, direct_first = make_registry(direct_root)
            manifest_second = manifest_root / "second.gtf"
            direct_second = direct_root / "second.gtf"
            for first, second in (
                (manifest_first, manifest_second),
                (direct_first, direct_second),
            ):
                second.write_text(
                    first.read_text(encoding="utf-8").replace("caller", "second"),
                    encoding="utf-8",
                )

            manifest = manifest_root / "manifest.tsv"
            manifest.write_text(
                "input\tsample\ttool\tannotation_name\n"
                f"{manifest_second}\ts2\tIsoQuant\tref-v1\n"
                f"{manifest_first}\ts1\tIsoQuant\tref-v1\n",
                encoding="utf-8",
            )
            manifest_results = batch_import(
                manifest_db,
                manifest=manifest,
                output_dir=manifest_root / "output",
                fuzzy_splice_tolerance=5,
                fuzzy_end_tolerance=15,
            )
            direct_results = multi_import(
                direct_db,
                inputs=[direct_second, direct_first],
                samples=["s2", "s1"],
                tool="IsoQuant",
                annotation_name="ref-v1",
                output_dir=direct_root / "output",
                fuzzy_splice_tolerance=5,
                fuzzy_end_tolerance=15,
            )

            self.assertEqual(
                [item["sample"] for item in manifest_results],
                [item["sample"] for item in direct_results],
            )
            for manifest_result, direct_result in zip(
                manifest_results, direct_results
            ):
                for key in ("output_gtf", "mapping"):
                    manifest_path = Path(manifest_result[key])
                    direct_path = Path(direct_result[key])
                    self.assertEqual(manifest_path.name, direct_path.name)
                    self.assertEqual(manifest_path.read_bytes(), direct_path.read_bytes())
            with Registry(direct_db, read_only=True) as registry:
                self.assertEqual(registry.validate(), [])

    def test_multi_import_requires_one_explicit_sample_per_input(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, caller = make_registry(root)
            with self.assertRaisesRegex(TxIDError, "exactly one sample"):
                multi_import(
                    database,
                    inputs=[caller, caller],
                    samples=["only-one"],
                    tool="IsoQuant",
                    annotation_name="ref-v1",
                    output_dir=root / "output",
                )

    def test_batch_fuzzy_clustering_runs_once_and_updates_every_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, first_input = make_registry(root)
            second_input = root / "caller-second.gtf"
            second_input.write_text(
                first_input.read_text(encoding="utf-8").replace("caller", "second"),
                encoding="utf-8",
            )
            manifest = root / "manifest.tsv"
            manifest.write_text(
                "input\tsample\ttool\tannotation_name\n"
                f"{first_input}\ts1\tIsoQuant\tref-v1\n"
                f"{second_input}\ts2\tFLAIR\tref-v1\n",
                encoding="utf-8",
            )

            original_run_fuzzy = Registry.run_fuzzy
            fuzzy_calls: list[tuple[int, int]] = []

            def counted_run_fuzzy(registry, *, splice_tolerance, end_tolerance, software_version):
                fuzzy_calls.append((splice_tolerance, end_tolerance))
                return original_run_fuzzy(
                    registry,
                    splice_tolerance=splice_tolerance,
                    end_tolerance=end_tolerance,
                    software_version=software_version,
                )

            with mock.patch.object(Registry, "run_fuzzy", new=counted_run_fuzzy):
                results = batch_import(
                    database,
                    manifest=manifest,
                    output_dir=root / "batch-output",
                    fuzzy_splice_tolerance=5,
                    fuzzy_end_tolerance=15,
                )

            self.assertEqual(fuzzy_calls, [(5, 15)])
            self.assertEqual(len(results), 2)
            form_to_cluster: list[dict[str, str]] = []
            for result in results:
                with Path(result["mapping"]).open(
                    encoding="utf-8", newline=""
                ) as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                self.assertTrue(rows)
                self.assertTrue(all(row["fuzzy_cluster"] for row in rows))
                form_to_cluster.append(
                    {row["txid_form"]: row["fuzzy_cluster"] for row in rows}
                )
            self.assertEqual(form_to_cluster[0], form_to_cluster[1])

    def test_locus_cache_is_discarded_after_transaction_rollback(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            with Registry(database) as registry:
                with self.assertRaisesRegex(RuntimeError, "abort staged locus"):
                    with registry.transaction():
                        registry.allocate_locus(
                            contig="chr1",
                            strand="+",
                            start=700,
                            end=800,
                            status="new_locus",
                            gene_candidates=(),
                        )
                        raise RuntimeError("abort staged locus")
                with registry.transaction():
                    accession, public_id = registry.allocate_locus(
                        contig="chr1",
                        strand="+",
                        start=700,
                        end=800,
                        status="new_locus",
                        gene_candidates=(),
                    )
                self.assertEqual(accession, 1)
                self.assertEqual(public_id, "txid:GL1.000001")
                self.assertEqual(
                    registry.connection.execute(
                        "SELECT COUNT(*) FROM locus"
                    ).fetchone()[0],
                    1,
                )

    def test_repeated_gtf_attributes_survive_registry_and_rewrite(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            caller = root / "repeated-tags.gtf"
            output = root / "rewritten.gtf"
            mapping = root / "mapping.tsv"
            caller.write_text(
                'chr1\tcaller\ttranscript\t100\t399\t.\t+\t.\tgene_id "up"; transcript_id "tagged"; tag "basic"; tag "MANE_Select";\n'
                'chr1\tcaller\texon\t100\t199\t.\t+\t.\tgene_id "up"; transcript_id "tagged"; tag "basic"; tag "Ensembl_canonical";\n'
                'chr1\tcaller\texon\t300\t399\t.\t+\t.\tgene_id "up"; transcript_id "tagged"; tag "basic"; tag "Ensembl_canonical";\n',
                encoding="utf-8",
            )

            first = import_file(
                database,
                input_path=caller,
                sample="sample-tags",
                tool="GENCODE-like",
                annotation_name="ref-v1",
                output_gtf=output,
                mapping_path=mapping,
            )
            self.assertTrue(first.created)
            first_bytes = output.read_bytes()
            text = first_bytes.decode("utf-8")
            self.assertEqual(text.count('tag "basic";'), 3)
            self.assertEqual(text.count('tag "MANE_Select";'), 1)
            self.assertEqual(text.count('tag "Ensembl_canonical";'), 2)

            second = import_file(
                database,
                input_path=caller,
                sample="sample-tags",
                tool="GENCODE-like",
                annotation_name="ref-v1",
                output_gtf=output,
                mapping_path=mapping,
            )
            self.assertFalse(second.created)
            self.assertEqual(first_bytes, output.read_bytes())
            reparsed = parse_annotation(output)
            self.assertEqual(
                [value for key, value in reparsed[0].attributes if key == "tag"],
                ["basic", "MANE_Select"],
            )

    def test_end_to_end_classification_retry_and_validation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, caller = make_registry(root)
            output = root / "rewritten.gtf"
            mapping = root / "mapping.tsv"
            first = import_file(
                database,
                input_path=caller,
                sample="sample-A",
                tool="IsoQuant",
                annotation_name="ref-v1",
                output_gtf=output,
                mapping_path=mapping,
                fuzzy_splice_tolerance=5,
                fuzzy_end_tolerance=15,
            )
            self.assertTrue(first.created)
            with mapping.open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 9)
            by_id = {row["original_transcript_id"]: row for row in rows}
            self.assertEqual(by_id["alt_known"]["classification"], "known")
            self.assertEqual(by_id["alt_known"]["txid_transcript_id"], "tx_ref1")
            self.assertEqual(by_id["alt_known"]["txid_gene_id"], "g1")
            self.assertEqual(by_id["alt_minus"]["txid_transcript_id"], "tx_ref3")
            self.assertEqual(by_id["end_variant"]["classification"], "novel_in_known_gene")
            self.assertEqual(by_id["retained_intron"]["txid_gene_id"], "g2")
            self.assertEqual(by_id["ambiguous_AB"]["classification"], "ambiguous_gene")
            self.assertEqual(by_id["readthrough_1"]["classification"], "ambiguous_gene")
            self.assertEqual(by_id["antisense_single"]["classification"], "new_locus")
            self.assertEqual(by_id["new_chr2"]["classification"], "new_locus")
            self.assertEqual(by_id["alt_known"]["txid_sc"], by_id["end_variant"]["txid_sc"])
            self.assertNotEqual(by_id["alt_known"]["txid_form"], by_id["end_variant"]["txid_form"])
            self.assertTrue(all(row["fuzzy_cluster"].startswith("txid:FC1.") for row in rows))
            first_gtf = output.read_bytes()
            first_mapping = mapping.read_bytes()

            second = import_file(
                database,
                input_path=caller,
                sample="sample-A",
                tool="IsoQuant",
                annotation_name="ref-v1",
                output_gtf=output,
                mapping_path=mapping,
                fuzzy_splice_tolerance=5,
                fuzzy_end_tolerance=15,
            )
            self.assertFalse(second.created)
            self.assertEqual(first_gtf, output.read_bytes())
            self.assertEqual(first_mapping, mapping.read_bytes())
            with Registry(database, read_only=True) as registry:
                self.assertEqual(registry.validate(), [])
                summary = registry.summary()
            self.assertEqual(summary["imports"], 1)
            self.assertEqual(summary["observations"], 9)
            reparsed = parse_annotation(output)
            self.assertEqual(len(reparsed), 9)

    def test_failed_reference_validation_leaves_no_partial_import(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            invalid = root / "unknown.gtf"
            invalid.write_text(
                'unknown\tt\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t";\n',
                encoding="utf-8",
            )
            with self.assertRaises(ReferenceError):
                import_file(
                    database,
                    input_path=invalid,
                    sample="bad",
                    tool="tool",
                    annotation_name="ref-v1",
                    output_gtf=root / "bad.gtf",
                    mapping_path=root / "bad.tsv",
                )
            with Registry(database, read_only=True) as registry:
                self.assertEqual(registry.summary()["imports"], 0)

            transaction_database = root / "transaction.sqlite"
            with Registry(transaction_database) as registry:
                registry.connection.execute("CREATE TABLE item (value TEXT)")
                with self.assertRaisesRegex(RuntimeError, "original failure"):
                    with registry.transaction():
                        registry.connection.execute(
                            "INSERT INTO item(value) VALUES ('staged')"
                        )
                        registry.connection.execute("ROLLBACK")
                        raise RuntimeError("original failure")
                self.assertEqual(
                    registry.connection.execute(
                        "SELECT COUNT(*) FROM item"
                    ).fetchone()[0],
                    0,
                )

    def test_second_annotation_context_retains_structure_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            release2 = root / "release2.gtf"
            release2.write_text(
                REFERENCE_GTF.replace('gene_id "g1"', 'gene_id "g1_release2"').replace(
                    'transcript_id "tx_ref1"', 'transcript_id "tx_release2"'
                ),
                encoding="utf-8",
            )
            add_annotation_context(
                database,
                annotation=release2,
                annotation_name="ref-v2",
            )
            exact = root / "exact.gtf"
            exact.write_text(
                'chr1\tt\texon\t300\t399\t.\t+\t.\tgene_id "x"; transcript_id "x";\n'
                'chr1\tt\texon\t100\t199\t.\t+\t.\tgene_id "x"; transcript_id "x";\n',
                encoding="utf-8",
            )
            v1 = import_file(
                database,
                input_path=exact,
                sample="s1",
                tool="t",
                annotation_name="ref-v1",
                output_gtf=root / "v1.gtf",
                mapping_path=root / "v1.tsv",
            )
            v2 = import_file(
                database,
                input_path=exact,
                sample="s2",
                tool="t",
                annotation_name="ref-v2",
                output_gtf=root / "v2.gtf",
                mapping_path=root / "v2.tsv",
            )
            self.assertEqual(v1.assignments[0].identities.form.public_id, v2.assignments[0].identities.form.public_id)
            self.assertEqual(v1.assignments[0].output_transcript_id, "tx_ref1")
            self.assertEqual(v2.assignments[0].output_transcript_id, "tx_release2")
            self.assertEqual(v2.assignments[0].output_gene_id, "g1_release2")


if __name__ == "__main__":
    unittest.main()
