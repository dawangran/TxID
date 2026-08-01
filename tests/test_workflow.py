from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from txid.errors import ReferenceError
from txid.parser import parse_annotation
from txid.registry import Registry
from txid.workflow import add_annotation_context, import_file

from tests.helpers import REFERENCE_GTF, make_registry


class WorkflowTests(unittest.TestCase):
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
