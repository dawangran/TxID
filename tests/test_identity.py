from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from txid.errors import IdentityCollisionError
from txid.identity import identify
from txid.models import ContigRecord, Exon, SequenceCollection, TranscriptModel
from txid.registry import Registry

ASSEMBLY = "sha256:" + "0" * 64


class IdentityTests(unittest.TestCase):
    def model(self, exons, strand="+"):
        return TranscriptModel("chr1", strand, tuple(Exon(*pair) for pair in exons), "t1")

    def test_coordinate_semantics_on_both_strands(self):
        plus = self.model([(100, 199), (300, 399)], "+")
        minus = self.model([(100, 199), (300, 399)], "-")
        self.assertEqual(plus.introns, ((200, 299),))
        self.assertEqual((plus.tss, plus.tes), (100, 399))
        self.assertEqual(minus.introns, ((200, 299),))
        self.assertEqual((minus.tss, minus.tes), (399, 100))
        self.assertNotEqual(identify(plus, ASSEMBLY).form.public_id, identify(minus, ASSEMBLY).form.public_id)

    def test_exon_order_is_irrelevant(self):
        left = self.model([(300, 399), (100, 199)])
        right = self.model([(100, 199), (300, 399)])
        self.assertEqual(identify(left, ASSEMBLY), identify(right, ASSEMBLY))

    def test_one_base_splice_change_changes_both_exact_ids(self):
        a = identify(self.model([(100, 199), (300, 399)]), ASSEMBLY)
        b = identify(self.model([(100, 200), (300, 399)]), ASSEMBLY)
        self.assertNotEqual(a.splice_chain.public_id, b.splice_chain.public_id)
        self.assertNotEqual(a.form.public_id, b.form.public_id)

    def test_end_variant_retains_splice_chain_but_changes_form(self):
        a = identify(self.model([(100, 199), (300, 399)]), ASSEMBLY)
        b = identify(self.model([(110, 199), (300, 410)]), ASSEMBLY)
        self.assertEqual(a.splice_chain.public_id, b.splice_chain.public_id)
        self.assertNotEqual(a.form.public_id, b.form.public_id)

    def test_single_exon_uses_se_family(self):
        identity = identify(self.model([(100, 399)]), ASSEMBLY)
        self.assertIsNone(identity.splice_chain)
        self.assertTrue(identity.form.public_id.startswith("txid:SE1."))

    def test_public_digest_collision_is_a_hard_error(self):
        reference = SequenceCollection(
            "test",
            ASSEMBLY,
            (ContigRecord("chr1", 1000, "f" * 64),),
        )
        with tempfile.TemporaryDirectory() as temporary:
            database = Path(temporary) / "collision.sqlite"
            with Registry.create(
                database,
                reference,
                fasta_checksum="sha256:" + "e" * 64,
                software_version="test",
            ) as registry:
                first = identify(
                    self.model([(10, 20)]),
                    ASSEMBLY,
                    digest_function=lambda _: "a" + "1" * 63,
                    public_digest_length=1,
                ).form
                second = identify(
                    self.model([(30, 40)]),
                    ASSEMBLY,
                    digest_function=lambda _: "a" + "2" * 63,
                    public_digest_length=1,
                ).form
                with registry.transaction():
                    registry.register_structural(first)
                with self.assertRaises(IdentityCollisionError):
                    with registry.transaction():
                        registry.register_structural(second)


if __name__ == "__main__":
    unittest.main()

