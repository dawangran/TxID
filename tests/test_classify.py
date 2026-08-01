from __future__ import annotations

import unittest

from txid.classify import ReferenceEntry, build_reference_index, classify
from txid.identity import identify
from txid.models import Exon, TranscriptModel


class ClassificationTests(unittest.TestCase):
    def test_gene_assignment_uses_aggregate_reference_gene_span(self):
        assembly = "sha256:" + "0" * 64
        candidate = TranscriptModel("chr1", "+", (Exon(220, 250),), "candidate")
        references = [
            ReferenceEntry("gene", "tx1", "chr1", "+", 100, 199, "form1", None),
            ReferenceEntry("gene", "tx2", "chr1", "+", 300, 399, "form2", None),
        ]
        result = classify(candidate, identify(candidate, assembly), references)
        self.assertEqual(result.label, "novel_in_known_gene")
        self.assertEqual(result.output_gene_id, "gene")

    def test_prebuilt_reference_index_preserves_exact_and_overlap_semantics(self):
        assembly = "sha256:" + "0" * 64
        known = TranscriptModel(
            "chr1",
            "+",
            (Exon(100, 150), Exon(201, 250)),
            "preferred",
        )
        known_identity = identify(known, assembly)
        references = [
            ReferenceEntry(
                "gene-b",
                "alternate",
                "chr1",
                "+",
                100,
                250,
                known_identity.form.public_id,
                known_identity.splice_chain.public_id,
            ),
            ReferenceEntry(
                "gene-a",
                "preferred",
                "chr1",
                "+",
                100,
                250,
                known_identity.form.public_id,
                known_identity.splice_chain.public_id,
            ),
            ReferenceEntry(
                "gene-c",
                "remote",
                "chr1",
                "+",
                1000,
                1200,
                "remote-form",
                None,
            ),
        ]
        index = build_reference_index(references)
        list_result = classify(known, known_identity, references)
        index_result = classify(known, known_identity, index)
        self.assertEqual(index_result, list_result)
        self.assertEqual(index_result.output_transcript_id, "preferred")
        self.assertEqual(index_result.gene_candidates, ("gene-a", "gene-b"))

        novel = TranscriptModel("chr1", "+", (Exon(175, 180),), "novel")
        list_result = classify(novel, identify(novel, assembly), references)
        index_result = classify(novel, identify(novel, assembly), index)
        self.assertEqual(index_result, list_result)
        self.assertEqual(index_result.label, "ambiguous_gene")
        self.assertEqual(index_result.gene_candidates, ("gene-a", "gene-b"))

        antisense = TranscriptModel("chr1", "-", (Exon(175, 180),), "antisense")
        index_result = classify(antisense, identify(antisense, assembly), index)
        self.assertEqual(index_result.label, "new_locus")


if __name__ == "__main__":
    unittest.main()
