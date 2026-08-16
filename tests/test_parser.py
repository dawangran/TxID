from __future__ import annotations

import unittest

from txid.errors import AnnotationParseError
from txid.parser import parse_annotation_lines, parse_gtf_attributes


class ParserTests(unittest.TestCase):
    def test_gtf_without_meta_features_and_escaped_attribute(self):
        models = parse_annotation_lines(
            [
                'chr1\tt\texon\t20\t30\t.\t+\t.\tgene_id "g"; transcript_id "t"; note "a\\\"b";\n',
                'chr1\tt\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t"; note "a\\\"b";\n',
            ],
            fmt="gtf",
        )
        self.assertEqual(len(models), 1)
        self.assertEqual([(e.start, e.end) for e in models[0].exons], [(1, 10), (20, 30)])
        self.assertEqual(dict(models[0].exons[0].attributes)["note"], 'a"b')

    def test_gff3_with_and_without_transcript_meta_feature(self):
        models = parse_annotation_lines(
            [
                "##gff-version 3\n",
                "chr1\tr\tgene\t1\t100\t.\t+\t.\tID=g1\n",
                "chr1\tr\tmRNA\t1\t100\t.\t+\t.\tID=t1;Parent=g1;Name=hello%20world\n",
                "chr1\tr\texon\t1\t10\t.\t+\t.\tID=e1;Parent=t1\n",
                "chr1\tr\texon\t20\t30\t.\t+\t.\tID=e2;Parent=t1,ghost\n",
                "chr1\tr\texon\t40\t50\t.\t+\t.\tParent=ghost\n",
            ],
            fmt="gff3",
        )
        by_id = {model.original_transcript_id: model for model in models}
        self.assertEqual(by_id["t1"].original_gene_id, "g1")
        self.assertEqual(dict(by_id["t1"].attributes)["Name"], "hello world")
        self.assertEqual(len(by_id["ghost"].exons), 2)
        self.assertIsNone(by_id["ghost"].original_gene_id)

    def test_repeated_nonidentity_gtf_attributes_are_preserved(self):
        models = parse_annotation_lines(
            [
                'chr1\tt\ttranscript\t1\t30\t.\t+\t.\tgene_id "g"; transcript_id "t"; tag "basic"; tag "MANE_Select";\n',
                'chr1\tt\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t"; tag "basic"; tag "Ensembl_canonical";\n',
                'chr1\tt\texon\t20\t30\t.\t+\t.\tgene_id "g"; transcript_id "t";\n',
            ],
            fmt="gtf",
        )
        transcript_tags = [
            value for key, value in models[0].attributes if key == "tag"
        ]
        exon_tags = [
            value for key, value in models[0].exons[0].attributes if key == "tag"
        ]
        self.assertEqual(transcript_tags, ["basic", "MANE_Select"])
        self.assertEqual(exon_tags, ["basic", "Ensembl_canonical"])

    def test_gtf_attribute_fast_path_matches_quoted_delimiter_semantics(self):
        attributes = parse_gtf_attributes(
            'gene_id "g"; transcript_id "t"; '
            'description "alpha; beta"; note "escaped\\;semicolon"; '
            'tag "basic"; tag "MANE_Select";',
            7,
        )
        self.assertEqual(
            attributes,
            (
                ("description", "alpha; beta"),
                ("gene_id", "g"),
                ("note", "escaped;semicolon"),
                ("tag", "basic"),
                ("tag", "MANE_Select"),
                ("transcript_id", "t"),
            ),
        )

    def test_malformed_and_duplicate_identity_attributes_are_rejected(self):
        with self.assertRaises(AnnotationParseError):
            parse_annotation_lines(["chr1\tt\texon\t1\t2\t.\t+\t.\tgene_id \"g\";\n"], fmt="gtf")
        with self.assertRaises(AnnotationParseError):
            parse_gtf_attributes('gene_id "a"; gene_id "b";', 1)
        with self.assertRaises(AnnotationParseError):
            parse_annotation_lines(
                ['chr1\tt\texon\t1\t2\t.\t+\t.\tgene_id "g"; transcript_id "t";\n',
                 'chr2\tt\texon\t4\t5\t.\t+\t.\tgene_id "g"; transcript_id "t";\n'],
                fmt="gtf",
            )


if __name__ == "__main__":
    unittest.main()
