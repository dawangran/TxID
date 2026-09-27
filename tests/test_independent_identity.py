"""Checks for the separately coded GTF/identity oracle (no production imports)."""

from __future__ import annotations

import ast
import gzip
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from benchmarks.independent_identity import (
    IndependentIdentityError,
    compute_identity,
    read_gtf_models,
    sha256_file,
    verify_conformance,
)


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "benchmarks/independent_identity.py"
VECTORS = ROOT / "tests/conformance/exact-v1.json"
ASSEMBLY = "sha256:" + "0" * 64


def row(start=100, end=199, *, feature="exon", strand="+", contig="chr1",
        transcript="t", gene="g", extra=""):
    attrs = f'transcript_id "{transcript}";'
    if gene is not None:
        attrs = f'gene_id "{gene}"; ' + attrs
    return f"{contig}\tcaller\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs} {extra}\n"


class IndependentCanonicalTests(unittest.TestCase):
    def identify(self, exons, *, strand="+", contig="chr1", assembly=ASSEMBLY):
        return compute_identity(assembly_fingerprint=assembly, contig=contig,
                                strand=strand, exons=exons)

    def test_all_existing_golden_canonical_bytes_and_digests(self):
        report = verify_conformance(VECTORS)
        self.assertTrue(report["valid"])
        self.assertEqual(report["vectors"], 3)
        self.assertEqual(report["identities_checked"], 5)
        self.assertEqual(report["failures"], [])
        self.assertEqual(report["inputs"][0]["sha256"], sha256_file(VECTORS))

    def test_transcription_order_end_policy_and_permutations_on_both_strands(self):
        intervals = [(2, 9), (21, 30), (41, 70)]
        for strand, expected_gaps, tss, tes in (
            ("+", [[10, 20], [31, 40]], 2, 70),
            ("-", [[31, 40], [10, 20]], 70, 2),
        ):
            expected = self.identify(intervals, strand=strand)
            form = json.loads(expected["form"]["canonical_json"])
            self.assertEqual(form["introns"], expected_gaps)
            self.assertEqual((form["tss"], form["tes"]), (tss, tes))
            for permutation in itertools.permutations(intervals):
                with self.subTest(strand=strand, permutation=permutation):
                    self.assertEqual(self.identify(permutation, strand=strand), expected)

    def test_terminal_variants_and_one_base_splice_change(self):
        for strand in ("+", "-"):
            reference = self.identify([(100, 199), (300, 399)], strand=strand)
            ends = self.identify([(90, 199), (300, 420)], strand=strand)
            splice = self.identify([(100, 200), (300, 399)], strand=strand)
            self.assertEqual(reference["splice_chain"], ends["splice_chain"])
            self.assertNotEqual(reference["form"], ends["form"])
            self.assertNotEqual(reference["splice_chain"], splice["splice_chain"])

    def test_single_exon_and_reference_contexts(self):
        original = self.identify([(1, 1)])
        self.assertIsNone(original["splice_chain"])
        canonical = json.loads(original["form"]["canonical_json"])
        self.assertEqual(canonical["algorithm"], "SE1")
        self.assertEqual((canonical["start"], canonical["end"]), (1, 1))
        for change in ({"strand": "-"}, {"contig": "chr2"}, {"assembly": "sha256:" + "f" * 64}):
            self.assertNotEqual(self.identify([(1, 1)], **change)["form"], original["form"])
        self.assertNotEqual(self.identify([(1, 2)])["form"], original["form"])

    def test_unicode_canonicalization_and_digest_truncation(self):
        result = self.identify([(4, 8)], contig="染色体1")["form"]
        self.assertIn("染色体1", result["canonical_json"])
        self.assertNotIn("\\u", result["canonical_json"])
        independent_digest = hashlib.sha256(result["canonical_json"].encode()).hexdigest()
        self.assertEqual(result["full_digest"], independent_digest)
        self.assertEqual(result["public_id"], "txid:SE1." + independent_digest[:24])

    def test_invalid_coordinates_and_context_are_errors(self):
        invalid = [[], [(0, 5)], [(9, 3)], [(True, 4)], [(1.0, 4)],
                   [(1, 4, 6)], [(1, 4), (1, 4)], [(1, 4), (4, 9)], [(1, 4), (5, 9)]]
        for intervals in invalid:
            with self.subTest(intervals=intervals), self.assertRaises(IndependentIdentityError):
                self.identify(intervals)
        for fields in ({"strand": "."}, {"contig": ""}, {"contig": "chr 1"}, {"contig": "chr\x001"},
                       {"assembly": "GRCh38"}, {"assembly": "sha256:" + "A" * 64}):
            with self.subTest(fields=fields), self.assertRaises(IndependentIdentityError):
                self.identify([(1, 4)], **fields)

    def test_oracle_has_no_txid_import(self):
        imports = []
        for node in ast.walk(ast.parse(SCRIPT.read_text())):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append(node.module or "")
        self.assertFalse(any(name == "txid" or name.startswith("txid.") for name in imports))


class IndependentGTFTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def parse(self, text, **kwargs):
        path = self.root / "input.gtf"
        path.write_text(text, encoding="utf-8")
        return read_gtf_models(path, **kwargs)

    def test_interleaved_records_missing_transcript_features_and_negative_strand(self):
        text = row(300, 350, strand="-") + row(12, 22, transcript="u", gene=None)
        text += row(100, 150, strand="-") + row(100, 350, feature="transcript", strand="-")
        models = self.parse(text)
        self.assertEqual(list(models), ["t", "u"])
        self.assertEqual(models["t"].exons, ((100, 150), (300, 350)))
        self.assertIsNone(models["u"].gene_id)
        payload = json.loads(models["t"].identity(ASSEMBLY)["form"]["canonical_json"])
        self.assertEqual((payload["tss"], payload["tes"]), (350, 100))
        self.assertEqual(models, self.parse("".join(reversed(text.splitlines(keepends=True)))))

    def test_quoted_escaping_semicolons_and_repeated_nonidentity_attributes(self):
        text = 'chr1\ttool\texon\t5\t10\t.\t+\t.\tgene_id "g"; transcript_id "t;\\\"x\\\\z"; tag "a"; tag "b"; cov 1.5;\n'
        models = self.parse(text)
        self.assertEqual(list(models), ['t;"x\\z'])
        self.assertEqual(models['t;"x\\z'].exons, ((5, 10),))

    def test_gzip_magic_and_explicit_contig_alias(self):
        text = row(contig="1") + row(300, 399, contig="chr1")
        compressed = self.root / "compressed.data"
        compressed.write_bytes(gzip.compress(text.encode(), mtime=0))
        models = read_gtf_models(compressed, contig_aliases={"1": "chr1"}, contig_lengths={"chr1": 500})
        self.assertEqual(models["t"].contig, "chr1")
        self.assertEqual(models["t"].exons, ((100, 199), (300, 399)))
        self.assertEqual(models, self.parse(text, contig_aliases={"1": "chr1"}, contig_lengths={"chr1": 500}))

    def test_reference_validation_never_guesses_aliases(self):
        cases = [
            (row(contig="1"), {"contig_lengths": {"chr1": 500}}),
            (row(), {"contig_lengths": {"chr1": 198}}),
            (row(), {"contig_aliases": {"1": "chr1"}}),
            (row(), {"contig_aliases": {"1": "chr2"}, "contig_lengths": {"chr1": 500}}),
            (row(), {"contig_aliases": {"chr1": "chr2"}, "contig_lengths": {"chr1": 500, "chr2": 500}}),
            (row(), {"contig_lengths": {"chr1": True}}),
        ]
        for text, kwargs in cases:
            with self.subTest(kwargs=kwargs), self.assertRaises(IndependentIdentityError):
                self.parse(text, **kwargs)

    def test_conflicting_transcripts_and_exons_fail_with_location(self):
        cases = [
            row() + row(300, 399, contig="chr2"),
            row() + row(300, 399, strand="-"),
            row() + row(300, 399, gene="other"),
            row() + row(),
            row() + row(190, 220),
            row() + row(200, 220),
            row() + row(99, 199, feature="transcript"),
            row(feature="transcript"),
            row() + row(feature="transcript") + row(feature="transcript"),
        ]
        for text in cases:
            with self.subTest(text=text), self.assertRaisesRegex(IndependentIdentityError, r"input.gtf:\d+:"):
                self.parse(text)

    def test_malformed_attributes_and_records_are_rejected(self):
        bad_attributes = [
            'gene_id "g";',
            'transcript_id "t"; transcript_id "t";',
            'gene_id "g"; gene_id "g"; transcript_id "t";',
            'transcript_id "unterminated;',
            'transcript_id "t" gene_id "g";',
            'transcript_id "t\\q";',
            'transcript_id "t\\n";',
            'transcript_id "";',
            'transcript_id;',
            'ID=t;Parent=g;',
        ]
        for attributes in bad_attributes:
            with self.subTest(attributes=attributes), self.assertRaises(IndependentIdentityError):
                self.parse('chr1\tcaller\texon\t1\t9\t.\t+\t.\t' + attributes + '\n')
        for text in [row(start=0), row(start="1.0"), row(strand="."), row().replace('\t.\t+\t.', '\tnan\t+\t.'),
                     row().replace('\t+\t.\t', '\t+\t3\t'), 'a\tb\tc\n', '# only comments\n']:
            with self.subTest(text=text), self.assertRaises(IndependentIdentityError):
                self.parse(text)

    def test_invalid_nonexon_records_are_not_silently_skipped(self):
        with self.assertRaises(IndependentIdentityError):
            self.parse(row() + row(start=0, feature="CDS"))


class IndependentCLITests(unittest.TestCase):
    def run_cli(self, *args, cwd=None):
        return subprocess.run([sys.executable, "-I", str(SCRIPT), *map(str, args)],
                              cwd=cwd, capture_output=True, text=True, check=False)

    def test_isolated_conformance_execution_and_exclusive_report(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "conformance.json"
            first = self.run_cli("conformance", "--vectors", VECTORS, "--output", output, cwd=directory)
            self.assertEqual(first.returncode, 0, first.stderr)
            original = output.read_bytes()
            second = self.run_cli("conformance", "--vectors", VECTORS, "--output", output, cwd=directory)
            self.assertEqual(second.returncode, 2)
            self.assertEqual(output.read_bytes(), original)
            self.assertFalse(list(Path(directory).glob(".independent-identity-*")))

    def test_gtf_report_checksums_and_repeat_stdout_are_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf = root / "models.gtf.gz"
            gtf.write_bytes(gzip.compress((row() + row(300, 399)).encode(), mtime=0))
            lengths = root / "lengths.json"
            lengths.write_text('{"chr1": 500}\n')
            args = ("gtf", "--input", gtf, "--assembly-fingerprint", ASSEMBLY, "--contig-lengths", lengths)
            first = self.run_cli(*args)
            second = self.run_cli(*args)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(first.stdout, second.stdout)
            report = json.loads(first.stdout)
            self.assertEqual(report["inputs"][0]["sha256"], sha256_file(gtf))
            self.assertEqual(report["transcripts"], 1)
            self.assertEqual(report["models"][0]["form"]["public_id"],
                             json.loads(VECTORS.read_text())["vectors"][0]["form"]["public_id"])

    def test_bad_conformance_and_bad_gtf_never_report_success(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            vectors = json.loads(VECTORS.read_text())
            vectors["vectors"][0]["form"]["public_id"] = "txid:TF1.wrong"
            tampered = root / "tampered.json"
            tampered.write_text(json.dumps(vectors))
            result = self.run_cli("conformance", "--vectors", tampered)
            self.assertEqual(result.returncode, 1)
            self.assertFalse(json.loads(result.stdout)["valid"])
            gtf, output = root / "bad.gtf", root / "bad-output.json"
            gtf.write_text(row() + row())
            result = self.run_cli("gtf", "--input", gtf, "--assembly-fingerprint", ASSEMBLY, "--output", output)
            self.assertEqual(result.returncode, 2)
            self.assertFalse(output.exists())
            self.assertEqual(result.stdout, "")


if __name__ == "__main__":
    unittest.main()
