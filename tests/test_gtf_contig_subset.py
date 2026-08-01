from __future__ import annotations

import gzip
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "workflows"
    / "publication_benchmark"
    / "scripts"
    / "stage_gtf_contig_subset.py"
)
SPEC = importlib.util.spec_from_file_location("stage_gtf_contig_subset", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
SUBSET = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SUBSET
SPEC.loader.exec_module(SUBSET)


class GtfContigSubsetTests(unittest.TestCase):
    def test_filter_is_byte_deterministic_and_retains_comments(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "models.gtf.gz"
            with gzip.open(source, "wt", encoding="utf-8") as handle:
                handle.write(
                    "# source\n"
                    'chr22\tt\ttranscript\t1\t9\t.\t+\t.\tgene_id "g"; '
                    'transcript_id "keep";\n'
                    'chr22\tt\texon\t1\t9\t.\t+\t.\tgene_id "g"; '
                    'transcript_id "keep";\n'
                    'chr1\tt\texon\t1\t9\t.\t+\t.\tgene_id "x"; '
                    'transcript_id "drop";\n'
                )
            first = root / "first.gtf.gz"
            second = root / "second.gtf.gz"
            record = SUBSET.filter_gtf(source, first, {"chr22"})
            SUBSET.filter_gtf(source, second, {"chr22"})
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(record["records"], 2)
            self.assertEqual(record["features"], {"exon": 1, "transcript": 1})
            with gzip.open(first, "rt", encoding="utf-8") as handle:
                output = handle.read()
            self.assertTrue(output.startswith("# source\n"))
            self.assertIn('transcript_id "keep"', output)
            self.assertNotIn('transcript_id "drop"', output)


if __name__ == "__main__":
    unittest.main()
