from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from txid.simulation import run_simulation


class SimulationTests(unittest.TestCase):
    def test_small_truth_defined_run(self):
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "results"
            summary = run_simulation(
                output,
                seed=17,
                transcript_count=20,
                sample_count=2,
                tool_count=2,
                artifact_rate=0.1,
            )
            exact = summary["exact_structural_identity"]
            self.assertEqual(exact["false_split_rate"], 0.0)
            self.assertEqual(exact["false_merge_rate"], 0.0)
            self.assertTrue(summary["annotation_release_check"]["exact_txid_invariant"])
            self.assertGreater(
                summary["sample_by_transcript_matrix"]["raw_upstream_columns"],
                summary["sample_by_transcript_matrix"]["txid_exact_columns"],
            )
            self.assertTrue((output / "simulation-summary.json").is_file())
            self.assertTrue((output / "sample-by-txid.tsv").is_file())


if __name__ == "__main__":
    unittest.main()
