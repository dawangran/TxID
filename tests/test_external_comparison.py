from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "benchmarks" / "evaluate_external_comparison.py"
SPEC = importlib.util.spec_from_file_location("txid_external_evaluator", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
EVALUATOR = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(EVALUATOR)


class ExternalComparisonEvaluatorTests(unittest.TestCase):
    def test_pairwise_rates_count_one_split_and_one_merge(self):
        result = EVALUATOR._pairwise_error_rates(
            ["truth-a", "truth-a", "truth-b"],
            ["group-1", "group-2", "group-2"],
        )
        self.assertEqual(result["false_split_pairs"], 1)
        self.assertEqual(result["false_merge_pairs"], 1)
        self.assertEqual(result["false_split_rate"], 1.0)
        self.assertEqual(result["false_merge_rate"], 1.0)

    def test_gffcompare_tracking_maps_query_observations(self):
        with tempfile.TemporaryDirectory() as temporary:
            tracking = Path(temporary) / "combined.tracking"
            tracking.write_text(
                "TCONS_1|2|200\tXLOC_1\t-\tu\t"
                "q1:G1|caller.sample.1|2|0|0|0|200\t"
                "q2:G1|caller.sample.2|2|0|0|0|200\n",
                encoding="utf-8",
            )
            self.assertEqual(
                EVALUATOR._gffcompare_mapping(tracking),
                {"caller.sample.1": "TCONS_1", "caller.sample.2": "TCONS_1"},
            )

    def test_order_check_ignores_renamed_but_equal_partitions(self):
        first = {"a": "1", "b": "1", "c": "2"}
        renamed = {"a": "9", "b": "9", "c": "8"}
        result = EVALUATOR._invariance(first, renamed)
        self.assertTrue(result["partition_invariant"])
        self.assertFalse(result["literal_label_invariant"])
        self.assertEqual(result["partition_changed_observations"], 0)
        self.assertEqual(result["literal_label_changed_observations"], 3)


if __name__ == "__main__":
    unittest.main()
