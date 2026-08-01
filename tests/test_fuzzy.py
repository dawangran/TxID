from __future__ import annotations

import json
import unittest

from txid.fuzzy import cluster_forms, compatible
from txid.models import StructuralIdentity


def fake_identity(public_id: str, start: int, end: int) -> StructuralIdentity:
    canonical = json.dumps(
        {
            "algorithm": "SE1",
            "assembly": "sha256:" + "0" * 64,
            "contig": "chr1",
            "end": end,
            "start": start,
            "strand": "+",
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return StructuralIdentity("SE1", public_id, "0", "0" * 64, canonical)


class FuzzyTests(unittest.TestCase):
    def test_single_exon_tolerance(self):
        a = fake_identity("txid:SE1.a", 100, 200)
        b = fake_identity("txid:SE1.b", 103, 204)
        self.assertTrue(compatible(a, b, splice_tolerance=0, end_tolerance=4))
        self.assertFalse(compatible(a, b, splice_tolerance=0, end_tolerance=3))

    def test_complete_linkage_exposes_ambiguous_bridge(self):
        left = fake_identity("txid:SE1.a", 100, 200)
        right = fake_identity("txid:SE1.c", 110, 210)
        bridge = fake_identity("txid:SE1.z", 105, 205)
        clusters = cluster_forms(
            [bridge, right, left], splice_tolerance=0, end_tolerance=5
        )
        self.assertEqual([cluster.members for cluster in clusters], [(left.public_id,), (right.public_id,), (bridge.public_id,)])
        self.assertEqual(clusters[-1].bridge_status, "ambiguous_bridge")


if __name__ == "__main__":
    unittest.main()

