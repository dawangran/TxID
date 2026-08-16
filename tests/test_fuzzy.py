from __future__ import annotations

import json
import random
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


def fake_tf_identity(
    public_id: str,
    tss: int,
    tes: int,
    introns: list[list[int]],
    *,
    contig: str = "chr1",
    strand: str = "+",
) -> StructuralIdentity:
    canonical = json.dumps(
        {
            "algorithm": "TF1",
            "assembly": "sha256:" + "0" * 64,
            "contig": contig,
            "introns": introns,
            "strand": strand,
            "tes": tes,
            "tss": tss,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return StructuralIdentity("TF1", public_id, "0", "0" * 64, canonical)


def slow_cluster_forms(identities, *, splice_tolerance, end_tolerance):
    ordered = sorted(identities, key=lambda item: item.public_id)
    by_id = {item.public_id: item for item in ordered}
    clusters: list[list[str]] = []
    statuses: list[str] = []
    for identity in ordered:
        matches = [
            index
            for index, members in enumerate(clusters)
            if all(
                compatible(
                    identity,
                    by_id[member],
                    splice_tolerance=splice_tolerance,
                    end_tolerance=end_tolerance,
                )
                for member in members
            )
        ]
        if len(matches) == 1:
            clusters[matches[0]].append(identity.public_id)
        elif len(matches) > 1:
            clusters.append([identity.public_id])
            statuses.append("ambiguous_bridge")
        else:
            clusters.append([identity.public_id])
            statuses.append("unambiguous")
    return sorted(
        (tuple(sorted(members)), status)
        for members, status in zip(clusters, statuses)
    )


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

    def test_candidate_index_matches_exhaustive_complete_linkage(self):
        rng = random.Random(20260814)
        identities: list[StructuralIdentity] = []
        for index in range(180):
            base = rng.randrange(100, 5000)
            if index % 3 == 0:
                identities.append(
                    fake_identity(
                        f"txid:SE1.{index:04d}",
                        base + rng.randrange(-4, 5),
                        base + 100 + rng.randrange(-4, 5),
                    )
                )
            else:
                identities.append(
                    fake_tf_identity(
                        f"txid:TF1.{index:04d}",
                        base + rng.randrange(-4, 5),
                        base + 400 + rng.randrange(-4, 5),
                        [
                            [
                                base + 100 + rng.randrange(-2, 3),
                                base + 299 + rng.randrange(-2, 3),
                            ]
                        ],
                        contig="chr1" if index % 5 else "chr2",
                        strand="+" if index % 7 else "-",
                    )
                )
        rng.shuffle(identities)

        expected = slow_cluster_forms(
            identities, splice_tolerance=2, end_tolerance=5
        )
        observed = sorted(
            (cluster.members, cluster.bridge_status)
            for cluster in cluster_forms(
                identities, splice_tolerance=2, end_tolerance=5
            )
        )
        self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
