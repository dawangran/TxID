from __future__ import annotations

import hashlib
import json
import unittest
from pathlib import Path

from txid.identity import identify
from txid.models import Exon, TranscriptModel


class ConformanceTests(unittest.TestCase):
    def test_exact_v1_golden_vectors(self):
        path = Path(__file__).parent / "conformance" / "exact-v1.json"
        document = json.loads(path.read_text(encoding="utf-8"))
        assembly = document["assembly_fingerprint"]
        for vector in document["vectors"]:
            with self.subTest(vector=vector["name"]):
                model = TranscriptModel(
                    vector["contig"],
                    vector["strand"],
                    tuple(Exon(*pair) for pair in vector["exons"]),
                    vector["name"],
                )
                observed = identify(model, assembly)
                expected_sc = vector["splice_chain"]
                self.assertEqual(observed.splice_chain is None, expected_sc is None)
                identities = [(observed.form, vector["form"])]
                if observed.splice_chain is not None:
                    identities.append((observed.splice_chain, expected_sc))
                for identity, expected in identities:
                    self.assertEqual(identity.canonical_json, expected["canonical_json"])
                    independently_hashed = hashlib.sha256(
                        expected["canonical_json"].encode("utf-8")
                    ).hexdigest()
                    self.assertEqual(independently_hashed, expected["full_digest"])
                    self.assertEqual(identity.full_digest, expected["full_digest"])
                    self.assertEqual(identity.public_id, expected["public_id"])


if __name__ == "__main__":
    unittest.main()

