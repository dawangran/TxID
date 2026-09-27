"""Focused integrity and boundary tests for the local release archive builder."""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

import build_release_archive as release


class ReleaseArchiveTests(unittest.TestCase):
    def test_reordered_inputs_produce_identical_zip_and_correct_hashes(self):
        first = {"README.md": b"small example\n", "src/example.py": b"print(1)\n"}
        second = dict(reversed(list(first.items())))
        data, manifest = release.render_archive(first)
        self.assertEqual(data, release.render_archive(second)[0])
        with zipfile.ZipFile(io.BytesIO(data)) as archive:
            for info in archive.infolist():
                self.assertEqual(info.date_time, (1980, 1, 1, 0, 0, 0))
                self.assertEqual(info.compress_type, zipfile.ZIP_STORED)
            prefix = release.STEM + "/"
            for line in archive.read(prefix + "SHA256SUMS").decode().splitlines():
                digest, name = line.split("  ", 1)
                self.assertEqual(digest, hashlib.sha256(archive.read(prefix + name)).hexdigest())
            self.assertEqual(json.loads(archive.read(prefix + "RELEASE-MANIFEST.json")), manifest)
            self.assertFalse(manifest["public_deposit_verified"])

    def test_excludes_unapproved_sensitive_and_traversal_paths(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name in ("not-allowed.txt", "../outside.txt", "/absolute.txt", ".ssh/id_rsa", "data/sample.bam"):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    release.collect_files(root, [name], set() if name == "not-allowed.txt" else {name})

    def test_excludes_symlinks_large_files_and_private_keys(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            (root / "real.txt").write_text("fixture")
            (root / "link.txt").symlink_to(root / "real.txt")
            with (root / "large.txt").open("wb") as handle:
                handle.truncate(release.MAX_FILE_BYTES + 1)
            (root / "secret.txt").write_text("-----BEGIN OPENSSH " + "PRIVATE KEY-----\nfixture\n")
            for name in ("link.txt", "large.txt", "secret.txt"):
                with self.subTest(name=name), self.assertRaises(ValueError):
                    release.collect_files(root, [name], {name})

    def test_tiny_test_fixture_is_accepted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            path = root / "tests/fixtures/reference.fa"
            path.parent.mkdir(parents=True)
            path.write_bytes(b">test\nACGT\n")
            name = "tests/fixtures/reference.fa"
            self.assertEqual(release.collect_files(root, [name], {name})[name], b">test\nACGT\n")


if __name__ == "__main__":
    unittest.main(verbosity=2)
