from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "workflows"
    / "publication_benchmark"
    / "scripts"
    / "prepare_isoseql_compat.py"
)
SPEC = importlib.util.spec_from_file_location("prepare_isoseql_compat", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class IsoSeQLCompatibilityTests(unittest.TestCase):
    def test_patch_is_narrow_and_preserves_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            output = root / "compat"
            source.mkdir()
            for name in MODULE.SOURCE_FILES:
                content = "placeholder\n"
                if name == "isoSeQL_db.py":
                    content = f"before\n\t\t\t{MODULE.ORIGINAL} # comment\nafter\n"
                (source / name).write_text(content, encoding="utf-8")
            original = (source / "isoSeQL_db.py").read_bytes()

            audit = MODULE.prepare(source, output)

            patched = (output / "isoSeQL_db.py").read_text(encoding="utf-8")
            self.assertIn(MODULE.REPLACEMENT, patched)
            self.assertNotIn(MODULE.ORIGINAL, patched)
            self.assertEqual((source / "isoSeQL_db.py").read_bytes(), original)
            self.assertFalse(
                audit["compatibility_change"]["structural_queries_changed"]
            )
            self.assertEqual(len(audit["source_files"]), 3)
            self.assertEqual(len(audit["output_files"]), 3)

    def test_refuses_an_unrecognized_upstream_source(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source"
            source.mkdir()
            for name in MODULE.SOURCE_FILES:
                (source / name).write_text("unexpected\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "not found exactly once"):
                MODULE.prepare(source, root / "compat")


if __name__ == "__main__":
    unittest.main()
