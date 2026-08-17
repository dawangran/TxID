from __future__ import annotations

import contextlib
import io
import json
import sqlite3
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from txid.cli import main
from txid.workflow import import_file

from tests.helpers import make_registry


class GenePlotTests(unittest.TestCase):
    def invoke(self, arguments: list[str]) -> tuple[int, dict[str, object], str]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = main(arguments)
        result = json.loads(stdout.getvalue()) if stdout.getvalue() else {}
        return status, result, stderr.getvalue()

    def import_fuzzy_fixture(self, root: Path) -> Path:
        database, caller = make_registry(root)
        import_file(
            database,
            input_path=caller,
            sample="sample-A",
            tool="IsoQuant",
            annotation_name="ref-v1",
            output_gtf=root / "rewritten.gtf",
            mapping_path=root / "mapping.tsv",
            fuzzy_splice_tolerance=5,
            fuzzy_end_tolerance=15,
        )
        return database

    def test_gene_plot_shows_exact_new_upstream_old_and_fuzzy_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = self.import_fuzzy_fixture(root)
            first = root / "g1.first.svg"
            second = root / "g1.second.svg"

            status, result, stderr = self.invoke(
                ["plot", "--db", str(database), "--gene", "g1", "--output", str(first)]
            )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(result["gene_ids"], ["g1"])
            self.assertEqual(result["reference_tracks"], 1)
            self.assertEqual(result["exact_forms"], 3)
            self.assertEqual(result["observations"], 3)

            document = ET.parse(first)
            namespace = {"svg": "http://www.w3.org/2000/svg"}
            observed = document.findall(".//svg:g[@class='observed-track']", namespace)
            self.assertEqual(len(observed), 3)
            form_ids = {track.attrib["data-form-id"] for track in observed}
            fuzzy_ids = {track.attrib["data-fuzzy-cluster"] for track in observed}
            self.assertTrue(all(value.startswith("txid:TF1.") for value in form_ids))
            self.assertEqual(len(fuzzy_ids), 1)
            self.assertTrue(next(iter(fuzzy_ids)).startswith("txid:FC1."))

            text = first.read_text(encoding="utf-8")
            for old_id in ("alt_known", "end_variant", "one_bp_splice"):
                self.assertIn(old_id, text)
            for form_id in form_ids:
                self.assertIn(form_id, text)
            self.assertIn("Splice-chain TxID: txid:SC1.", text)
            self.assertIn("tx_ref1", text)
            self.assertIn("splice ±5 bp", text)
            self.assertIn("ends ±15 bp", text)
            self.assertIn("Upstream (old):", text)

            status, _, stderr = self.invoke(
                ["plot", "--db", str(database), "--gene", "g1", "--output", str(second)]
            )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_gene_name_and_upstream_gene_id_are_supported(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database = self.import_fuzzy_fixture(root)

            status, symbol, stderr = self.invoke(
                [
                    "plot",
                    "--db",
                    str(database),
                    "--gene",
                    "GENE1",
                    "--output",
                    str(root / "symbol.svg"),
                ]
            )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(symbol["resolved_by"], "gene_name")
            self.assertEqual(symbol["gene_ids"], ["g1"])

            status, upstream, stderr = self.invoke(
                [
                    "plot",
                    "--db",
                    str(database),
                    "--gene",
                    "up_g1",
                    "--output",
                    str(root / "upstream.svg"),
                ]
            )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(upstream["resolved_by"], "original_gene_id")
            self.assertEqual(upstream["exact_forms"], 3)

    def test_unknown_gene_is_an_error_and_writes_no_partial_svg(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            output = root / "missing.svg"

            status, result, stderr = self.invoke(
                [
                    "plot",
                    "--db",
                    str(database),
                    "--gene",
                    "missing-gene",
                    "--output",
                    str(output),
                ]
            )
            self.assertEqual(status, 2)
            self.assertEqual(result, {})
            self.assertIn("gene not found in registry", stderr)
            self.assertFalse(output.exists())

    def test_reference_only_gene_can_be_plotted(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            output = root / "reference-only.svg"

            status, result, stderr = self.invoke(
                ["plot", "--db", str(database), "--gene", "gSE", "--output", str(output)]
            )
            self.assertEqual(status, 0, stderr)
            self.assertEqual(result["reference_tracks"], 1)
            self.assertEqual(result["exact_forms"], 0)
            ET.parse(output)

    def test_ambiguous_gene_name_is_not_silently_combined(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            database, _ = make_registry(root)
            with sqlite3.connect(database) as connection:
                connection.execute(
                    """
                    UPDATE reference_transcript
                    SET attributes_json = '[["gene_name","DUPLICATE"]]'
                    WHERE reference_gene_id IN ('gA', 'gB')
                    """
                )
            output = root / "ambiguous.svg"

            status, _, stderr = self.invoke(
                [
                    "plot",
                    "--db",
                    str(database),
                    "--gene",
                    "DUPLICATE",
                    "--output",
                    str(output),
                ]
            )
            self.assertEqual(status, 2)
            self.assertIn("matching gene IDs: gA, gB", stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
