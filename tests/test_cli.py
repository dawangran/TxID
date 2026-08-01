from __future__ import annotations

import contextlib
import io
import json
import tempfile
import unittest
from pathlib import Path

from txid.cli import main

from tests.helpers import write_inputs


class CliTests(unittest.TestCase):
    def invoke(self, arguments: list[str]) -> tuple[int, dict]:
        stdout = io.StringIO()
        stderr = io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            status = main(arguments)
        self.assertEqual(status, 0, stderr.getvalue())
        return status, json.loads(stdout.getvalue())

    def test_command_surface(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fasta, reference, caller = write_inputs(root)
            database = root / "cli.sqlite"
            _, initialized = self.invoke(
                [
                    "init", "--db", str(database), "--fasta", str(fasta),
                    "--assembly", "synthetic-v1", "--annotation", str(reference),
                    "--annotation-name", "ref-v1",
                ]
            )
            self.assertEqual(initialized["reference_transcripts"], 6)
            _, added = self.invoke(
                [
                    "add", "--db", str(database), "--input", str(caller),
                    "--sample", "s", "--tool", "tool", "--annotation-name", "ref-v1",
                    "--output-gtf", str(root / "out.gtf"), "--mapping", str(root / "map.tsv"),
                ]
            )
            self.assertEqual(added["transcripts"], 9)
            _, exported = self.invoke(
                ["export", "--db", str(database), "--catalog", str(root / "catalog.tsv")]
            )
            self.assertGreater(exported["rows"], 0)
            _, validation = self.invoke(["validate", "--db", str(database)])
            self.assertTrue(validation["valid"])
            first_form = (root / "map.tsv").read_text(encoding="utf-8").splitlines()[1].split("\t")[7]
            _, inspected = self.invoke(["inspect", "--db", str(database), first_form])
            self.assertEqual(inspected["query"], first_form)
            _, plotted = self.invoke(
                ["plot", "--db", str(database), "--output", str(root / "overview.svg")]
            )
            self.assertTrue(Path(plotted["plot"]).is_file())


if __name__ == "__main__":
    unittest.main()
