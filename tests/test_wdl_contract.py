from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WDL_PATH = ROOT / "workflows" / "txid_batch.wdl"
MULTI_ADD_WDL_PATH = ROOT / "workflows" / "txid_multi_add.wdl"


class WdlContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = WDL_PATH.read_text(encoding="utf-8")

    def test_batch_wdl_uses_platform_compatible_array_types(self):
        self.assertTrue(self.source.startswith("version 1.1\n"))
        self.assertIn("Array[File] gtfs", self.source)
        self.assertNotIn("sample_names", self.source)
        self.assertNotIn("tool_names", self.source)
        self.assertNotRegex(self.source, r"Array\[[^]]+\]\+")

    def test_manifest_metadata_is_derived_without_parallel_arrays(self):
        script = self.source.split("<<'PY'", 1)[1].split("\n    PY", 1)[0]
        script = textwrap.dedent(script).lstrip()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            gtf_list = root / "gtfs.list"
            annotation = root / "annotation.list"
            assembly = root / "assembly.list"
            manifest = root / "manifest.tsv"
            gtf_list.write_text(
                "/localized/S1.isoquant.gtf\n/localized/S2.flair.gtf.gz\n",
                encoding="utf-8",
            )
            annotation.write_text("GENCODE-v49\n", encoding="utf-8")
            assembly.write_text("GRCh38\n", encoding="utf-8")

            result = subprocess.run(
                [
                    sys.executable,
                    "-",
                    str(gtf_list),
                    str(annotation),
                    str(assembly),
                    str(manifest),
                ],
                input=script,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with manifest.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual([row["sample"] for row in rows], ["S1.isoquant", "S2.flair"])
        self.assertEqual([row["tool"] for row in rows], ["unspecified", "unspecified"])

    def test_task_command_is_posix_sh_syntax(self):
        command = self.source.split("command <<<", 1)[1].split(">>>", 1)[0]
        command = textwrap.dedent(command).strip() + "\n"
        result = subprocess.run(
            ["/bin/sh", "-n"],
            input=command,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        for bash_only in (
            "set -eu",
            "mapfile",
            "[[",
            "pipefail",
            "fuzzy_args=(",
            "init_args=(",
        ):
            self.assertNotIn(bash_only, command)

    def test_published_runtime_image_is_the_default(self):
        self.assertIn(
            'String docker_image = "dawang02/txid:0.1.2-jupyter"',
            self.source,
        )


class MultiAddWdlContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = MULTI_ADD_WDL_PATH.read_text(encoding="utf-8")

    def test_interface_uses_direct_multi_add_without_manifest(self):
        self.assertTrue(self.source.startswith("version 1.1\n"))
        self.assertIn("workflow TxIDMultiAdd", self.source)
        self.assertIn("Array[File] gtfs", self.source)
        self.assertIn('"multi-add"', self.source)
        self.assertNotIn("manifest", self.source.lower())
        self.assertNotRegex(self.source, r"Array\[[^]]+\]\+")

    def test_task_command_is_posix_sh_syntax(self):
        command = self.source.split("command <<<", 1)[1].split(">>>", 1)[0]
        command = textwrap.dedent(command).strip() + "\n"
        result = subprocess.run(
            ["/bin/sh", "-n"],
            input=command,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_embedded_runner_passes_all_gtfs_to_multi_add(self):
        script = self.source.split("<<'PY'", 1)[1].split("\n    PY", 1)[0]
        script = textwrap.dedent(script).lstrip()
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bin_dir = root / "bin"
            bin_dir.mkdir()
            log = root / "txid-commands.jsonl"
            fake_txid = bin_dir / "txid"
            fake_txid.write_text(
                "#!/usr/bin/env python3\n"
                "import json, os, sys\n"
                "from pathlib import Path\n"
                "with Path(os.environ['TXID_TEST_LOG']).open('a', encoding='utf-8') as h:\n"
                "    h.write(json.dumps(sys.argv[1:]) + '\\n')\n"
                "if sys.argv[1] == 'init':\n"
                "    Path(sys.argv[sys.argv.index('--db') + 1]).touch()\n"
                "elif sys.argv[1] == 'export':\n"
                "    Path(sys.argv[sys.argv.index('--catalog') + 1]).touch()\n"
                "print('{}')\n",
                encoding="utf-8",
            )
            fake_txid.chmod(0o755)

            fasta = root / "reference.fa"
            reference_gtf = root / "reference.gtf"
            first = root / "sample one.gtf"
            second = root / "sample-two.gtf"
            for path in (fasta, reference_gtf, first, second):
                path.touch()
            gtf_list = root / "gtfs.list"
            assembly = root / "assembly.list"
            annotation = root / "annotation.list"
            tool = root / "tool.list"
            gtf_list.write_text(f"{first}\n{second}\n", encoding="utf-8")
            assembly.write_text("GRCh38\n", encoding="utf-8")
            annotation.write_text("GENCODE-v49\n", encoding="utf-8")
            tool.write_text("IsoQuant\n", encoding="utf-8")

            environment = dict(os.environ)
            environment["PATH"] = f"{bin_dir}:{environment['PATH']}"
            environment["TXID_TEST_LOG"] = str(log)
            result = subprocess.run(
                [
                    sys.executable,
                    "-",
                    str(gtf_list),
                    str(fasta),
                    str(reference_gtf),
                    str(assembly),
                    str(annotation),
                    str(tool),
                    "false",
                    "",
                    "false",
                    "",
                    "false",
                    "",
                ],
                input=script,
                text=True,
                capture_output=True,
                cwd=root,
                env=environment,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            commands = [json.loads(line) for line in log.read_text().splitlines()]

        multi_add = next(command for command in commands if command[0] == "multi-add")
        input_index = multi_add.index("--input")
        tool_index = multi_add.index("--tool")
        self.assertEqual(multi_add[input_index + 1 : tool_index], [str(first), str(second)])
        self.assertEqual(multi_add[tool_index + 1], "IsoQuant")


if __name__ == "__main__":
    unittest.main()
