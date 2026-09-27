from __future__ import annotations

import ast
import os
from pathlib import Path
import re
import subprocess
import sys
import unittest

from txid import __version__


ROOT = Path(__file__).resolve().parents[1]


class PackagingVersionTests(unittest.TestCase):
    def setUp(self):
        self.recipe = (ROOT / "conda-recipe" / "meta.yaml").read_text(encoding="utf-8")
        # The metadata uses single-line string literals. Check these directly
        # without adding a TOML/Jinja dependency to tests on Python 3.10.
        project = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
        project_section = re.search(r"(?ms)^\[project\]\n(.*?)(?=^\[|\Z)", project)
        self.assertIsNotNone(project_section, "pyproject.toml must declare [project]")
        version = re.search(r"(?m)^version\s*=\s*(.+)$", project_section.group(1))
        self.assertIsNotNone(version, "the project must declare its release version")
        self.project_version = ast.literal_eval(version.group(1))
        recipe_version = re.search(
            r"(?m)^{%\s*set version\s*=\s*(.+?)\s*%}$", self.recipe
        )
        self.assertIsNotNone(recipe_version, "the recipe must declare its release version")
        self.recipe_version = ast.literal_eval(recipe_version.group(1))

    def test_release_versions_match(self):
        self.assertEqual(self.project_version, __version__)
        self.assertEqual(self.recipe_version, self.project_version)
        package_section = re.search(r"(?ms)^package:\n(.*?)(?=^\S|\Z)", self.recipe)
        self.assertIsNotNone(package_section)
        package_version = re.search(
            r"(?m)^\s+version:\s*(.+)$", package_section.group(1)
        )
        self.assertIsNotNone(package_version)
        rendered = package_version.group(1).replace("{{ version }}", self.recipe_version)
        self.assertEqual(rendered.strip("\"'"), self.project_version)

    def test_conda_version_assertion_accepts_current_source(self):
        assertion = re.search(r'(?m)^\s+- python -c "(.+)"$', self.recipe)
        self.assertIsNotNone(assertion, "the recipe must check the installed TxID version")
        code = assertion.group(1).replace("{{ version }}", self.recipe_version)
        environment = dict(os.environ)
        paths = [str(ROOT / "src"), environment.get("PYTHONPATH", "")]
        environment["PYTHONPATH"] = os.pathsep.join(path for path in paths if path)
        result = subprocess.run(
            [sys.executable, "-B", "-c", code],
            cwd=ROOT,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
