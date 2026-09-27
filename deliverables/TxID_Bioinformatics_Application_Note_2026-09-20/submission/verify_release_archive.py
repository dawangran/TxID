#!/usr/bin/env python3
"""Verify the local ZIP using an isolated venv, synthetic inputs and no downloads."""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile
import zipfile

import build_release_archive as release


def verify(submission: Path) -> dict:
    archive = submission / (release.STEM + ".zip")
    original = archive.read_bytes()
    report = {"archive": archive.name, "archive_sha256": hashlib.sha256(original).hexdigest(),
              "python": sys.version.split()[0], "public_deposit_verified": False, "checks": []}

    def run(label: str, command: list[str], cwd: Path, env: dict | None = None) -> str:
        result = subprocess.run(command, cwd=cwd, env=env, text=True, capture_output=True)
        output = result.stdout + result.stderr
        report["checks"].append({"name": label, "command": command, "returncode": result.returncode,
                                 "output_tail": output[-5000:]})
        print(f"{label}: {'PASS' if result.returncode == 0 else 'FAIL'}", flush=True)
        if result.returncode:
            raise RuntimeError(label + "\n" + output[-5000:])
        return output

    try:
        with tempfile.TemporaryDirectory(prefix="txid-release-validation-") as temporary:
            temp = Path(temporary)
            with zipfile.ZipFile(archive) as handle:
                assert handle.testzip() is None
                names = handle.namelist()
                assert len(names) == len(set(names))
                assert all(name.startswith(release.STEM + "/") and ".." not in Path(name).parts for name in names)
                for line in handle.read(release.STEM + "/SHA256SUMS").decode().splitlines():
                    digest, name = line.split("  ", 1)
                    assert hashlib.sha256(handle.read(release.STEM + "/" + name)).hexdigest() == digest
                handle.extractall(temp)
            root = temp / release.STEM
            report["checks"].append({"name": "archive payload hashes and CRC", "returncode": 0,
                                     "members": len(names)})
            run("rebuild extracted snapshot", [sys.executable, "-B", str(root / release.PACKAGE / "build_release_archive.py"),
                "--repo", str(root), "--output-dir", str(temp / "rebuilt")], root)
            assert (temp / "rebuilt" / archive.name).read_bytes() == original
            report["checks"].append({"name": "byte-identical rebuild from extracted snapshot", "returncode": 0})
            run("wheel build offline", [sys.executable, "-B", "-m", "pip", "wheel", "--no-index", "--no-cache-dir",
                "--no-deps", "--no-build-isolation", "--wheel-dir", str(temp / "wheels"), str(root)], temp)
            wheel = next((temp / "wheels").glob("*.whl"))
            environment = temp / "venv"
            run("create isolated venv", [sys.executable, "-B", "-m", "venv", str(environment)], temp)
            interpreter = str(environment / "bin/python")
            env = dict(os.environ)
            env.pop("PYTHONPATH", None)
            env["PYTHONDONTWRITEBYTECODE"] = "1"
            env["PATH"] = str(environment / "bin") + os.pathsep + env.get("PATH", "")
            run("wheel install offline", [interpreter, "-B", "-m", "pip", "install", "--no-index", "--no-cache-dir",
                "--no-deps", str(wheel)], temp, env)
            run("installed package version", [interpreter, "-B", "-c",
                'import txid; assert txid.__version__ == "0.1.3"; print(txid.__version__, txid.__file__)'], temp, env)
            output = run("included software and benchmark tests", [interpreter, "-B", "-m", "unittest",
                "discover", "-s", "tests", "-v"], root, env)
            match = re.search(r"Ran (\d+) tests?", output)
            report["included_test_count"] = int(match.group(1)) if match else None
            run("standalone independent golden conformance", [interpreter, "-I", "-B",
                str(root / "benchmarks/independent_identity.py"), "conformance", "--vectors",
                str(root / "tests/conformance/exact-v1.json"), "--output", str(temp / "conformance.json")], temp, env)
            smoke = temp / "smoke"
            smoke.mkdir()
            run("create embedded synthetic tutorial inputs", [interpreter, "-B", "-c",
                "from pathlib import Path; from tests.helpers import write_inputs; write_inputs(Path(" + repr(str(smoke)) + "))"], root, env)
            cli = [interpreter, "-B", "-m", "txid"]
            run("tutorial init", cli + ["init", "--db", str(smoke / "cohort.sqlite"), "--fasta",
                str(smoke / "reference.fa"), "--assembly", "synthetic-v1", "--annotation",
                str(smoke / "reference.gtf"), "--annotation-name", "ref-v1"], smoke, env)
            for sample in ("A", "B"):
                run("tutorial add " + sample, cli + ["add", "--db", str(smoke / "cohort.sqlite"), "--input",
                    str(smoke / "caller.gtf"), "--sample", sample, "--tool", "synthetic-fixture", "--annotation-name",
                    "ref-v1", "--output-gtf", str(smoke / (sample + ".gtf")), "--mapping", str(smoke / (sample + ".tsv"))], smoke, env)
            run("tutorial export", cli + ["export", "--db", str(smoke / "cohort.sqlite"), "--catalog",
                str(smoke / "catalog.tsv")], smoke, env)
            run("tutorial validate", cli + ["validate", "--db", str(smoke / "cohort.sqlite")], smoke, env)
            rows = []
            for sample in ("A", "B"):
                with (smoke / (sample + ".tsv")).open() as handle:
                    rows.append(list(csv.DictReader(handle, delimiter="\t")))
            assert len(rows[0]) == len(rows[1]) == 9
            assert sorted(row["txid_form"] for row in rows[0]) == sorted(row["txid_form"] for row in rows[1])
            report["checks"].append({"name": "two sample mappings share all 9 exact form keys", "returncode": 0})
            report["status"] = "passed"
    except Exception as error:
        report["status"] = "failed"
        report["error"] = str(error)
    (submission / "release-validation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    result = verify(Path(__file__).resolve().parent)
    print(json.dumps({key: value for key, value in result.items() if key != "checks"}, indent=2))
    sys.exit(0 if result["status"] == "passed" else 1)
