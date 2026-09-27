# Reproducibility checklist

- Identity specification: normative draft 1.0.1 in `docs/superpowers/specs/`;
  the 2026-09-23 clarification records existing fuzzy bridge behavior without
  changing identity algorithms or cluster memberships.
- Registry schema: `schemas/registry-v1.sql`.
- Machine-readable canonical and mapping schemas: `schemas/*.schema.json`.
- Golden vectors: `tests/conformance/exact-v1.json`.
- Test command: `PYTHONPATH=src python -m unittest discover -s tests -v`.
- Synthetic evaluation command and fixed seed: `benchmarks/run_simulation.py`.
- External-input generator: `benchmarks/prepare_external_comparison.py`.
- Checked external runners: `benchmarks/run_gffcompare_comparison.py` and
  `benchmarks/run_isoseql_comparison.py`; TALON commands are recorded in
  `docs/benchmark.md`.
- External evaluator and fixed order seed: `benchmarks/evaluate_external_comparison.py`,
  seed 20260731.
- Checked summary and matrix: `benchmarks/results/`.
- Runtime dependencies: Python standard library only; Python >=3.10.
- License: BSD-3-Clause.
- Conda recipe: `conda-recipe/meta.yaml`, aligned with the current 0.1.3 package.
  `tests/test_packaging.py` checks recipe, project and runtime version consistency.

The checked external comparison uses synthetic truth-labelled inputs. isoSeQL's
SQANTI3-compatible fixtures were generated directly; SQANTI3 was not executed.
Large transient comparator databases and SAM/BAM files are not committed and must
be deposited with the benchmark data archive together with their recorded
checksums.

Before a publication release, record the software and benchmark archive identifiers,
source commit, input checksums, and external comparator versions. Record the
immutable digest of the published final container for an archival run. The
Dockerfile pins an existing TxID/JupyterLab
base image by digest and installs the current TxID wheel on top; that base digest
does not identify the final image. The WDL defaults to the current image tag and
accepts an explicit `docker_image` value, which should use the final published
image digest for archival execution. Historical benchmark results labelled TxID
0.1.0 retain that version and are not relabelled by packaging updates.
