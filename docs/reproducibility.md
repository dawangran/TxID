# Reproducibility checklist

- Identity specification: normative draft 1.0.0 in `docs/superpowers/specs/`.
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

The checked external comparison uses synthetic truth-labelled inputs. isoSeQL's
SQANTI3-compatible fixtures were generated directly; SQANTI3 was not executed.
Large transient comparator databases and SAM/BAM files are not committed and must
be deposited with the manuscript data archive together with their recorded
checksums.

Before a publication release, replace manuscript placeholders with author,
repository, RRID/bio.tools, archive DOI, GigaDB/Zenodo DOI, funding, and external
comparator version information. Build and publish a digest-pinned container; the
development Dockerfile currently pins the Python version but not an immutable
base-image digest.
