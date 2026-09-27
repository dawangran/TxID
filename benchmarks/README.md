# Benchmarks and verification

These scripts evaluate identity behavior and interoperability. They do not change
the production identity rules. Start with the [benchmark guide](../docs/benchmark.md)
for evaluated datasets, parameters and limitations.

| Task | Entry point |
| --- | --- |
| Defined synthetic structures and boundary perturbations | `run_simulation.py` |
| Independent exact-ID calculation and conformance | [independent_identity.py](independent_identity.py), [guide](../docs/independent-identity.md) |
| Fixed-model comparisons | `prepare_external_comparison.py`, `run_gffcompare_comparison.py`, `run_isoseql_comparison.py`, `evaluate_external_comparison.py` |
| Bounded real-read IsoQuant–StringTie case | `run_real_interop_callers.py`, `evaluate_real_interop.py` |
| Multi-caller workflow with manifests | [Snakemake workflow](../workflows/publication_benchmark/README.md) |

Run the independent conformance check from the repository root:

```bash
python -I benchmarks/independent_identity.py conformance --vectors tests/conformance/exact-v1.json
```

The real-read runner requires Python 3.11+ and separately installed alignment and
caller tools. Each script documents its arguments through `--help`.

`results/` contains an explicit allowlist of small, recorded reports. Their
versions, checksums and provenance describe the original runs; new code or a
documentation move does not retroactively update those results. Generate new
runs in separate directories, and record their inputs and commands. Large input
data, transient databases and generated outputs belong outside Git.
