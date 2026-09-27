# TxID documentation

TxID assigns exact structural identifiers to GTF/GFF3 transcript models within a
validated reference context. Start with installation and the tutorial; consult
the specification when implementing or checking identifier compatibility.

## Using TxID

| Guide | Contents |
| --- | --- |
| [Installation](installation.md) | Python environment, Conda recipe and container setup. |
| [Minimal example](../examples/minimal/README.md) | A runnable synthetic example with bundled inputs. |
| [Tutorial](tutorial.md) | Initialize a registry, import models, export a catalog and add an annotation release. |
| [CLI reference](cli.md) | Commands, input formats and output fields. |
| [Identity and reference context](identity.md) | Identifier families, reference matching, provenance and interpretation limits. |
| [Multi-add WDL guide (中文)](txid-multi-add-wdl-guide.zh-CN.md) | Import multiple inputs from one tool and annotation context. |
| [Batch WDL guide (中文)](txid-batch-wdl-guide.zh-CN.md) | Import a cohort with per-input provenance. |

## Specification and implementation

- [Identity v1 specification](spec/identity-v1.md): normative draft 1.0.1 for
  SC1, TF1 and SE1; also defines registry and fuzzy-grouping requirements.
- [Machine-readable schemas](../schemas/): canonical objects, mappings and the
  SQLite registry schema.
- [Conformance vectors](../tests/conformance/exact-v1.json): expected canonical
  bytes, digests and exact identifiers.
- [Independent identity calculator](independent-identity.md): a separately coded
  structural verifier, its supported input domain and its limitations.

The specification version, identity algorithm families and software release
number are distinct. A software update does not change an existing identifier
family's meaning.

## Evaluation and reproducibility

- [Benchmark guide](benchmark.md): synthetic evaluation and external comparison
  commands, metrics and interpretation.
- [Benchmark protocol](design/benchmark-protocol.md): separate discovery and
  identity tracks, dataset tiers and required provenance.
- [Benchmark workflow](../workflows/publication_benchmark/README.md): Snakemake
  configuration and execution.
- [Reproducibility checklist](reproducibility.md): source versions, input
  checksums, environment records and archival requirements.

Historical execution reports retain the software versions and evidence limits of
their original runs:

| Report | Scope |
| --- | --- |
| [Tier-0 integration](benchmark/publication-smoke-2026-07-31.md) | Actual caller and comparator executables on a small synthetic read fixture. |
| [ENCODE full-genome audit](benchmark/encode-gtf-full-2026-07-31.md) | Fixed-model identity comparison using six released WTC11 TALON annotations. |
| [ENCODE pilot](benchmark/encode-gtf-pilot-2026-07-31.md) | Superseded chromosome-restricted preflight record. |

Large inputs, transient databases and generated execution bundles are kept
outside the source repository. Compact checked reports remain under
[`benchmarks/results/`](../benchmarks/results/).
