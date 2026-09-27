<p align="center">
  <img src="docs/assets/txid-logo.svg" width="360" alt="TxID — transcript identity">
</p>

# TxID

**Deterministic identities for transcript models.**

[![Tests](https://github.com/dawangran/TxID/actions/workflows/tests.yml/badge.svg)](https://github.com/dawangran/TxID/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-315D80)](pyproject.toml)
[![Research alpha](https://img.shields.io/badge/version-0.1.3%20%7C%20research%20alpha-54766E)](pyproject.toml)
[![BSD 3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-54766E)](LICENSE)

TxID assigns the same exact identifier to the same transcript structure within a
shared reference sequence context. It imports GTF/GFF3 models, preserves matching
reference identifiers, and records structures and their provenance in a SQLite
registry. Mapping tables provide common keys across samples, callers and runs.

TxID is an identity and interoperability layer. Transcript discovery, read
alignment and abundance estimation are performed by upstream tools.

[Quick start](#quick-start) · [Documentation](docs/README.md) ·
[CLI reference](docs/cli.md) · [Specification](docs/spec/identity-v1.md) ·
[Contributing](CONTRIBUTING.md)

## Install

Requires **Python 3.10+**; no third-party runtime dependencies. Version **0.1.3**
is a research alpha. Install from source:

```bash
git clone https://github.com/dawangran/TxID.git
cd TxID
python -m venv .venv
source .venv/bin/activate
python -m pip install .
txid --version
```

For Windows activation, the Conda recipe and container usage, see
[installation](docs/installation.md).

## Quick start

Run the included synthetic example after installation:

```bash
python examples/minimal/run.py
```

It initializes a registry, imports nine models, exports a catalog, validates the
stored identities and draws a gene view. Results go into a fresh temporary
directory whose path is printed. The script checks the expected outputs and a
repeated import; no genome download or external caller is needed.

Open the [example guide](examples/minimal/README.md) to inspect the inputs and
expected results, or follow the [tutorial](docs/tutorial.md) with your own data.

## What you get

| Output | Purpose |
| --- | --- |
| Rewritten GTF | Reference or TxID identifiers with classifications and upstream attributes |
| Mapping TSV | Original identifiers linked to exact forms and splice chains |
| Cohort catalog TSV | Shared structures with observation, sample and tool counts |
| SQLite registry | Canonical objects, full digests, annotation contexts and provenance |
| SVG plots | Registry summaries and gene-level transcript structures |

Use `txid_form` to join exact transcript forms and `txid_sc` to compare complete
multi-exon splice chains. Counts describe model observations, not read counts or
TPM.

## Identity contract

- **SC1** identifies an exact multi-exon splice chain; **TF1** also includes its
  transcript ends. **SE1** identifies an exact single-exon interval.
- Exact IDs include the reference sequence fingerprint, contig and strand.
  Annotation versions affect classification and aliases, not exact identity.
- Exact matching is the default. Optional fuzzy **FC1** groups retain the exact
  IDs; **FC1** groups and **GL1** loci are registry-managed accessions.
- Different reference contexts cannot be silently combined. Identical assembly
  labels alone do not establish the same context.

See [identity and output semantics](docs/identity.md) for coordinate rules,
reference-name handling, gene assignment and the limits of fuzzy grouping. The
[versioned specification](docs/spec/identity-v1.md), [schemas](schemas/README.md)
and [conformance vectors](tests/conformance/exact-v1.json) define the contract.

## Project layout

| Directory | Contents |
| --- | --- |
| [`src/txid/`](src/txid/) | Python package and CLI |
| [`tests/`](tests/) | Unit, integration and conformance tests |
| [`examples/`](examples/minimal/README.md) | Small runnable example with expected results |
| [`docs/`](docs/README.md) | User guides, identity specification and evaluation documentation |
| [`schemas/`](schemas/README.md) | Canonical objects, mapping tables and SQLite schema |
| [`benchmarks/`](benchmarks/README.md) | Evaluation scripts and selected compact reports |
| [`workflows/`](workflows/README.md) | WDL integrations and the Snakemake benchmark workflow |

## Development and validation

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
python -I benchmarks/independent_identity.py conformance --vectors tests/conformance/exact-v1.json
```

CI builds source and wheel distributions, checks the installed package outside the
checkout, and runs the tests and independent conformance check on Python 3.10 and
3.13. See [Contributing](CONTRIBUTING.md) for the development workflow.

The [benchmark guide](docs/benchmark.md) distinguishes synthetic, integration and
fixed-model evidence. Exact identity does not establish biological correctness or
transcript discovery accuracy. Historical results retain their original software
versions; see the [reproducibility checklist](docs/reproducibility.md).

## Citation and license

Record the software version and source commit when reporting results.
[CITATION.cff](CITATION.cff) provides software citation metadata.
TxID is distributed under the [BSD 3-Clause license](LICENSE).

Use [GitHub issues](https://github.com/dawangran/TxID/issues/new/choose) for bug
reports and feature requests.
