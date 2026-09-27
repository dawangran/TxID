# Contributing to TxID

Bug reports and focused pull requests are welcome. Include a small reproducible
example and describe the behavior you expect.

## Development setup

From a clone of this repository, using Python 3.10 or later:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
txid --version
```

On Windows, activate the environment with `.venv\Scripts\Activate.ps1` in
PowerShell. The core package and its unittest suite use the standard library;
external benchmark callers require separate environments.

## Where changes belong

- `src/txid/` contains the CLI, shared transcript model, identity logic and
  registry implementation. All supported input formats use the same identity
  rules.
- `tests/` contains unit and integration tests, small synthetic fixtures and
  exact-identity conformance vectors.
- `schemas/` defines canonical objects, mapping records and the SQLite schema.
- `docs/` contains user documentation and the versioned specification.
- `benchmarks/` and `workflows/` contain reproducible evaluations and workflow
  integrations; generated inputs and run outputs belong outside the source tree.

## Identity compatibility

For an identity-semantics change, first update the
[identity specification](docs/spec/identity-v1.md). Keep the implementation,
schemas, documentation and conformance vectors consistent. Changes to canonical
inputs, serialization, coordinates, digest rules or identifier meaning require a
new identity algorithm family; a package version bump alone is insufficient.

Preserve exact identifiers across input order, sample names and repeated imports.
Keep fuzzy similarity separate from exact identity. Registry changes must retain
provenance, use explicit schema versions and leave no partial import on failure.

## Verification

Run these commands from the repository root:

```bash
python -m unittest discover -s tests -v
python -I benchmarks/independent_identity.py conformance \
  --vectors tests/conformance/exact-v1.json
```

Add tests for changed behavior using small synthetic cases with explicit expected
results. Check both strands when changing coordinate handling. For export changes,
check deterministic bytes as well as parsed values. Documentation-only changes
need working links and accurate commands.

CI also builds the source distribution and wheel, installs a wheel built from the
source distribution in a fresh environment, and runs the minimal example outside
the checkout to check the installed package.

The source distribution contains the installable package, schemas and minimal
example. Use a Git checkout for the full development tests and benchmark tools.
The build follows the [PyPA source-to-wheel workflow](https://build.pypa.io/en/stable/how-to/basic-usage.html).

## Pull requests and repository contents

Keep each pull request focused. Explain the problem, the resulting behavior and
the checks you ran. Update the CLI documentation when user-visible behavior
changes. Preserve version labels and provenance in historical benchmark reports.

Use synthetic examples in issues and tests. Do not commit manuscripts, submission
materials, private sample identifiers, credentials, large reference genomes,
reads, alignments, working databases or generated benchmark runs. Record public
data accessions and checksums in manifests instead. Existing compact benchmark
summaries are retained as reproducibility evidence. Review `git status` and the
staged diff before committing; ignore rules do not protect files already tracked.
