# Minimal example

This example uses a 3 kb artificial reference, six reference transcripts and nine
caller models. The inputs are synthetic fixtures for learning and installation
checks; they do not represent biological samples or discovery performance.

After [installing TxID](../../docs/installation.md), run from the repository root:

```bash
python examples/minimal/run.py
```

The script uses the installed package through `python -m txid`. It creates a fresh
temporary directory and prints its path. No test modules, genome downloads or
external callers are required. To choose a location, supply a **new** directory:

```bash
python examples/minimal/run.py --output-dir /path/to/new/txid-example
```

It runs `init`, `add`, `export`, `validate` and `plot`, then repeats `add` to check
idempotence. The input files stay unchanged; an existing output directory is
rejected.

| File | Contents |
| --- | --- |
| `demo.txid.gtf` | Rewritten transcript models with reference or TxID identifiers |
| `demo.mapping.tsv` | Original identifiers, exact TxIDs and annotation classifications |
| `catalog.tsv` | Nine exact transcript forms |
| `cohort.sqlite` | Registry, structures and import provenance |
| `g1.svg` | Reference match, terminal variant and one-base splice change |
| `summary.json` | Results checked against [expected.json](expected.json) |

The expected classifications are two reference matches, three novel forms assigned
to known genes, two ambiguous gene assignments and two new loci. The end variant
shares SC1 with the reference match but has a different TF1. The one-base splice
change has a different SC1. Reimporting the same file reports `created: false`.

For the individual commands and your own annotations, follow the
[tutorial](../../docs/tutorial.md). See [identity and output semantics](../../docs/identity.md)
before using TxIDs to join measurement tables.
