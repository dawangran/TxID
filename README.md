<p align="center">
  <img src="docs/assets/txid-logo.svg" width="400" alt="TxID — transcript identity">
</p>

# TxID

**Reference-aware, reproducible identities for long-read transcript models.**

[![Tests](https://github.com/dawangran/TxID/actions/workflows/tests.yml/badge.svg)](https://github.com/dawangran/TxID/actions/workflows/tests.yml)
[![Python 3.10+](https://img.shields.io/badge/Python-3.10%2B-315D80)](pyproject.toml)
[![Research alpha](https://img.shields.io/badge/version-0.1.3%20%7C%20research%20alpha-54766E)](pyproject.toml)
[![BSD 3-Clause](https://img.shields.io/badge/license-BSD--3--Clause-54766E)](LICENSE)

TxID gives independently submitted GTF/GFF3 transcript models a common structural
key. Its versioned specification defines exact splice-chain, transcript-form and
single-exon identities within one reference sequence collection. Equivalent
structures receive the same exact identifiers across samples, tools and import
orders, including across independently initialized registries using that context.

A SQLite registry keeps those structures separate from observations, reference
annotation classifications and provenance. Existing reference identifiers are
preserved where a model exactly matches the reference; TxIDs provide a common key
for comparing models and joining externally produced measurement tables. Transcript
discovery, read alignment and abundance estimation remain upstream steps.

[Quick start](#quick-start) · [Identity rules](#identity-rules) ·
[CLI reference](docs/cli.md) · [Tutorial](docs/tutorial.md) ·
[Specification](docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md)

## Installation

TxID **0.1.3 is a research alpha**. The package requires Python **3.10 or later**
and uses only the Python standard library at runtime. Install from this repository:

```bash
git clone https://github.com/dawangran/TxID.git
cd TxID
python -m venv .venv
source .venv/bin/activate
python -m pip install .
txid --version
```

The version command prints `txid 0.1.3`. On Windows, activate the environment with
`.venv\Scripts\Activate.ps1` in PowerShell. A local
[Conda recipe](conda-recipe/meta.yaml) is included for building the package.

For reproducible analyses, record the source commit, software version, reference
FASTA, annotation release and command options. The registry records input checksums
and import provenance.

## Quick start

Run this example from the repository root after installation, using Bash or another
POSIX shell. It writes the small synthetic fixtures in
[`tests/helpers.py`](tests/helpers.py) into a fresh temporary directory. No genome
download or caller installation is needed; the nine input models are test cases,
not biological samples.

```bash
TXID_DEMO="$(mktemp -d "${TMPDIR:-/tmp}/txid-demo.XXXXXX")"
python - "$TXID_DEMO" <<'PY'
from pathlib import Path
import sys
from tests.helpers import write_inputs

write_inputs(Path(sys.argv[1]))
PY

txid init \
  --db "$TXID_DEMO/cohort.sqlite" \
  --fasta "$TXID_DEMO/reference.fa" \
  --assembly synthetic-v1 \
  --annotation "$TXID_DEMO/reference.gtf" \
  --annotation-name ref-v1

txid add \
  --db "$TXID_DEMO/cohort.sqlite" \
  --input "$TXID_DEMO/caller.gtf" \
  --sample demo \
  --tool synthetic-fixture \
  --annotation-name ref-v1 \
  --output-gtf "$TXID_DEMO/demo.txid.gtf" \
  --mapping "$TXID_DEMO/demo.mapping.tsv"

txid export --db "$TXID_DEMO/cohort.sqlite" \
  --catalog "$TXID_DEMO/cohort.catalog.tsv"
txid validate --db "$TXID_DEMO/cohort.sqlite"
txid inspect --db "$TXID_DEMO/cohort.sqlite" tx_ref1
txid plot --db "$TXID_DEMO/cohort.sqlite" \
  --gene g1 --output "$TXID_DEMO/g1.svg"
```

The import reports **9 transcripts**, the catalog contains **9 exact forms**, and
validation returns `"valid": true`. The gene plot shows a reference match, an end
variant and a one-base splice change. Inspect the files in `$TXID_DEMO`; the input
annotations remain unchanged. Repeating the same `add` command reuses the existing
import and reports `"created": false`.

## Outputs

| Output | Contents and use |
| --- | --- |
| Rewritten GTF | Reference or TxID transcript names, classifications, exact IDs and retained upstream attributes. Reserved attributes are namespaced where necessary. |
| Per-input mapping TSV | Original sample/tool/gene/transcript identifiers, `txid_form`, `txid_sc`, output names, annotation context, gene candidates and optional fuzzy assignments. |
| Cohort catalog TSV | Exact forms and splice chains with gene assignments, classifications and observation/sample/tool counts. |
| SQLite registry | Canonical objects, full digests, observations, annotation contexts, aliases, locus assignments and import manifests. |
| Optional SVG | A registry overview or a gene/locus view comparing reference and observed exon structures. |

Use **`txid_form`** as the exact-form join key in mapping tables. It contains TF1
for multi-exon models and SE1 for single-exon models. Use `txid_sc` when the intended
comparison is the complete multi-exon splice chain; single-exon models have no SC1.
Mapping and catalog counts describe model observations, not read counts or TPM.

## Identity rules

| Family | Object represented | Assignment and portability |
| --- | --- | --- |
| `txid:SC1.<digest>` | Exact complete intron chain of a multi-exon model; outer transcript ends are excluded. | Computed from structure and reference context; reproducible across registries sharing that context. |
| `txid:TF1.<digest>` | Exact multi-exon form: intron chain plus model-defined transcription start (TSS) and end (TES). | Computed exact key; changing an end changes TF1 even when SC1 is unchanged. |
| `txid:SE1.<digest>` | Exact single-exon interval and strand. | Computed exact key; single-exon models are classified separately. |
| `txid:FC1.<accession>` | Optional fuzzy cluster with recorded splice/end tolerances and membership. | Registry-managed; not a portable exact structural key. |
| `txid:GL1.<accession>` | New or unresolved gene/locus assignment. | Registry-managed; not a portable exact structural key. |

Exact objects include the algorithm family, reference fingerprint, primary contig
name and strand. Coordinates are **1-based closed**. Canonical JSON fixes key
ordering and serialization; SHA-256 supplies the full digest and its first 24
hexadecimal characters form the public suffix. Registries retain the complete
canonical object and digest. A conflicting public identifier is an error, never
repaired by adding a suffix.

Sample labels, tool names, annotation releases and import order are excluded from
exact hash inputs. Changing identity semantics requires a new algorithm-family
version. The [normative specification](docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md),
[schemas](schemas/) and [golden conformance vectors](tests/conformance/exact-v1.json)
define the contract independently of the CLI and registry.

### Reference and annotation context

The reference fingerprint includes normalized sequence digests, contig lengths
and **primary FASTA names**. Matching assembly labels alone are insufficient:
renaming primary contigs or changing the sequence collection changes the context.
Input contigs must match primary names or use an explicit alias table whose targets
exist in the selected FASTA. TxID does not infer aliases by stripping `chr`.
Registries with different reference fingerprints cannot be combined as exact
identity contexts; cross-assembly liftover is not implemented.

Annotation versions have separate normalized fingerprints. On the same reference
context, changing annotation release can change classification and reference
aliases while leaving exact structural IDs unchanged:

| Relationship to the selected annotation | Rewritten identifiers |
| --- | --- |
| Exact reference-form match | Preserve the selected reference `gene_id` and `transcript_id`; record TxIDs separately. |
| No exact match; overlap with one same-strand reference gene span | Retain that reference `gene_id`; use the TF1 or SE1 transcript ID. |
| Multiple candidate genes or no same-strand gene overlap | Keep ambiguous candidates explicit where present; assign a GL1 locus and a TF1 or SE1 transcript ID. |

Gene-span overlap is an assignment rule, not evidence of biological gene membership.
Novelty is relative to the annotation and is never part of a permanent exact ID.

## Import samples and annotation releases

TxID reads GTF and GFF3, including gzip-compressed inputs and exon models lacking
transcript meta-features. Inputs from callers such as IsoQuant, TALON, FLAIR and
StringTie enter the same model representation. Unknown contigs, inconsistent
models, overlapping exons and out-of-range coordinates are rejected. GFF3 inputs
are rewritten as GTF; this is not a byte-preserving format conversion.

For multiple files from one tool and annotation context:

```bash
txid multi-add --db cohort.sqlite \
  --input calls/sample-a.gtf calls/sample-b.gtf \
  --samples sample-a sample-b \
  --tool IsoQuant --annotation-name release-1 --output-dir results
```

For heterogeneous inputs, use `txid batch` with a tab-delimited manifest containing
`input`, `sample`, `tool` and `annotation_name`; `format` is optional. Register
another annotation on the same reference with:

```bash
txid annotation-add --db cohort.sqlite \
  --annotation reference-v2.gtf --annotation-name release-2
```

Each input is staged and validated before its atomic registry transaction. An
identical sample/tool/input retry is idempotent. A batch consists of separate
per-input transactions, so a later failed input does not undo earlier completed
imports. Exports from a fixed registry are sorted deterministically; exact-ID
invariance does not imply identical locus assignments or complete catalog bytes
across different import histories.

Exact mode is the default and needs no tolerance. Optional fuzzy grouping requires
both `--fuzzy-splice-tolerance` and `--fuzzy-end-tolerance` in base pairs. Every pair
in a cluster must meet the tolerances, and ambiguous bridges are marked. Exact IDs
remain available alongside FC1 assignments. See the [CLI reference](docs/cli.md)
for full command options and the [tutorial](docs/tutorial.md) for context changes.

## Containers and workflow integration

The repository includes a [Dockerfile](Dockerfile) for TxID and JupyterLab. Build
and inspect a local image with:

```bash
docker build -t txid:0.1.3-jupyter .
docker run --rm txid:0.1.3-jupyter txid --version
docker run --rm -p 8888:8888 -v "$PWD:/workspace" txid:0.1.3-jupyter
```

The final command starts JupyterLab and prints its access token. The container's
working user has UID 1000; the mounted directory must be writable by that user.
Record the **final image digest** for archival runs; the pinned base-image digest
does not identify the complete built image.

| Workflow | Use |
| --- | --- |
| [`txid_multi_add.wdl`](workflows/txid_multi_add.wdl) | Direct same-tool file imports; [example inputs](workflows/txid_multi_add.inputs.example.json) and [Chinese guide](docs/txid-multi-add-wdl-guide.zh-CN.md). |
| [`txid_batch.wdl`](workflows/txid_batch.wdl) | Cohort imports with per-file provenance; [example inputs](workflows/txid_batch.inputs.example.json) and [Chinese guide](docs/txid-batch-wdl-guide.zh-CN.md). |
| [Publication benchmark workflow](workflows/publication_benchmark/README.md) | Separate integration, discovery and fixed-model evaluation stages with pinned inputs and tools. |

WDL imports run within one task because isolated scatter tasks cannot safely
mutate one shared SQLite registry. Change the workflow image input to the desired
published digest for an archival run.

## Verification and scientific scope

From the repository root, run the source tests and separately implemented exact
identity calculator:

```bash
python -m unittest discover -s tests -v
python -I benchmarks/independent_identity.py conformance \
  --vectors tests/conformance/exact-v1.json
```

The [independent calculator](docs/independent-identity.md) has its own GTF reader
and canonicalization code and imports no TxID modules. It checks identity under a
supplied assembly fingerprint; it does not independently derive that fingerprint
from FASTA. This is an implementation-separation check within the project.

Evidence is reported at its evaluated scope:

- The [synthetic benchmark](docs/benchmark.md) tests identity against defined
  emitted structures, including deliberate boundary perturbations.
- The [Tier-0 report](docs/benchmark/publication-smoke-2026-07-31.md) tests workflow
  integration using actual caller and comparator executables.
- The [ENCODE fixed-model report](docs/benchmark/encode-gtf-full-2026-07-31.md)
  evaluates 538,804 observations from six released WTC11 TALON annotations.
- The [Application Note package](deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/README.md)
  includes the manuscript, supplement, figure sources and a bounded real-read
  IsoQuant–StringTie interoperability case. Its 110 shared exact forms are all
  reference-known; the case does not demonstrate shared novel forms or an advantage
  over matching shared reference aliases.

Historical evaluations labelled TxID 0.1.0 retain that version; the additional
real-caller case used 0.1.3. Its selection script,
`benchmarks/run_real_interop_callers.py`, requires **Python 3.11 or later** and
separately installed alignment/caller tools. This requirement does not change the
TxID package's Python 3.10 minimum.

Exact IDs preserve caller boundary differences. The evaluations do not establish
biological correctness, discovery sensitivity, abundance accuracy or the biological
validity of locus assignments and fuzzy groups. Large FASTA, reads, BAMs and working
databases are identified through manifests rather than bundled in the source tree.
See the [reproducibility checklist](docs/reproducibility.md) for provenance and
release requirements.

## Documentation, citation and contributions

Start with the [tutorial](docs/tutorial.md) and [CLI reference](docs/cli.md).
Implementation details are specified in the
[versioned identity specification](docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md)
and [machine-readable schemas](schemas/). Scientific comparisons and references
are included in the [manuscript package](deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/README.md).

When reporting work that uses TxID, record the software version and source commit;
[CITATION.cff](CITATION.cff) supplies repository citation metadata. The manuscript
package records the software and data availability statements for its version.

Report bugs or propose changes through
[GitHub issues](https://github.com/dawangran/TxID/issues). Include the TxID version,
command and a small synthetic reproducer where possible. Identity changes must
keep the specification, schemas, implementation and conformance tests aligned.

TxID is distributed under the [BSD 3-Clause license](LICENSE).
