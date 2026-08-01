# Synthetic evaluation and external benchmark plan

## Publication workflow

The current publication protocol is implemented as a separate Snakemake workflow
under `workflows/publication_benchmark/`. Unlike the legacy structure-only
simulation described below, its Tier-0 profile starts from sequence-backed reads
and actually invokes IsoQuant, FLAIR, StringTie, TALON, gffcompare, SQANTI3, and
isoSeQL. The checked execution and its evidence boundary are recorded in
`docs/benchmark/publication-smoke-2026-07-31.md`.

The public WTC11 PacBio/ONT triplicate manifest and the two-annotation-context
template are checked in but disabled. They become runnable only after local
staging, SHA-256 recording, and audited SIRV read-truth preparation. The Tier-0
numbers must not be substituted for that public experiment.

## Reproducible checked evaluation

Run:

```bash
PYTHONPATH=src python benchmarks/run_simulation.py \
  --seed 20260730 --transcripts 120 --samples 6 --tools 4 \
  --artifact-rate 0.05 --output-dir benchmarks/results
```

The simulator creates one synthetic assembly, two annotation releases, and
independent GTF/GFF3 outputs labelled IsoQuant-like, TALON-like, FLAIR-like, and
StringTie-like. These labels indicate serialization/provenance diversity only;
they are not discovery outputs from those packages. Five percent of observations
receive a 1-bp boundary perturbation. Real integration-tool executions are
described below.

Metrics distinguish three questions:

1. Structural truth asks whether identical emitted exon structures split or
   different emitted structures merge. These are the exact-identity false-split
   and false-merge rates.
2. Underlying biological truth treats perturbed outputs as observations of the
   same source transcript. Fragmentation here measures injected caller-model
   disagreement, which exact TxID is designed to preserve rather than hide.
3. The matrix result counts the union of upstream identifiers versus exact TxID
   columns after aggregating simulated counts.

The scaling rows reuse synthetic models and are explicitly labelled
`synthetic_resampled_scaling`; they are not independent biological samples.

## Checked result

With the default seed, 2,511 observations produced 0 exact false-split pairs and
0 false-merge pairs among 22,720 eligible pairs. Original upstream identifiers
split all 22,720 same-structure pairs. The matrix decreased from 2,511 raw ID
columns to 190 exact-structure columns. Of 109 observations checked across
annotation releases, every exact TxID was invariant and 22 annotation-relative
classifications changed. End-to-end execution took approximately 3.5 s with
5.47 MiB peak Python-traced memory on the development host; this single-host value
is not a cross-tool performance claim. Identity-only resampling processed 120,
240, and 480 models at approximately 11,191, 11,422, and 11,368 models per second,
respectively.

## Checked external comparison

The 31 July 2026 run used gffcompare 0.12.10, isoSeQL tag 1.0.1 at commit
`25d9366d8b236d3b912e62dcd2c267fe6df31a4c`, and TALON 6.0.1 with PyRanges
0.0.129 pinned for Python 3.7 compatibility. Inputs were generated with:

```bash
PYTHONPATH=src python benchmarks/prepare_external_comparison.py \
  --output-dir "$TXID_EXTERNAL_INPUTS" \
  --work-dir "$TXID_EXTERNAL_WORK"

cd "$TXID_EXTERNAL_INPUTS"
sha256sum --check SHA256SUMS
```

The bundle contains 24 datasets, 2,511 observations, 190 complete emitted forms,
and 162 multi-exon splice chains. The real packages consume these neutral,
truth-labelled inputs; the `-like` labels are never presented as real discovery
outputs.

gffcompare was run against the same reference and all 24 full GTF files. The
checked runner records the exact file order and repeats the run after a fixed
permutation:

```bash
python benchmarks/run_gffcompare_comparison.py \
  --bundle "$TXID_EXTERNAL_INPUTS" \
  --gffcompare "$TXID_GFFCOMPARE" \
  --output "$TXID_EXTERNAL_RUNS/gffcompare/sorted" \
  --order sorted

python benchmarks/run_gffcompare_comparison.py \
  --bundle "$TXID_EXTERNAL_INPUTS" \
  --gffcompare "$TXID_GFFCOMPARE" \
  --output "$TXID_EXTERNAL_RUNS/gffcompare/shuffled" \
  --order shuffled
```

isoSeQL was run once per dataset into one shared database. The checked fixtures
satisfy the SQANTI3 classification and genePred schemas, but SQANTI3 itself was
not executed or timed:

```bash
python benchmarks/run_isoseql_comparison.py \
  --bundle "$TXID_EXTERNAL_INPUTS" \
  --isoseql-source "$TXID_ISOSEQL_SOURCE" \
  --python "$TXID_ISOSEQL_PYTHON" \
  --output "$TXID_EXTERNAL_RUNS/isoseql/sorted" \
  --order sorted
```

TALON used a database freshly initialized from `reference.gtf`, then annotated
the 24 oriented exact-match SAM files in the supported shared-database workflow:

```bash
talon_initialize_database \
  --f "$TXID_EXTERNAL_INPUTS/reference.gtf" \
  --g TxID-simulation-v1 --a ref-v1 --idprefix TXSIM \
  --o "$TXID_TALON_RUN/shared"

talon \
  --f "$TXID_EXTERNAL_INPUTS/talon-config.csv" \
  --db "$TXID_TALON_RUN/shared.db" \
  --build TxID-simulation-v1 --threads 4 \
  --cov 0.9 --identity 0.8 --nsg \
  --tmpDir "$TXID_TALON_RUN/tmp" --verbosity 1 \
  --o "$TXID_TALON_RUN/run"

talon_abundance \
  --db "$TXID_TALON_RUN/shared.db" \
  -a ref-v1 -b TxID-simulation-v1 \
  --o "$TXID_TALON_RUN/abundance"
```

A second TALON database used `talon-config-shuffled.csv`. The evaluator compares
group membership rather than sequential internal identifiers and writes compact
checked results:

```bash
python benchmarks/evaluate_external_comparison.py \
  --bundle "$TXID_EXTERNAL_INPUTS" \
  --runs "$TXID_EXTERNAL_RUNS" \
  --output benchmarks/results
```

| Workflow | Exact-form columns | False-split rate | False-merge rate | Partition stable after shuffle |
|---|---:|---:|---:|---|
| TxID | 190 | 0 | 0 | yes |
| gffcompare | 180 | 0 | 0.015043 | yes |
| isoSeQL | 190 | 0 | 0 | yes |
| TALON | 180 | 0 | 0.015043 | yes |

gffcompare and TALON merged exact-end variants under their tested semantics; all
three external workflows recovered the 162 multi-exon splice-chain partition
without false splits or merges. Core-stage wall time and peak RSS were 0.03 s and
8.68 MiB for gffcompare, 1.23 s and 13.31 MiB for isoSeQL ingestion, and 2.68 s
and 87.94 MiB for TALON annotation. These stages perform different work and are
not an interchangeable speed ranking. TALON reference initialization and
abundance export and the omitted SQANTI3 prerequisite are reported separately in
`benchmarks/results/external-comparison-summary.json`.
