# TxID publication benchmark workflow

This Snakemake workflow implements the benchmark protocol in
`docs/superpowers/specs/2026-07-31-publication-benchmark-design.md`.
It keeps transcript discovery and transcript identity as separate evidence
tracks.

## Profiles

- `config/config.smoke.yaml` generates a tiny truth-defined reference and
  full-length reads, then invokes real executables. It is an integration smoke
  test, not publication evidence.
- `config/config.public.template.yaml` defines the intended two-platform,
  three-replicate public run. It remains non-runnable while manifest rows are
  disabled and reference checksums are `NA`.
- `config/config.public.pacbio.yaml` is the checked WTC11 PacBio triplicate run
  profile. It becomes runnable only after staging produces its checksum-bearing
  dataset and reference manifests; it is an interim platform-specific result,
  not a substitute for the two-platform analysis.
- `config/config.public.pacbio-preflight.yaml` runs the same checked pipeline on
  the first two fully staged PacBio replicates. It is a preflight for detecting
  workflow failures while the third transfer is pending, not manuscript
  evidence.
- `config/datasets.public.tsv` and `config/references.public.tsv` are the public
  data registry. Rows whose local inputs have not been staged remain disabled.
  Enabling a row without its accession, checksum, and local path is a validation
  error.

Generated work belongs under `benchmark-work/`, which is ignored by Git. Large
FASTQ, BAM, databases, and caller outputs must not be committed.

The checked 31 July 2026 Tier-0 execution is summarized in
`docs/benchmark/publication-smoke-2026-07-31.md`. It includes real SQANTI3
preprocessing, isoSeQL ingestion, and independent sorted/shuffled TALON
databases. Its numerical results must not be promoted to biological-performance
claims.

The separate six-file ENCODE WTC11 full-genome GTF experiment is summarized in
`docs/benchmark/encode-gtf-full-2026-07-31.md`. It avoids raw BAM transfer by
starting from official ENCODE TALON annotations, verifies every byte count and
MD5, and evaluates 538,804 fixed models through TxID, gffcompare, and the
complete SQANTI3-to-isoSeQL path. It is public fixed-model evidence, not a
multi-caller discovery benchmark.

## Install and inspect

Create a small workflow-controller environment:

```bash
conda env create -p /tmp/txid-snakemake \
  -f workflows/publication_benchmark/envs/snakemake.yaml

conda run -p /tmp/txid-snakemake snakemake \
  --snakefile workflows/publication_benchmark/Snakefile \
  --configfile workflows/publication_benchmark/config/config.smoke.yaml \
  --lint
```

Render the DAG and perform a dry run:

```bash
conda run -p /tmp/txid-snakemake snakemake \
  --snakefile workflows/publication_benchmark/Snakefile \
  --configfile workflows/publication_benchmark/config/config.smoke.yaml \
  --dag > benchmark-work/smoke/dag.dot

conda run -p /tmp/txid-snakemake snakemake \
  --snakefile workflows/publication_benchmark/Snakefile \
  --configfile workflows/publication_benchmark/config/config.smoke.yaml \
  --use-conda --cores 8 --dry-run
```

Run the checked smoke target:

```bash
conda run -p /tmp/txid-snakemake snakemake \
  --snakefile workflows/publication_benchmark/Snakefile \
  --configfile workflows/publication_benchmark/config/config.smoke.yaml \
  --use-conda --cores 8 publication_smoke
```

The publication profile is intentionally not enabled by default. Copy the
template to a run-specific configuration, stage and checksum its public inputs,
change only the required manifest rows to `enabled=true`, validate the manifest,
and archive the generated provenance bundle with the manuscript data. The public
workflow derives SIRV transcript sequences from the checked GTF and assigns only
uniquely mapping, near-full-length reads after predeclared platform-specific
identity, coverage, and score-margin filters. Assignment attrition is retained
in `truth/sirv-read-truth.run.json`.

For offline or rate-limited execution, `sources.sqanti3.cache` and
`sources.isoseql.cache` may name local Git repositories containing the requested
revisions. The fetch audit continues to record the declared upstream URL,
requested revision, resolved commit, retrieval source, and whether a cache was
used; a cache never changes the pinned upstream identity.

Stage ENCODE inputs through the checked downloader. It reads official ENCODE
JSON metadata, resumes partial transfers, checks ENCODE's byte count and MD5,
computes SHA-256, atomically promotes each completed file, and writes a new
manifest rather than editing the source template:

```bash
python workflows/publication_benchmark/scripts/stage_encode_files.py \
  --manifest workflows/publication_benchmark/config/datasets.public.tsv \
  --metadata-dir benchmark-work/public-metadata \
  --repo . \
  --accession ENCFF105WIJ --accession ENCFF212HLP \
  --accession ENCFF003QZT \
  --download \
  --output-manifest benchmark-work/public/config/datasets.pacbio.tsv \
  --audit benchmark-work/public/metadata/encode-pacbio-staging.json
```

Use `--refresh-metadata` when an online refresh is required. Without
`--download`, the command only validates and archives metadata. A failed or
interrupted transfer remains as `*.part` and can be resumed by repeating the
same command; it is never treated as an enabled input.

### GTF-first ENCODE identity benchmark

The identity/interoperability track can use released ENCODE transcript-model
GTFs without downloading read-level BAMs. The checked registries
`config/datasets.encode-gtf.public.tsv` and
`config/references.encode-gtf.public.tsv` declare three WTC11 PacBio and three
WTC11 Oxford Nanopore replicates, the ENCODE GRCh38 reference, and GENCODE v29
and v49 annotation contexts.
Their six compressed GTFs total approximately 303 MB, compared with the much
larger alignment files.

These ENCODE files combine the GENCODE reference annotation and sample-specific
TALON models. `extract_encode_talon_models.py` therefore performs a two-pass,
audited extraction: it first selects transcript rows whose source column is
`TALON`, then retains every exon belonging to those transcript IDs even when an
exon itself is labelled `HAVANA` or `ENSEMBL`. The six checked files contain
538,804 TALON transcript observations. `make_gtf_identity_view.py` can then
produce a deterministic structure-and-identifier-only view for tools that reject
legal repeated nonidentity GTF attributes; the full source GTF and its hashes
remain unchanged.

This route supports real biological comparisons of transcript identity,
cross-replicate registry behavior, sample-by-transcript matrices, gffcompare,
and the SQANTI3-to-isoSeQL path. It does **not** constitute a new local TALON run:
the models are official TALON-derived ENCODE outputs. Read-level TALON discovery,
alignment accuracy, and abundance claims still require BAM or FASTQ inputs and
must be reported as a separate evidence track.

`scripts/classify_fixed_models.py` is a read-only recovery route. It loads an
existing annotation context and calls the production parser, exact-identity,
reference-index, and classification code without persisting observations or
allocating `GL1` accessions. Its outputs must be labelled
`classification_only_no_observation_persistence` and must not be represented as
a completed registry import.

Stage the LRGASP/GENCODE references in the same way:

```bash
python workflows/publication_benchmark/scripts/stage_public_references.py \
  --sources workflows/publication_benchmark/config/references.public.sources.tsv \
  --references workflows/publication_benchmark/config/references.public.tsv \
  --repo . --download \
  --output-manifest benchmark-work/public/config/references.checked.tsv \
  --audit benchmark-work/public/metadata/reference-staging.json
```

The source registry pins Zenodo LRGASP files by byte count and MD5 and pins the
GENCODE v49 chromosome annotation using GENCODE's release `MD5SUMS`. The tool
constructs `gencode-v49-plus-sirv4.gtf.gz` with deterministic gzip metadata,
computes SHA-256 for both annotation contexts, and enables a new reference
manifest without changing the checked template.

## Result boundaries

The workflow emits:

- sequence-backed truth, FASTQ, BAM, alignment statistics, and checksums;
- real IsoQuant, FLAIR, StringTie, and sorted/shuffled shared-database TALON
  outputs;
- TxID registries for sorted, shuffled, incremental, and annotation-release
  comparisons;
- gffcompare and full SQANTI3-to-isoSeQL runs over fixed caller models;
- truth, identity, order, annotation, matrix, runtime, and memory tables;
- a software/input/output provenance manifest.

TALON 6.0.1 accepts SAM or BAM input. The checked runner passes the
coordinate-sorted BAMs directly and records this choice in provenance, avoiding
large redundant uncompressed SAM intermediates. TALON still creates its own
sorted and merged BAMs below the rule-specific work directory as required by its
supported workflow; source FASTQ, BAM, annotations, and result artifacts are
never removed.

Discovery precision, recall, and F1 include deterministic transcript-structure
bootstrap intervals. Pairwise false-merge intervals resample output clusters and
false-split intervals resample exact truth structures. Both profiles record 2,000
replicates and seed 20260731 in their configuration and JSON output.

No rule modifies an input annotation or reference in place.
