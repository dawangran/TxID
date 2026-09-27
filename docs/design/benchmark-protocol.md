# TxID publication benchmark design

Status: implementation protocol 0.1  
Fixed random seeds: 20260730 (simulation), 20260731 and 20260732 (orders)

## Question and evidence tracks

The benchmark must not confuse transcript discovery with transcript identity.
It therefore has two linked, separately reported tracks.

1. **Discovery track.** The same truth-defined or public long-read data are
   aligned once and supplied to IsoQuant, FLAIR, StringTie, and TALON. Exact
   transcript-form and splice-chain precision/recall, end error, redundancy,
   runtime, and peak memory are measured against truth where truth exists.
2. **Identity track.** Per-sample GTF outputs from the discovery tools are held
   fixed. TxID, gffcompare, and SQANTI3 followed by isoSeQL consume those fixed
   models. TALON is evaluated through its supported shared-database BAM workflow,
   not misrepresented as a GTF registry. False splits, false merges, order
   sensitivity, incremental behavior, annotation-release behavior, and effects
   on sample-by-transcript matrices are measured.

Runtime values from stages that perform different work are never presented as a
speed ranking. Tool setup, alignment, discovery, SQANTI3 preprocessing, registry
ingestion, and export are timed separately.

## Dataset tiers

### Tier 0: deterministic smoke test

The repository generates a small sequence-backed transcriptome, two annotation
releases, two samples, and error-bearing full-length reads. These reads pass
through minimap2 and the real caller executables. This tier catches installation,
CLI, format, and provenance failures. It is not biological publication evidence.

### Tier 1: public truth-defined controls

The primary truth tier uses LRGASP SIRV Set 4 controls and the associated public
long-read runs. The manifest records accessions, reference/annotation versions,
checksums, platform, read type, replicate, and source URLs. This tier estimates
exact false-split/false-merge and discovery accuracy with known structures.

### Tier 2: public biological data

At least two biological replicates from one PacBio and one Oxford Nanopore
experiment are processed by every compatible caller. K562/ENCODE is the initial
human cohort. This tier measures cross-caller overlap, incremental stability,
matrix dimensionality, sparsity, and annotation-relative classification. It does
not claim biological truth for unvalidated novel transcripts.

Resampled or duplicated inputs are labelled scaling data and never counted as
independent biological samples.

### Auxiliary public fixed-model audit

Released full-genome caller GTFs may be used to test identity interoperability
without transferring the much larger raw-read or alignment files. TxID,
gffcompare, and the complete SQANTI3-to-isoSeQL path consume the fixed models.
The source caller is recorded as provenance and is not rerun or presented as an
independent GTF registry. This audit can support claims about fixed-model
identity, input order, annotation context, and matrix consequences, but it does
not measure discovery or abundance and does not satisfy the Tier 1 or Tier 2
publication gates.

## Versioned tools and execution contracts

The initial environment contracts pin:

- minimap2 2.31 and samtools 1.21;
- IsoQuant 3.13.0;
- FLAIR 3.0.0;
- StringTie 3.0.3;
- TALON 6.0.1 with Python 3.7 and PyRanges 0.0.129;
- gffcompare 0.12.10;
- SQANTI3 6.0.1;
- isoSeQL 1.0.1 at commit
  `25d9366d8b236d3b912e62dcd2c267fe6df31a4c`;
- Snakemake 9.23.1.

Each rule writes a machine-readable execution record containing the exact
argument vector, versions, start/end times, exit status, input/output checksums,
wall time, and peak resident memory when available. Environment YAML files are
part of the workflow; publication runs must additionally archive explicit conda
package lists or immutable container digests.

Long-read alignments are coordinate sorted, indexed, and include `MD` tags so
the same BAM-derived SAM can be used by TALON. SIRV mapping disables splice-flank
priors. Platform-specific minimap2 presets and caller data-type flags are selected
from the dataset manifest and recorded.

## Required perturbations and comparisons

For every truth tier and for the biological tier where applicable:

- import the same files in lexicographic order and in at least two fixed random
  orders;
- compare final group membership, exact labels, and byte-level TxID catalogs;
- import incrementally and retain matrix/catalog snapshots after every dataset;
- run the same emitted structures against two annotation releases on one
  assembly;
- report exact forms separately from splice chains and from optional fuzzy
  clusters;
- compare raw upstream-ID, TxID-form, splice-chain, gffcompare, isoSeQL, and TALON
  matrix column counts and sparsity;
- preserve ambiguous bridges, unmatched models, malformed outputs, and tool
  failures instead of silently dropping them.

## Metrics

Truth-defined metrics include exact-form and splice-chain precision, recall, F1,
false-split pairs, false-merge pairs, unmatched predictions, missed truth forms,
end-coordinate error among splice-matched forms, and single-exon results reported
separately.

Order and release checks include changed exact identifiers, changed partitions,
changed classifications, deterministic catalog SHA-256, and snapshot consistency.
Matrix metrics include dimensions, nonzero fraction, duplicated columns,
sample-pair similarity, and the number of entries changed by identifier
harmonization.

Bootstrap confidence intervals are computed over truth transcripts or biological
samples, never over duplicated reads. Seeds and the resampling unit are recorded.

## Publication gate

Scope decision, 31 July 2026: the gate below governs promotion of a broad
raw-read discovery and caller-comparison claim. The current Technical Note is
narrowed to deterministic identity and interoperability of released transcript
models, supported by the six-sample full-genome ENCODE GTF experiment. It does
not promote or imply discovery sensitivity, abundance accuracy, or full-genome
caller superiority. The Tier 2 raw-read experiment therefore remains required
for those broader future claims, not for the fixed-model claim.

A manuscript result may be promoted from “pilot” only after:

- Tier 0 passes with every named executable actually invoked;
- Tier 1 includes at least two sequencing platforms or explicitly narrows the
  claim to one;
- Tier 2 includes independent biological replicates and at least three upstream
  callers on identical alignments;
- fixed-order and incremental tests pass;
- two real annotation releases on the same assembly are evaluated;
- SQANTI3 preprocessing is included in isoSeQL timing and provenance;
- all plotted values trace to checked TSV/JSON source data;
- raw data accessions, input checksums, environment locks, logs, and failed-run
  accounting are archived.

The website is a presentation and exploration layer over released artifacts. It
must not be used as the only archive of benchmark data or software.
