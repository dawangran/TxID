# Publication benchmark Tier-0 execution report

**Execution date:** 31 July 2026  
**Evidence level:** Tier 0 engineering smoke test  
**Publication status:** passed as an executable integration test; insufficient as
standalone manuscript evidence

This report records the first end-to-end execution of the publication benchmark
protocol. It used deterministic, sequence-backed synthetic reads but invoked the
real alignment, transcript-discovery, comparison, quality-control, registry, and
database executables. Generated data and databases remain under the ignored
`benchmark-work/smoke/` directory.

## Executed software

| Stage | Executed version |
|---|---|
| Alignment | minimap2 2.31-r1302; samtools 1.23.1 |
| Transcript discovery | IsoQuant 3.13.0; FLAIR 3.0.0; StringTie 3.0.3 |
| Shared-database discovery | TALON 6.0.1 |
| Fixed-model comparison | gffcompare 0.12.10 |
| isoSeQL prerequisite | SQANTI3 6.0.1 |
| Multi-dataset database | isoSeQL 1.0.1, commit `25d9366d8b236d3b912e62dcd2c267fe6df31a4c` |
| Identity registry | TxID 0.1.0 |

SQANTI3 was executed on all seven fixed caller annotations before isoSeQL
ingestion. TALON was executed twice using fresh shared databases: sample order
`smoke-pb-r1, smoke-pb-r2` and the fixed permutation
`smoke-pb-r2, smoke-pb-r1` (seed 20260731).

## Discovery-track results

The truth set contained eight exact transcript forms and five multi-exon splice
chains. These small results verify wiring and metric calculation; they are not
estimates of caller performance.

| Caller | Unit | Exact-form F1 | Splice-chain F1 |
|---|---|---:|---:|
| IsoQuant | sample 1 / sample 2 | 0.615 / 0.667 | 1.000 / 0.889 |
| FLAIR | sample 1 / sample 2 | 0.769 / 0.769 | 0.889 / 0.889 |
| StringTie | sample 1 / sample 2 | 0.857 / 0.769 | 1.000 / 0.889 |
| TALON | one shared database | 0.857 | 1.000 |

The machine-readable source is
`benchmark-work/smoke/reports/discovery-summary.tsv`.

## Fixed-model identity comparison

Seven real caller annotations contributed 36 fixed transcript observations. The
target was independently parsed equality of contig, strand, and all exon
intervals. False-merge rate is the fraction of within-output-cluster observation
pairs that have different target structures; false-split rate is the fraction of
same-target observation pairs assigned to different output clusters.

| Layer | Target groups | Output groups | False merges | False splits |
|---|---:|---:|---:|---:|
| TxID exact form | 7 | 7 | 0/95 (0) | 0/95 (0) |
| gffcompare tracking | 7 | 6 | 2/97 (0.0206) | 0/95 (0) |
| isoSeQL `isoform_ends` | 7 | 7 | 0/95 (0) | 0/95 (0) |
| isoSeQL common junction | 6 | 6 | 0/97 (0) | 0/97 (0) |

The single gffcompare mixed cluster contained transcript-end variants with one
common splice chain; no splice structures were mixed. This is a difference in
the tested equivalence relation, not evidence that gffcompare is defective.

Randomizing the seven fixed input files changed no pairwise group relation for
gffcompare or isoSeQL. gffcompare's tested cluster labels also remained fixed.
isoSeQL preserved both structural partitions but 14 of 36 observations received
different internal numeric identifiers at each of its exact-end and
common-junction layers. isoSeQL does not claim content-derived permanent
identifiers; TxID addresses that separate interoperability requirement.

## TALON shared-database comparison

TALON was evaluated on its supported BAM-to-shared-database workflow rather than
being represented as a fixed-GTF registry. All 194 aligned reads were present in
both fresh databases.

| Target | Truth groups | TALON groups | False merges | False splits |
|---|---:|---:|---:|---:|
| Exact transcript form | 8 | 6 | 1,500/3,867 (0.3879) | 0/2,367 (0) |
| Splice chain | 6 | 6 | 0/3,867 (0) | 0/3,867 (0) |

The two mixed exact-form groups are the two truth-defined transcript-end variant
pairs. Reversing sample order changed neither read-pair relations nor TALON
transcript identifiers in this smoke fixture. These read-level denominators are
not directly comparable with the 36-model fixed-GTF table.

## TxID invariance and matrix consequence

- All 36 fixed observations retained their exact form and splice-chain TxIDs
  after the fixed input permutation.
- Sorted and shuffled catalog exports were byte-identical with SHA-256
  `8c78f4b3ffa7538ce4061b97915b8131950ee9795123d393d13f53d0cdfd77d3`.
- Seven incremental snapshots had no identifier disappearance, structural
  metadata change, or classification regression; the last snapshot equalled the
  final catalog byte for byte.
- The second annotation context changed two classifications and zero exact TxIDs.
- Excluding TALON's pooled shared export from the biological matrix, 30
  per-sample caller observations formed 16 upstream identifier columns and seven
  TxID exact-form columns across two samples.

## Reproducibility artifacts

The workflow is defined by
[Snakefile](../../workflows/publication_benchmark/Snakefile) and its
[README](../../workflows/publication_benchmark/README.md). The ignored execution
directory contains command vectors, logs, wall time, peak RSS, input/output
checksums, databases, and these checked reports:

| Artifact | SHA-256 |
|---|---|
| `discovery-summary.json` | `e3303e21e3f86a7a11cf665c85d4367b0a07b886f62e8b8936261e5ddafef10b` |
| `identity-partition-summary.json` | `850b167c29c100ac8c6871710c1f2c516f790abfd8c0dda22f1587f2e3c21a73` |
| `talon-partition-summary.json` | `cd306b8bf7fcca333a5d272c79539ebbd691c8199ee7634327e3b807560d3834` |

The provenance bundle hash is intentionally not fixed in this document because
it changes whenever a new execution record is added. Publication runs must
archive the whole bundle and environment locks under a persistent identifier.

## Deterministic bootstrap addendum

The retained caller models and databases were re-evaluated without rerunning or
overwriting the original external-tool results. Point estimates were unchanged.
The addendum uses 2,000 percentile-bootstrap replicates and seed 20260731.
Discovery intervals resample unique truth and predicted structures;
false-merge intervals resample output clusters; false-split intervals resample
truth structural groups. These Tier-0 intervals quantify structural-unit
variation in this fixture, not patient-level biological variation.

| Evaluation | Point estimate | 95% bootstrap interval |
|---|---:|---:|
| gffcompare exact-form false-merge rate | 0.0206 | 0–0.0984 |
| TxID exact-form false-merge rate | 0 | 0–0 |
| isoSeQL exact-end false-merge rate | 0 | 0–0 |
| TALON exact-form false-merge rate | 0.3879 | 0–0.4810 |
| All tested false-split rates | 0 | 0–0 |

The additive outputs are retained under
`benchmark-work/smoke/reports-bootstrap/`. Their JSON SHA-256 values are
`998e4284f8dbdf4e380bf438c48e4bcb44c3f363c39a943d119cab1b3800d9b7`
(discovery),
`b6e9e769960c0a7f4db89922e802490218244e4ce1cd8b127ffb228875aac70f`
(fixed-model identity), and
`27631753c172f7ce54f7e873506e882c71ac4f1ac8bac320838c898ca55b0b34`
(TALON).

## Remaining publication gates

This run does not support submission by itself. Promotion from pilot requires:

1. staging the disabled LRGASP WTC11 PacBio and ONT triplicates in the public
   manifest and recording downloaded-file checksums;
2. building checked GRCh38+SIRV reference contexts for GENCODE v38 and v49;
3. deriving and auditing an unambiguous SIRV read-truth table for the read-level
   TALON analysis;
4. executing the complete workflow on the public reads, with failures retained;
5. bootstrap confidence intervals over truth transcripts or biological
   replicates, plus public archive, container digest, and source-data tables;
6. regenerating manuscript figures and the website only from the archived
   public TSV/JSON outputs.

The website is therefore a later presentation layer, not a substitute for the
missing public experiment.
