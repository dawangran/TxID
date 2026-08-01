# ENCODE WTC11 full-genome GTF identity experiment — 31 July 2026

## Scope

This experiment starts from six released ENCODE TALON GTF outputs instead of
downloading much larger BAM files. It evaluates whether fixed transcript
structures receive consistent identities across samples, tools, input orders,
and annotation contexts. It does not evaluate transcript discovery from reads,
abundance estimation, or agreement among independent callers.

Compact checked artifacts are in
[`benchmarks/results/encode-gtf-full/`](../../benchmarks/results/encode-gtf-full/).
Large source GTFs, rewritten annotations, SQANTI3 products, and SQLite databases
remain outside Git.

## Public inputs

Every transfer matched the ENCODE byte count and MD5 before a local SHA-256 was
recorded. The six compressed GTFs total 302,909,514 bytes.

| Platform / experiment | Replicate | Accession | Bytes | ENCODE MD5 |
|---|---:|---|---:|---|
| PacBio capped mRNA / ENCSR309IKK | 1 | ENCFF694HZY | 45,926,211 | `d3ad12969ad078f31783955ea5df9e1d` |
| PacBio capped mRNA / ENCSR309IKK | 2 | ENCFF978SRK | 44,538,435 | `52b814a519b2dddde9fbe38430d6b94e` |
| PacBio capped mRNA / ENCSR309IKK | 3 | ENCFF973ZQF | 48,066,325 | `402e5c1cd4a76a650f920141b7f832b2` |
| ONT direct RNA / ENCSR392BGY | 1 | ENCFF527AJC | 50,774,676 | `60fe3c44bd8d27c55494b22ef2ede216` |
| ONT direct RNA / ENCSR392BGY | 2 | ENCFF615FZJ | 69,001,520 | `1471d3a04c8b2822317eddf925bd0702` |
| ONT direct RNA / ENCSR392BGY | 3 | ENCFF798LYI | 44,602,347 | `091fc224abc2ba15cf951fd153478169` |

A deterministic two-pass extractor retained TALON transcript rows and every exon
linked to a retained `transcript_id`:

| Accession | Transcript observations | Genes | Exons |
|---|---:|---:|---:|
| ENCFF527AJC | 94,130 | 16,558 | 695,873 |
| ENCFF615FZJ | 189,280 | 21,011 | 1,699,489 |
| ENCFF798LYI | 58,600 | 13,984 | 369,152 |
| ENCFF694HZY | 62,906 | 13,111 | 484,193 |
| ENCFF978SRK | 56,211 | 13,372 | 400,682 |
| ENCFF973ZQF | 77,677 | 13,623 | 599,605 |
| **Total** | **538,804** | — | **4,248,994** |

The reference sequence collection is GRCh38. GENCODE v29 and the checksum-pinned
GENCODE v49 primary-assembly annotations define separate annotation contexts on
that assembly. An NCBI assembly-report audit generated 582 aliases, resolved 45
v49 contig aliases, and left no unresolved contigs. An all-regions v49 annotation
was quarantined after context validation and never entered the registry.

## Executed paths

| Stage | Version / source | Execution |
|---|---|---|
| Exact registry | TxID 0.1.0 | Fresh sorted and shuffled v29 mappings; five persisted v49 imports plus one disclosed read-only classification; six incremental snapshots |
| Fixed-model comparison | gffcompare 0.12.10 | Fresh sorted and shuffled prefixes |
| isoSeQL prerequisite | SQANTI3 6.0.1 | All six full GTFs, four CPUs each |
| Multi-dataset database | isoSeQL 1.0.1, commit `25d9366d8b236d3b912e62dcd2c267fe6df31a4c` | Fresh sorted and shuffled databases |

Unmodified isoSeQL failed its `ends_counts` uniqueness constraint for repeated
exact-end observations. The audited compatibility source changed only that
count-table insertion to `INSERT OR IGNORE`; PBID insertion and structural lookup
code were unchanged. Both completed databases retained all 538,804 PBIDs.

## Primary identity results

Rates were scored against structure parsed independently from the fixed GTFs.
Confidence intervals use 2,000 deterministic structural-unit bootstrap
replicates with seed 20260731.

| Layer / target | Emitted / input | Groups | False merges | False-merge rate (95% CI) | False splits | False-split rate (95% CI) |
|---|---:|---:|---:|---:|---:|---:|
| TxID exact form | 538,804 / 538,804 | 467,819 | 0 | 0 (0–0) | 0 | 0 (0–0) |
| gffcompare / exact form | 538,802 / 538,804 | 379,407 | 171,218 | 0.62972 (0.62646–0.63291) | 0 | 0 (0–0) |
| gffcompare / splice chain | 538,802 / 538,804 | 379,407 | 20,458 | 0.07524 (0.07299–0.07735) | 0 | 0 (0–0) |
| isoSeQL / exact ends | 538,804 / 538,804 | 467,908 | 0 | 0 (0–0) | 93 | 0.0009237 (0.0007240–0.0011357) |
| isoSeQL / common junction | 538,804 / 538,804 | 364,032 | 90,319 | 0.26893 (0.26231–0.27528) | 5,917 | 0.02353 (0.02220–0.02488) |

gffcompare omitted two single-exon observations. Its exact-form merge rate is
high largely because 55,494 mixed clusters join transcript-end variants; the
separate splice-chain row is therefore required for a fair semantic comparison.

## Input-order and matrix effects

TxID exact form, splice-chain, transcript identifier, classification, and
gene-candidate fields had zero changes after permutation. gffcompare changed
4,438 pair relations. Both isoSeQL structural partitions were unchanged, while
their sequential internal labels changed for all 538,804 observations. The
stateful TxID `GL1` layer changed 418,224 pair relations and is not a
cross-registry structural join key.

| Layer | Matrix columns | Sparsity |
|---|---:|---:|
| TxID exact form | 467,819 | 80.804% |
| isoSeQL exact ends | 467,908 | 80.808% |
| gffcompare tracking | 379,407 | 76.331% |
| isoSeQL common junction | 364,032 | 75.731% |

The identity relation materially changes matrix dimensionality, but fewer
columns do not imply greater biological correctness.

## Annotation-release and incremental invariance

All 538,804 observations were paired between the GENCODE v29 and v49 contexts.
The exact structural identifiers were invariant while annotation-relative
interpretation changed:

| Metric | Result |
|---|---:|
| Exact form identifier changes | 0 / 538,804 |
| Exact splice-chain identifier changes | 0 / 538,804 |
| Classification changes | 29,435 / 538,804 (5.463%) |
| Gene-candidate string changes | 527,676 / 538,804 (97.935%) |
| Incremental snapshots | 6 |
| Exact identifiers disappearing across snapshots | 0 |
| Structural metadata changes across snapshots | 0 |
| Aggregate classification regressions | 0 |
| Final snapshot equals final catalog | yes, byte-for-byte |

The v49 SQLite registry persisted five sample imports. A sixth import was
automatically rolled back after the shared filesystem exhausted its quota, and a
later retry ended with a direct `disk I/O error` during its read-only staging
phase. The sixth sample was therefore evaluated with a classification-only run
that loaded the same registered v49 reference index and called the same TxID
staging, exact-identity, and classification functions without writing an
observation or allocating a locus accession. This completes the six-sample
structural and classification comparison but does not constitute a six-sample
persisted v49 registry. A later retry left a hot journal that could not be
recovered while the shared filesystem was returning `disk I/O error`; the
generated v49 database and journal were therefore removed. The five compressed
mapping outputs, reference context, classification-only audit, and compact
six-sample results were retained. The database can be regenerated from the
checked inputs when storage is stable.

The first four incremental catalogs reconstructed from checked component
mappings were byte-identical to the corresponding direct SQLite exports. The
same deterministic aggregation generated snapshots five and six, and snapshot
six was byte-identical to the final sorted catalog.

## Storage recovery note

Shared storage reached its quota during the long registry outputs. Only
reproducible generated intermediates were removed: redundant combined
gffcompare GTFs, uncompressed copies with checked gzip counterparts, SQANTI3
corrected FASTA/GTF and junction intermediates, duplicate observation tables,
and completed comparator/order SQLite databases whose compact mappings and
scored results had already been retained. The unrecoverable v49 SQLite file and
its hot journal were also removed only after verifying the five persisted
mappings and six-sample compact comparison. Source ENCODE GTFs, reference inputs,
checked compressed mappings, classifications, primary metrics, and compact
audits were preserved.

## Interpretation boundary

This experiment supports deterministic identity and interoperability claims for
released transcript models. It does not evaluate raw-read discovery sensitivity,
abundance accuracy, or full-genome caller superiority. A benchmark in which the
same raw reads are processed independently by IsoQuant, FLAIR, StringTie, and
TALON would address those separate questions and remains future work outside the
claims of the current fixed-model Technical Note.
