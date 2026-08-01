# ENCODE WTC11 GTF identity pilot — 31 July 2026

> Superseded for primary reporting by the full-genome experiment in
> [`encode-gtf-full-2026-07-31.md`](encode-gtf-full-2026-07-31.md). This file is
> retained as the chromosome-restricted preflight record.

## Scope

This is a public fixed-model identity experiment over chr21 and chr22. It starts
from official ENCODE TALON GTF outputs to avoid the much larger BAM transfers.
It does not measure transcript discovery accuracy, full-genome performance, or
agreement among independent callers.

The checked compact artifacts are in
[`benchmarks/results/encode-gtf-pilot/`](../../benchmarks/results/encode-gtf-pilot/).
Large source GTFs, SQANTI3 products, rewritten annotations, and SQLite databases
remain outside Git.

## Public inputs

All six files were released ENCODE GTFs. The checked downloader matched the
official byte count and MD5, then recorded SHA-256.

| Platform / experiment | Replicate | Accession | Bytes | ENCODE MD5 | SHA-256 |
|---|---:|---|---:|---|---|
| PacBio capped mRNA / ENCSR309IKK | 1 | ENCFF694HZY | 45,926,211 | `d3ad12969ad078f31783955ea5df9e1d` | `c34122b42d1b7c7b30fcc97edc545db99c42f01ad3cb6d2f2c4e32e144a66859` |
| PacBio capped mRNA / ENCSR309IKK | 2 | ENCFF978SRK | 44,538,435 | `52b814a519b2dddde9fbe38430d6b94e` | `2a4135b10364a7375c0e0584de884b01f9a5cc69cc371274288d79490c0af794` |
| PacBio capped mRNA / ENCSR309IKK | 3 | ENCFF973ZQF | 48,066,325 | `402e5c1cd4a76a650f920141b7f832b2` | `d359f0426d74f87d7b8204a0e5a894d3fafa388a6d373b3464e780fd50afbed8` |
| ONT direct RNA / ENCSR392BGY | 1 | ENCFF527AJC | 50,774,676 | `60fe3c44bd8d27c55494b22ef2ede216` | `446295ff1c911a41a98a5ebedd4bf77bc509ba68aee936de4a3f37b3b660a8b8` |
| ONT direct RNA / ENCSR392BGY | 2 | ENCFF615FZJ | 69,001,520 | `1471d3a04c8b2822317eddf925bd0702` | `9a8496858f4473b3168c084931c8361128d5a19fc45dcb034b8d67a5b4260ac0` |
| ONT direct RNA / ENCSR392BGY | 3 | ENCFF798LYI | 44,602,347 | `091fc224abc2ba15cf951fd153478169` | `8bea22c31380df842c5ab452af0bbe7f66281816d6728d63fcf525924196af22` |

The reference was GRCh38 with GENCODE v29 primary-assembly annotation using UCSC
contig names. The deterministic chr21+chr22 subsets contained:

| Accession | Transcript observations |
|---|---:|
| ENCFF527AJC | 10,374 |
| ENCFF615FZJ | 14,057 |
| ENCFF694HZY | 9,228 |
| ENCFF798LYI | 9,125 |
| ENCFF973ZQF | 9,739 |
| ENCFF978SRK | 8,986 |
| **Total** | **61,509** |

## Audited comparison view

TxID v1 rejects repeated keys on one GTF feature. GENCODE v29 contains legitimate
multi-valued non-identity attributes, so the benchmark generated a structure-only
view with:

- transcript and exon features only;
- coordinate columns preserved;
- original `gene_id` and `transcript_id` preserved;
- every omitted repeated attribute occurrence counted in an audit.

The GENCODE view recorded 50,644 repeated occurrences: 50,305 `tag` and 339
`ont`. The six TALON views contained no repeated attributes. Raw files were never
modified.

## Executed software paths

| Stage | Version / source | Execution |
|---|---|---|
| Exact registry | TxID 0.1.0 | Fresh sorted and shuffled registries |
| Fixed-model comparison | gffcompare 0.12.10 | Fresh sorted and shuffled prefixes |
| isoSeQL prerequisite | SQANTI3 6.0.1 | All six GTFs, 4 CPUs each |
| Multi-dataset database | isoSeQL 1.0.1, commit `25d9366d8b236d3b912e62dcd2c267fe6df31a4c` | Fresh sorted and shuffled databases |

Observed aggregate wall times on this host were 74.39 s for the sorted TxID run,
1.59 s for sorted gffcompare, 325.96 s for all six SQANTI3 runs, and 4.70 s for
sorted isoSeQL ingestion. Stage scope differs, concurrent SQANTI3 execution
affects aggregation, and these values are not a speed ranking.

## Unmodified isoSeQL failure and compatibility evaluation

Unmodified isoSeQL failed on the first public input:

```text
sqlite3.IntegrityError: UNIQUE constraint failed:
ends_counts.ends_id, ends_counts.exp, ends_counts.read_count
```

Different TALON transcript names can share one exact end structure within an
experiment, contrary to the assumption in isoSeQL's count insertion. The failed
databases and logs were retained. To evaluate isoSeQL's structural partitions,
an audited compatibility source changed only the non-single-cell `ends_counts`
statement from `INSERT` to `INSERT OR IGNORE`. PBID insertion and structural
queries were unchanged. Both resulting databases retained all 61,509 PBIDs.

## Identity results

Exact truth was equality of contig, strand, and all exon intervals parsed from
the fixed input GTFs. Confidence intervals used 2,000 structural-unit bootstrap
replicates and seed 20260731.

| Layer | Emitted / input | Groups | False merges | False-merge rate (95% CI) | False splits | False-split rate (95% CI) |
|---|---:|---:|---:|---:|---:|---:|
| TxID exact form | 61,509 / 61,509 | 23,338 | 0 | 0 (0–0) | 0 | 0 (0–0) |
| gffcompare tracking | 61,087 / 61,509 | 20,493 | 5,475 | 0.04779 (0.04511–0.05086) | 0 | 0 (0–0) |
| isoSeQL exact ends | 61,509 / 61,509 | 23,338 | 0 | 0 (0–0) | 0 | 0 (0–0) |
| isoSeQL common junction | 61,509 / 61,509 | 19,831 | 3,546 | 0.02935 (0.02528–0.03367) | 119 | 0.001014 (0.000660–0.001431) |

gffcompare omitted 422 fixed models from tracking. Its rates use only emitted
members and the omission is not imputed.

## Input-order behavior

| Layer | Pair relations changed | Literal identifiers changed |
|---|---:|---:|
| TxID exact form | 0 | 0 |
| gffcompare tracking | 142 | 60,919 |
| isoSeQL exact ends | 0 | 61,492 |
| isoSeQL common junction | 0 | 61,492 |
| TxID GL gene/locus | 11,413 | 6,009 |

The TxID exact form, splice-chain, transcript identifier, classification, and
gene-candidate fields were invariant. `GL1` is a persistent registry accession,
not a content-derived exact ID, and must not be used as a cross-registry join key.

## Matrix consequence

| Layer | Columns | Occupied cells | Sparsity | Columns in multiple samples |
|---|---:|---:|---:|---:|
| TxID exact form | 23,338 | 61,503 | 56.1% | 8,994 |
| gffcompare tracking | 20,493 | 61,087 | 50.3% | 10,160 |
| isoSeQL exact ends | 23,338 | 61,503 | 56.1% | 8,994 |
| isoSeQL common junction | 19,831 | 60,763 | 48.9% | 10,242 |

The smaller junction/end-collapsed matrices are not automatically more
biologically correct; they implement different identity relations.

## Interpretation boundary

This experiment supports the claim that TxID exact identifiers reproduce fixed
public exon structures across samples and input orders. It does not support
claims about transcript discovery sensitivity, full-genome resource use,
abundance, disease biology, or cross-caller consensus. Those require public raw
reads processed independently through multiple callers.
