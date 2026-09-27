# Supplementary tables

## Table S1. Public fixed-model inputs

| Accession | Platform | Replicate | Transcripts | Exons |
| --- | --- | --- | --- | --- |
| ENCFF694HZY | pacbio | 1 | 62906 | 484193 |
| ENCFF978SRK | pacbio | 2 | 56211 | 400682 |
| ENCFF973ZQF | pacbio | 3 | 77677 | 599605 |
| ENCFF527AJC | ont | 1 | 94130 | 695873 |
| ENCFF615FZJ | ont | 2 | 189280 | 1699489 |
| ENCFF798LYI | ont | 3 | 58600 | 369152 |

All inputs are WTC11 TALON annotations. PacBio replicates belong to ENCSR309IKK; ONT replicates to ENCSR392BGY. Full download URLs and SHA-256 checksums are in table-s1-inputs.tsv.

## Table S2. Structural pair comparisons

| Layer / target | Emitted / input | Groups | False-merge pairs / eligible pairs | False-split pairs / eligible pairs |
| --- | --- | --- | --- | --- |
| TxID exact form / Exact form | 538804 / 538804 | 467819 | 0 / 100680 | 0 / 100680 |
| gffcompare tracking / Exact form | 538802 / 538804 | 379407 | 171218 / 271895 | 0 / 100677 |
| gffcompare tracking / Complete splice chain | 538802 / 538804 | 379407 | 20458 / 271895 | 0 / 251437 |
| isoSeQL exact ends / Exact form | 538804 / 538804 | 467908 | 0 / 100587 | 93 / 100680 |
| isoSeQL common junction / Complete splice chain | 538804 / 538804 | 364032 | 90319 / 335842 | 5917 / 251440 |

Rates and separate denominators are included in table-s2-structural-comparison.tsv. Targets are coordinate-based comparisons using the shared parser. Missing observations are excluded from eligible pairs and explicitly reported. The isoSeQL comparison used the compatibility patch described in Supplementary Section S4. Different equivalence relations are not a general performance ranking.

## Table S3. Input-order effects

| Layer | Shared observations | Changed observation labels | Changed pair relations |
| --- | --- | --- | --- |
| TxID exact form | 538804 | 0 | 0 |
| TxID gene/locus | 538804 | 137658 | 418224 |
| gffcompare tracking | 538802 | 538768 | 4438 |
| isoSeQL exact ends | 538804 | 538804 | 0 |
| isoSeQL common junction | 538804 | 538804 | 0 |

Changed observation labels and changed pair relations have different units. The TxID locus layer is separate from exact identity. isoSeQL label changes coexist with unchanged partitions.

## Table S4. Sample-by-group dimensions

| Layer | Rows | Columns | Occupied cells | Sparsity (%) |
| --- | --- | --- | --- | --- |
| TxID exact form | 6 | 467819 | 538804 | 80.804 |
| gffcompare tracking | 6 | 379407 | 538802 | 76.331 |
| isoSeQL exact ends | 6 | 467908 | 538804 | 80.808 |
| isoSeQL common junction | 6 | 364032 | 530082 | 75.731 |

Each matrix has six sample rows. Occupancy represents model observations, not RNA abundance. These layers define different kinds of columns.
