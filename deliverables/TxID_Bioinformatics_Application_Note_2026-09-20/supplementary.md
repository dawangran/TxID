# Supplementary methods for TxID

## S1 Reference context and canonical structures

TxID v1 defines identity within one reference sequence collection. FASTA sequences are upper-cased and whitespace-normalized. A SHA-256 digest is calculated for each contig. The collection fingerprint hashes records containing the primary contig name, sequence length and sequence digest, sorted by primary name. Consequently, changing the included sequences or their primary names can change the fingerprint even when an assembly label such as GRCh38 is unchanged. Contig aliases resolve input annotations to the registered primary names; they do not make arbitrary independently named FASTA collections identical. The reference files and aliases used in the evaluation are recorded in the benchmark manifests.

External and internal exon coordinates are 1-based closed intervals. Exons are sorted by increasing genomic position, and a transcript must have a consistent contig and strand. For adjacent exons [a,b] and [c,d], the intron interval is [b+1,c−1]. Intron arrays follow transcriptional order, which is reversed on the negative strand. The transcription start site is the minimum exon start on the positive strand and the maximum exon end on the negative strand; the transcription end site is the opposite outer boundary. These sites describe submitted model boundaries, without asserting experimental support for transcription initiation or termination.

Canonical JSON uses lexicographically sorted keys, integer coordinates, UTF-8 encoding, unescaped Unicode and no insignificant whitespace. SC1 objects contain algorithm, assembly fingerprint, primary contig, strand and ordered introns. TF1 adds the exact start and end sites. SE1 instead contains the exact single-exon start and end. The public suffix is the first 24 hexadecimal characters of the SHA-256 digest, giving 96 bits; the full 256-bit digest and canonical object are retained. A public-ID collision with a different digest or canonical object is a hard import error. No suffix is added to repair a collision. Changing these semantics requires a new identity-family version.

For the schematic in Fig. 1a, models M1, M2 and M3 use exon intervals [100,200]/[300,400], [90,200]/[300,420] and [100,201]/[300,400], respectively, on the same positive-strand contig and reference context. M1 and M2 share the intron [201,299] but have different outer boundaries. M3 has intron [202,299]. The donor boundary change is labelled directly as 200 to 201 in the workflow schematic. The labels in the schematic denote equality classes, not measured biological isoforms or literal public digest strings.

## S2 Annotation, provenance and optional similarity

An annotation fingerprint is derived from a normalized, sorted collection of reference transcript structures and their gene and transcript identifiers. The selected annotation context controls classification and aliases but is excluded from exact structural identity. An exact form match uses a reference transcript identifier; where several reference records share a form, selection is deterministic and preferentially retains a matching upstream transcript name. Without an exact match, the current gene-assignment rule uses overlap of aggregate reference gene spans on the same contig and strand. Exactly one candidate supplies the output gene identifier. Multiple candidates remain explicit; a model with no same-strand candidate is classified as a new locus. Antisense-only overlap does not assign a model to the opposite-strand gene.

The SQLite schema separates structural objects, observations, annotation contexts, aliases, locus assignments, fuzzy groupings and import manifests. It stores original gene and transcript names, sample and tool labels, input checksums, options and software version. Rewritten GTF records preserve upstream attributes where their roles remain valid and namespace replaced values. Each input file has an atomic import transaction after staging and validation. Batch execution consists of separate input-file transactions, not one transaction for the entire cohort. Source annotations are not edited in place.

Optional fuzzy grouping uses explicit non-negative splice and end tolerances. Candidate forms must share reference context, primary contig, strand, exon count and family. Corresponding boundaries must satisfy the tolerances for every pair within a cluster. A candidate compatible with more than one established cluster starts a new cluster with ambiguous-bridge status, without merging the existing clusters. In the current implementation, later mutually compatible models can join that new cluster. FC1 membership is recorded in addition to exact IDs. FC1 and GL1 accessions depend on registry-managed entities and must not be used as cross-registry exact structural keys. The present results do not establish the biological validity of fuzzy clusters or gene assignments.

## S3 Public fixed-model inputs

The public evaluation used six released ENCODE WTC11 TALON annotations: three PacBio capped-mRNA outputs from ENCSR309IKK and three Oxford Nanopore direct-RNA outputs from ENCSR392BGY. File accessions and input checksums are provided in Table S1 and its companion TSV. The historical download audit compared transferred bytes and MD5 values with ENCODE metadata and recorded local SHA-256 digests. The 538,804 observations are transcript records across six annotations from one cell line, not 538,804 biological replicates or independently validated transcripts.

Released GTFs also contained reference records. A two-pass extractor selected transcript records with TALON in the source column and retained all exon records linked to those transcript identifiers. It retained 4,248,994 exon records. The historical run then projected the reference and selected sample annotations to transcript/exon records containing gene_id and transcript_id attributes. Coordinates and record order were preserved, while other attributes, including repeated non-identity keys, were counted in the audit and omitted from these benchmark inputs. These projected inputs support evaluation of structural identity, not full-attribute round-trip preservation. The resulting six model sets were evaluated against the GRCh38 sequence collection, using GENCODE v29 and v49 as separate annotation contexts. An NCBI assembly-report audit supplied explicit contig aliases; an incompatible all-regions v49 annotation was excluded from the registry context. Input lists, projection checksums and alias audits are retained with the benchmark artifacts.

This experiment begins after transcript discovery. It assesses identity of released structures and does not compare independent caller reconstructions of the same real reads. The software runs reported here used TxID 0.1.0. The current source checkout identifies itself as 0.1.3; the historical numerical results have not been reassigned to that version.

## S4 Comparator execution and scoring

The fixed-model comparators were gffcompare 0.12.10 and isoSeQL 1.0.1 at commit 25d9366d8b236d3b912e62dcd2c267fe6df31a4c. SQANTI3 6.0.1 was executed on all six model files before isoSeQL ingestion. Fresh sorted and shuffled executions were retained. The unmodified isoSeQL run encountered an ends_counts uniqueness conflict when an experiment contained repeated exact-end observations with the same read-count value. An audited compatibility change replaced only that table's insertion with INSERT OR IGNORE. Structural lookup and PBID insertion were unchanged, and all 538,804 source transcript observations were retained. These results concern structural mapping, not the validity of abundance aggregation after that modification.

The scoring script defined the exact target by direct equality of contig, strand and all ordered exon intervals on the selected assembly, rather than by comparing TxID digests. It used the shared TxID annotation parser, so this is a coordinate-level consistency check, not validation with an independent parser. A separate splice-chain target used contig, strand and the complete ordered intron chain; single-exon records retained their exact intervals. Let T(i) denote the target group of observation i and P(i) its reported output group. False-merge rate is the number of unordered pairs with P(i)=P(j) and T(i)≠T(j), divided by all pairs sharing P. False-split rate is the number of pairs with T(i)=T(j) and P(i)≠P(j), divided by all pairs sharing T. Each rate uses the eligible emitted observations; missing observations are reported separately. Zero pair counts are descriptive consistency results and do not estimate an upper population error bound.

Table S2 keeps the identity layers and targets explicit. The isoSeQL exact-end layer produced 93 separated pairs out of 100,680 equal-target pairs, or 0.0924%, and no false-merged pair in this run. gffcompare omitted two observations; 55,494 of its output clusters joined only terminal variants of a common splice chain. Scoring those groups against full exon-interval equality measures a difference in equivalence rules as well as any other grouping differences. It does not by itself identify a software defect. The separate splice-chain comparison is therefore retained. isoSeQL common-junction results likewise use the complete-chain target rather than being treated as full-form equality. No single overall ranking is inferred from these layers.

The retained reports also contain percentile intervals from 2,000 structural-unit bootstrap replicates. The false-merge calculation resampled output groups with seed 20260731; the false-split calculation resampled target structures with seed 20260732. These intervals describe resampling of the evaluated structural groups, not independent biological sampling or uncertainty in caller discovery. The short manuscript and its tables use raw pair counts, denominators and descriptive rates. We do not attach inferential significance tests to the deterministic equality checks.

## S5 Order, annotation and incremental audits

Order comparisons distinguish numeric or public label changes from changes in pairwise co-membership (Table S3). Under the recorded permutation, TxID exact-form and splice-chain identifiers did not change. The isoSeQL exact-end and common-junction partitions also did not change, although their database labels changed for all 538,804 observations. This is consistent grouping with different assigned labels. gffcompare changed 4,438 pair relations among its 538,802 shared observations. The TxID locus layer changed 137,658 output gene labels and 418,224 pair relations; the complete sorted and shuffled catalogs therefore had different bytes. Exact-ID invariance must not be extended to those locus assignments or to the entire exported catalog.

Reclassification between GENCODE v29 and v49 paired all 538,804 observations and found no exact form or splice-chain identifier changes. Classification changed for 29,435 observations (5.463%). The v49 registry persisted five input files. A storage quota failure rolled back the sixth import, and a later retry encountered an I/O failure. That sixth sample, ENCFF973ZQF, was processed with the same reference index and classification functions in a read-only run without observation persistence or locus allocation. Thus the result is a six-input structure and classification comparison, not a completed six-input persistent v49 registry. The unrecoverable generated database was removed in the historical recovery process; checked mappings and compact audits were retained.

Six incremental catalog snapshots showed no disappearance of exact identifiers, changes to existing form-to-splice-chain associations or aggregate classification regressions. The association check compared the stored splice-chain identifier for each existing form; it did not revalidate every canonical-object field at every snapshot. The first four catalogs reconstructed from checked component mappings matched direct SQLite exports byte for byte. The same deterministic aggregation generated snapshots five and six, with the final reconstruction matching the final sorted catalog. This checks consistency of accumulated structural records and exports; it is not evidence that every snapshot survived as a persisted database. Historical recovery audits document an interrupted sorted run, journal recovery and validated resumption.

The sample-by-group matrices in Table S4 use structural occupancy. An occupied cell indicates at least one model observation for that sample and identity group. Matrix sparsity is one minus occupied cells divided by the product of sample rows and identity columns. The matrices are not read-count or TPM matrices. The exact-form matrix has six rows, 467,819 columns and 538,804 occupied cells (80.804% sparsity). Different grouping rules change these dimensions; fewer columns alone do not indicate a more biologically accurate transcriptome. Read-count tables generated elsewhere can be joined through the mapping outputs, but TxID itself does not perform quantification.

## S6 Synthetic and integration checks

The retained structural simulation used seed 20260730, 120 latent transcript templates, six simulated samples and four explicitly labelled caller-like serializers, with a configured 5% boundary-perturbation rate. These labels do not represent executions of the named transcript callers. The run produced 2,511 model observations, including 128 perturbed observations, and 190 distinct emitted exact forms. Relative to emitted structure, there were zero false splits and merges across 22,720 eligible same-target pairs. Relative to the latent templates, 70 of 120 templates were fragmented and 2,432 of 25,152 same-template pairs were separated. This is the expected cost of retaining coordinate differences rather than correcting them. The 2,511-to-190 matrix-column change uses an intentionally unique local-name baseline and is not an estimate of identifier fragmentation in real datasets.

A separate Tier-0 integration test used sequence-backed synthetic reads and executed minimap2, IsoQuant, FLAIR, StringTie, TALON, gffcompare, SQANTI3 and isoSeQL. Its eight truth forms, seven caller annotations and 36 model observations establish that the workflow stages can run together. They do not provide estimates of caller accuracy in biological samples. TALON was evaluated through its supported shared-database read workflow; its read-pair denominators are not pooled with fixed-model pair counts. Timing and memory from these unmatched stages are not used for a performance ranking.

## S7 Reproducibility files

The source repository contains the versioned identity specification, schema definitions, conformance vectors, tests and workflow scripts. The supplementary TSVs and Fig. 1 source data are derived directly from the compact retained JSON reports in `benchmarks/results/encode-gtf-full/`. Within the manuscript package, `reproducibility/source-manifest.json` records source checksums, `reproducibility/build_artifacts.py` generates the tables and figure, and `figures/source-data.json` contains the plot inputs. The historical full benchmark has not been rerun for this draft. Large source annotations, databases and intermediate outputs are not bundled with this manuscript; public input accessions and checked download manifests identify the source models.

The illustrative structures and empirical numerical values are distinguished in the figure data manifest. Figure and table generation can be repeated from the retained reports without rerunning the full benchmark.

Section S9 and Table S6 report a separate real-read interoperability case executed with TxID 0.1.3. The accompanying local source archive contains the tested software snapshot; the scientific materials archive includes the new caller GTFs, chr22 reference annotation, per-observation coordinate evidence, caller-presence matrices, execution records and SHA-256 manifests. Original failure records and the successful IsoQuant execution correction are retained together. Raw reads, the complete reference FASTA, BAMs and SQLite databases are excluded from this compact package. The public accession and reference checksum identify the files needed to repeat the new computation. Public permanent deposition of the submission version remains distinct from constructing and verifying these local archives.

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

## S8 Identity semantics and related workflows

Table S5 compares the objects and matching rules described by each tool, rather than ranking performance. Isosceles provides a direct precedent for stable transcript-path hashes (Kabza et al., 2024). Its same-genome-build scope does not establish equivalence to TxID's fingerprint of a selected FASTA sequence collection and its primary contig names. Consequently, the present comparison supports TxID's role as a separate versioned exact-identity layer, without claiming that transcript hashing, incremental aggregation or intron-chain/end separation originated with TxID.

## Table S5. Identity semantics and naming contexts

| Tool | Input object | Equivalence or end treatment | Naming context | Evidence |
| --- | --- | --- | --- | --- |
| TxID | GTF/GFF3 models and reference FASTA; annotation for classification | SC: exact intron chain; TF: chain and exact ends; SE: exact single-exon interval | Versioned hashes include the reference collection fingerprint; annotation classification is separate | TxID v1 specification; Supplementary S1–S2 |
| Isosceles | Aligned long reads and reference annotation for discovery and quantification | Annotation preparation merges matching intron structures and matching start/end bins; default bin width 50 bp | Stable hashes for novel transcript paths, described for one genome build; serialization and FASTA-fingerprint policy not established from the cited paper | Kabza et al. (2024), pp. 2, 7 |
| isoSeQL | SQANTI3 classification and genePred files, with sample/experiment metadata | Intron-chain groups linked to separate start/end coordinate combinations | Database isoform numbers and linked end coordinates; additional samples can enter the same SQLite database | Liu and Chun (2026), Sections 2.1 and 3.2 |
| TALON | SAM read alignments and a database initialized from reference GTF | Read annotation uses splice junctions and configurable end distances; documented defaults: 5′ 500 bp, 3′ 300 bp | Database-issued transcript IDs under selected build and annotation names; later datasets reuse the database | Wyman et al. (2020); official TALON documentation |
| GffCompare | Transcript models in GTF/GFF files; optional reference annotation | Multi-exon tracking matches introns and permits terminal-length differences | Tracking output assigns structure and locus labels across supplied files | Pertea and Pertea (2020); official GffCompare tracking documentation |

The rows describe different operations: TALON's read annotation, Isosceles' annotation preparation and transcript-path naming, and GffCompare's multi-file tracking are not interchangeable with TxID exact-form equality. GffCompare's accuracy-estimation options should not be substituted for its tracking definition. Unknown implementation details marked “not established” are evidence limits, not claims that a feature is absent. The cited material does not establish a common versioned serialization or FASTA-collection fingerprint contract across these tools. TxID GL1 and FC1 accessions remain registry-managed and are outside the exact-key comparison.

The Isosceles statements were checked in the [publisher PDF](https://www.nature.com/articles/s41467-024-51584-3.pdf), with stable hashes described on page 2 and annotation end bins on page 7. The operational descriptions for [TALON](https://github.com/mortazavilab/TALON#initializing-a-talon-database) and [GffCompare tracking](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) were checked against official documentation on 23 September 2026. Defaults are reported as documented there, rather than retrospectively assigned to every historical release. The companion `tables/table-s5-identity-semantics.tsv` contains the same comparison.

## S9 Real-read interoperability and a separately implemented identity check

The additional case uses TxID 0.1.3 and public WTC11 PacBio CCS reads from ENCODE accession ENCFF105WIJ (ENCSR309IKK). Its selection protocol was recorded before caller execution. The complete retained FASTQ matched SHA-256 `a335dc75d3356e6e16dc1097f856a050b1f3d978732935ab5b6b0f14848a3427`. The first 100,000 complete records in file order contained 112,715,597 sequenced bases. This bounded selection was not chosen by mapping success, caller agreement or TxID outcomes. It is a subset of one public sample, not 100,000 biological replicates.

The subset was aligned to the complete GRCh38 no-alt analysis sequence collection, GCA_000001405.15, with minimap2 2.31-r1302 (`-ax splice:hq -uf --secondary=no --MD -t 8`) and sorted with samtools 1.23.1. Chromosome chr22 alignments were then selected, excluding unmapped, secondary and supplementary records with mask 2308; original genomic coordinates were retained. The resulting 2,765 primary alignments, rather than all 100,000 reads, were supplied to each caller. Caller execution used the chr22 transcript/exon records from the retained GENCODE v29 identity-view annotation (32,468 records). This reference projection retains gene and transcript identifiers, with other attributes omitted as described in S3. Both callers received exactly this same annotation and alignment file. The reference FASTA and complete input file checksums are recorded in the selection manifest.

IsoQuant 3.13.0 ran in `pacbio_ccs` mode, and StringTie 3.0.3 ran with `-L -G`; both used four threads. IsoQuant initially rejected `--complete_genedb` because the projected annotation lacks gene records. Its documented inference of missing annotation meta-features was therefore enabled by omitting that option in a fresh working directory; the failed attempt and execution correction are retained. This adjustment did not change the selected reads, coordinates, reference or caller matching thresholds. The experiment runner records full argument lists, versions and output checksums. An IsoQuant entry wrapper redirects mutable configuration and cache paths to the new working directory; it does not change reconstruction code. Caller GTFs were passed to TxID without coordinate or attribute projection. Neither caller was tuned to increase the exact-form or splice-chain intersection. The two outputs are alternative reconstructions of the same selected sample; their matrix rows denote caller outputs, not biological replicates. Runtime records provide provenance and are not a ranking of the callers.

Two empty registries were initialized independently using the same full FASTA and chr22 annotation. One first received IsoQuant models and the other StringTie models. Their initial mapping tables were joined on exact form identifiers and, separately, on multi-exon SC1 identifiers. Each registry subsequently received the other caller's unchanged GTF. The audit checks existing exact keys and form-to-splice-chain associations after this increment, and compares exact structures across the two final import orders. It also validates stored canonical objects and full digests. GL1 locus accessions and FC1 fuzzy clusters are outside these cross-registry comparisons.

The oracle in `benchmarks/independent_identity.py` uses only the Python standard library and imports no TxID modules. It independently parses the two GTFs, orders their exon intervals and calculates SC1/TF1/SE1 canonical JSON, full SHA-256 digests and public IDs from the published specification. It first reproduces all three existing conformance vectors, comprising five non-null identities. Assembly identity is supplied from the selected registry context, with contig lengths checked; the oracle does not independently derive the FASTA fingerprint. This is a separately implemented parser and identity calculator within the project, not an external replication or a biological truth set.

For the real-case join audit, direct contig, strand and exon-interval tuples define exact-form equality; contig, strand and complete ordered introns define multi-exon chain equality within the shared reference context. A false join is a cross-caller pair sharing a reported ID but not the corresponding coordinate target. A missed join is a cross-caller coordinate-equal pair receiving different IDs. Both denominators and missing observations are reported; a zero denominator gives an undefined rate. Chain-only matches retain their distinct exact forms. Binary presence matrices and per-observation evidence tables accompany the results, so the consequences of form versus chain equality can be inspected without inferring abundance. This case addresses interoperability on one real read subset and two callers; it does not measure transcript discovery accuracy, differential expression or generalization to other tissues.

### Table S6. Real-read interoperation and structural join audit

| Quantity | IsoQuant | StringTie | Cross-caller result |
| --- | ---: | ---: | --- |
| Transcript observations / distinct exact forms | 126 / 126 | 185 / 185 | 110 shared forms; 201-form union |
| Distinct multi-exon splice chains | 122 | 176 | 109 shared chains; 189-chain union |
| Forms without an exact counterpart in the other output | 16 | 75 | Unequal coordinates, not identity false splits |
| Exact-form join | — | — | 110/110 equal-coordinate pairs recovered; 0/110 false merges and 0/110 false splits |
| Multi-exon splice-chain join | — | — | 115/115 equal-chain pairs recovered; 0/115 false merges and 0/115 false splits |
| Original-ID string join (descriptive) | — | — | 0 joined pairs; 110/110 equal-form pairs missed |
| Shared chains containing distinct forms | — | — | 6 chains; 7 same-chain/different-end cross-caller pairs |
| Changes after incremental addition | 0/126 existing observation keys | 0/185 existing observation keys | Canonical objects and full digests also unchanged |
| Changes between final registries after opposite import orders | — | — | 0/311 observation keys; 0/8,684 structural objects |

Pair counts refer to cross-caller observation pairs, whereas shared-form and shared-chain counts refer to distinct structural groups. False-merge denominators are joined pairs and false-split denominators are coordinate-equal pairs under the stated identity layer. The original-ID join emits no pairs, so its false-merge rate is undefined. The chain analysis excludes single-exon forms. Its 115 matched observation pairs belong to 109 shared chains; seven pairs have the same complete intron chain but different transcript ends, distributed across six chains. These end differences remain distinct exact forms. The two binary occupancy matrices have caller outputs as rows: the form matrix contains 201 columns and 311 occupied cells, and the multi-exon chain matrix contains 189 columns and 298 occupied cells. They show exact structural reconciliation of two outputs from one sample and do not constitute a sample-by-transcript expression matrix or an abundance analysis.

The classification audit is essential to interpreting the matches. All 110 shared exact forms (108 TF1 and two SE1) were classified as known against the selected annotation in both outputs. Their preserved reference transcript identifiers also agree in all 110 pairs. The zero intersection of original caller-ID strings is therefore only a descriptive naming baseline; it does not show an advantage over joining shared reference aliases. This selected case demonstrates that independently assigned structural keys reproduce coordinate equality in actual caller outputs, but it provides no example of a shared novel exact form. The 16 IsoQuant-only forms comprise eight known, six novel-in-known-gene and two ambiguous-gene classifications; the 75 StringTie-only forms comprise 69 known, four novel-in-known-gene and two ambiguous-gene classifications. Independent-reader coordinate tuples confirm that none has an exact counterpart in the other output. These are differences between caller model structures rather than erroneous splitting of identical structures, and neither their biological validity nor the relative accuracy of the callers was evaluated.

The compact evidence package at `reproducibility/real-interop/` retains unchanged copies of the caller GTFs, initial and incremental mappings, catalogs, matrices, per-observation canonical objects and digests, execution manifests and success/failure logs. `evidence-manifest.json` records and verifies the SHA-256 of every copied file against its original; `supplement-classification-audit.json` adds the classification analysis without changing the original evaluation summary. Table S6 (`tables/table-s6-real-interop.tsv`) records the metric units and denominators. Registry-managed locus and fuzzy accessions were excluded from the cross-registry comparison. The test addresses implementation conformance and structural interoperability; it does not establish independent biological truth or replace the separately acknowledged shared-database TALON comparison.
