# TxID: a reference-aware identity specification for long-read transcript models

[AUTHOR NAMES AND AFFILIATION INDICES]

[AFFILIATIONS: DEPARTMENT, INSTITUTION, CITY, POSTCODE, COUNTRY]

Corresponding author: [CORRESPONDING AUTHOR NAME]

## Abstract

### Summary

TxID provides a versioned, reference-aware identity specification and registry for GTF/GFF3 transcript models, enabling independent analyses to recompute common structural keys for exact splice chains, transcript forms and single-exon models while recording annotation-dependent classifications and provenance separately. Evaluation across 538,804 ENCODE WTC11 observations showed stable exact identifiers despite input permutation and 29,435 annotation-classification changes; a real-read IsoQuant–StringTie case joined 110 shared reference-known forms across independently initialized registries, with a separately implemented parser and identity calculator confirming the joins.

### Availability and Implementation

TxID is implemented in Python 3.10 or later, has no third-party runtime dependencies and is distributed under the BSD 3-Clause license. Source code, the versioned specification, tests and installation instructions are provided at [https://github.com/dawangran/TxID](https://github.com/dawangran/TxID).

### Contact

[CONTACT EMAIL]

### Supplementary Information

Supplementary methods, input accessions, comparison tables and figure source data accompany this manuscript.

## 1 Introduction

Joining independently generated long-read RNA-seq catalogs requires a common definition of which transcript models represent the same structure. In the SQANTI3 workflow addressed by isoSeQL, sample-specific isoform identifiers complicate comparison, while jointly recollapsing an expanding collection can become costly (Liu and Chun, 2026). Structural equivalence also depends on the comparison: two models may have identical introns yet differ at their transcript ends.

Existing approaches support reconciliation at different stages. GffCompare tracks models across annotation files using structural comparisons that can group terminal variants (Pertea and Pertea, 2020). TALON tracks transcripts and abundance within a shared database (Wyman et al., 2020), whereas isoSeQL incrementally aggregates SQANTI3 outputs, linking intron chains, transcript ends, metadata and abundances (Liu and Chun, 2026). Isosceles assigns stable hashes to novel transcript paths for cross-study matching on the same genome build; its annotation preparation merges transcripts with matching intron structures and start/end bins, with a default bin size of 50 bp (Kabza et al., 2024).

Content-derived naming also has precedents in gene annotation through UniqTag (Jackman et al., 2015) and genomic variation through the GA4GH Variation Representation Specification (Wagner et al., 2021). TxID provides a versioned, reference-aware identity specification for supplied GTF/GFF3 transcript models, with an accompanying registry. It specifies exact splice-chain, transcript-form and single-exon identities, the reference sequence-collection fingerprint and canonical serialization. This explicit contract allows separate implementations to recompute and verify structural keys without shared assignment history. Annotation-dependent classifications, registry-managed locus and fuzzy-cluster assignments, and observation provenance remain separate from exact identity.

## 2 Implementation

TxID implements a versioned identity specification for GTF/GFF3 models against a declared reference FASTA and annotation. The reference fingerprint hashes sorted primary contig names, lengths and SHA-256 digests of normalized sequences. Input contigs resolve through primary names or an explicit alias table; renaming FASTA contigs changes the context. Models use 1-based closed coordinates. Exons are sorted genomically; adjacent exons [a,b] and [c,d] define intron [b+1,c−1], with introns ordered by transcriptional direction. Validation rejects unknown contigs, out-of-range coordinates and overlapping exons. Exact identity is confined to this sequence collection; cross-assembly equivalence is not established.

Canonical JSON objects contain the algorithm family, reference fingerprint, primary contig and strand. For multi-exon models, `txid:SC1` encodes the complete intron chain; `txid:TF1` adds the exact, strand-aware outer boundaries defining model transcription start and end. Identical chains with different ends therefore share SC1 but have distinct TF1 identifiers (Fig. 1a). Single-exon models receive `txid:SE1` identifiers from exact start and end coordinates. The specification fixes lexicographically sorted JSON keys, integer coordinates, UTF-8 encoding, unescaped Unicode and no insignificant whitespace, allowing the canonical bytes and identifiers to be recomputed independently. SHA-256 hashes these bytes; public identifiers use its first 24 hexadecimal characters. The registry stores the full digest and canonical object and rejects conflicting public identifiers. Samples, upstream names, annotation releases and import order are excluded from the hash inputs. Changing these rules requires a new identity-family version.

Annotations receive separate normalized fingerprints and supply classifications and reference aliases. An exact reference-form match preserves the selected reference `gene_id` and `transcript_id` in rewritten annotations, with TxIDs recorded separately. Otherwise, overlap with exactly one same-strand reference gene span assigns that gene identifier and a TxID transcript identifier. Multiple overlapping genes remain explicit candidates; ambiguous assignments and models without same-strand gene overlap receive persistent `txid:GL1` locus accessions. Gene-span overlap records an assignment rather than biological gene membership. Optional fuzzy grouping uses recorded splice and end tolerances, requires mutual compatibility within each cluster and marks ambiguous bridges. Its `txid:FC1` accessions preserve exact identities. FC1 and GL1 are registry-managed accessions without the portability guarantees of exact structural identifiers.

A versioned SQLite registry separates structural objects from observations, annotation contexts and provenance. Each input is staged and validated before an atomic transaction; checksums and manifests support idempotent retries. Provenance records sample, upstream tool, original identifiers and attributes, software version and options. Outputs from a fixed registry are deterministically sorted and comprise rewritten GTFs, per-sample mapping tables and a cohort catalog. Preserved reference names and TxID transcript names coexist, while `txid_form` provides a common exact key for joining external sample-by-transcript matrices. TxID operates on supplied models and measurements and does not estimate abundance.

## 3 Evaluation and discussion

Evaluation tested conformance to the identity specification and portability of its structural keys. Historical experiments used TxID 0.1.0. A simulation generated 2,511 observations from 120 latent templates through boundary perturbations and four caller-style serializers. TxID recovered 190 emitted forms without false-merge or false-split pairs against exact emitted structure. Perturbations nevertheless fragmented 70 latent templates into multiple exact forms, illustrating how faithful structural identification preserves differences introduced upstream.

For fixed models, six released ENCODE WTC11 TALON annotations supplied 538,804 observations from three PacBio capped-mRNA and three ONT direct-RNA outputs on GRCh38 (ENCODE Project Consortium et al., 2020). Equality of contig, strand and exon intervals defined the structural target. The scoring script reused the TxID parser, making this a coordinate-consistency check. TxID assigned 467,819 exact forms, with no false splits among 100,680 same-target pairs or false merges among the same number of same-ID pairs.

Comparisons used gffcompare 0.12.10 and isoSeQL 1.0.1 after SQANTI3 6.0.1 characterization (Pardo-Palacios et al., 2024). isoSeQL exact ends yielded 467,908 groups, no false-merge pairs and 93 false-split pairs. Its disclosed compatibility patch used `INSERT OR IGNORE` for duplicate count records, retaining source names and structural queries. gffcompare produced 379,407 tracking groups and omitted two observations. Its end-tolerant rule combined terminal variants in 55,494 groups; disagreement with exact-form equality is therefore not an overall performance ranking. Table S2 separates form and chain targets and pair denominators.

These equivalence rules determine matrix columns (Fig. 1c). Six-row binary occupancy matrices contained 467,819 TxID exact-form columns, 467,908 isoSeQL exact-end columns and 379,407 gffcompare tracking columns. Their entries record model presence. TxID matrix sparsity was 80.804%; column counts and sparsity describe structural organization rather than quantification accuracy.

Identity stability was tested across input order and annotation contexts. The recorded permutation preserved TxID exact identifiers and partitions. isoSeQL preserved its exact-end partition, although all numeric labels changed (Table S3). TxID locus accessions remained stateful, and full catalogs were not byte-identical. Between GENCODE v29 and v49 (Mudge et al., 2025), exact IDs remained unchanged while 29,435 classifications changed (5.463%; Fig. 1b). Following storage failure, the sixth v49 input was classified without persisting observations. Incremental snapshots, partly reconstructed from mappings, showed no lost exact keys or changed form-to-splice-chain associations; the final reconstruction matched the final catalog byte for byte.

A real-caller case tested TxID 0.1.3. Whole-genome alignment of the first 100,000 reads of ENCFF105WIJ yielded 2,765 chr22 primary alignments. IsoQuant 3.13.0 and StringTie 3.0.3 produced 126 and 185 forms, respectively (Prjibelski et al., 2023; Kovaka et al., 2019). Independently initialized registries joined 110 shared forms into a two-row, 201-column caller-occupancy matrix despite disjoint original transcript names. All 110 were reference-known and shared reference aliases; this case demonstrates neither shared novel forms nor an advantage over reference-alias matching. A separately implemented GTF reader and identity calculator confirmed zero false or missed exact joins (0/110 pairs each). Six shared chains contained seven cross-caller pairs with different ends, retaining distinct form IDs. Reciprocal incremental imports preserved existing keys, canonical objects and form-to-chain associations. Table S6 reports chain comparisons and denominators; matrix rows represent caller outputs of one sample.

The evidence covers one cell line, fixed TALON models and a bounded two-caller case. Exact identities retain caller boundary errors; discovery accuracy, expression estimates, gene-span assignment and fuzzy grouping were not biologically validated. The comparisons did not evaluate the shared-database TALON workflow. Within this scope, TxID makes structural equivalence explicit and independently checkable, providing portable keys under a shared reference context while retaining annotation-dependent interpretation and provenance.

## Acknowledgements

### Funding

[AUTHOR FUNDING STATEMENT]

### Author contributions

[AUTHOR CONTRIBUTIONS]

### Use of AI tools

OpenAI Codex assisted with drafting, translation, code revision and verification workflows. The installed Codex CLI version recorded during final preparation was 0.142.2. [AUTHOR REVIEW STATEMENT]

## Conflict of interest

[AUTHOR COMPETING-INTEREST STATEMENT]

## Data availability

Public ENCODE input accessions are listed in Table S1; the real-read case uses ENCFF105WIJ. Numerical source data, input manifests, analysis scripts and unmodified caller annotations are provided in the [versioned reproducibility materials](https://github.com/dawangran/TxID/tree/35c124be9b0e76cd4d71b39c3ec394d035e5be8d/deliverables/TxID_Bioinformatics_Application_Note_2026-09-20).

## Figure 1

![Reference-aware identity workflow, annotation-classification transitions and matrix column counts.](figures/figure1.png)

Figure 1. Structural identity, annotation context and catalog integration. (a) GTF/GFF3 models receive exact identifiers within a selected FASTA context; annotation-dependent classifications and provenance are recorded separately. Schematic models share a reference, contig and positive strand: M1 and M2 share a splice chain (SC) but differ in transcript form (TF); M3 shifts the M1 donor boundary from 200 to 201. Coordinates are 1-based closed; labels denote equality classes. (b) Classifications changed for 29,435 of 538,804 paired observations between GENCODE v29 and v49 (5.463%), with no exact-form identifier changes. Ribbon widths encode transition counts and colours indicate the v29 class; unchanged classifications are omitted. Node counts include changed observations only. “Novel in gene” denotes a novel transcript in a known gene. The sixth v49 input was classified without database persistence. (c) Column counts for six-row binary occupancy matrices from sorted-input runs under TxID exact form, isoSeQL exact ends, gffcompare tracking and isoSeQL common junction. Occupancy records model presence; column counts reflect the grouping relation and do not measure abundance or biological accuracy. Panels b and c use fixed models; the isoSeQL count-insertion compatibility patch is described in Supplementary Section S4. Source values accompany the figure.

Alt text: Reference-aware transcript identity workflow with separate annotation classification. A flow diagram shows changed classifications between annotation releases, while exact IDs remain unchanged. Four bars compare structural matrix column counts.

## References

ENCODE Project Consortium, Moore JE, Purcaro MJ et al. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 2020;583:699–710. [doi:10.1038/s41586-020-2493-4](https://doi.org/10.1038/s41586-020-2493-4).

Jackman SD, Bohlmann J, Birol I. UniqTag: content-derived unique and stable identifiers for gene annotation. PLoS One 2015;10:e0128026. [doi:10.1371/journal.pone.0128026](https://doi.org/10.1371/journal.pone.0128026).

Kabza M, Ritter A, Byrne A et al. Accurate long-read transcript discovery and quantification at single-cell, pseudo-bulk and bulk resolution with Isosceles. Nat Commun 2024;15:7316. [doi:10.1038/s41467-024-51584-3](https://doi.org/10.1038/s41467-024-51584-3).

Kovaka S, Zimin AV, Pertea GM et al. Transcriptome assembly from long-read RNA-seq alignments with StringTie2. Genome Biol 2019;20:278. [doi:10.1186/s13059-019-1910-1](https://doi.org/10.1186/s13059-019-1910-1).

Liu CS, Chun J. isoSeQL: comparing long-read isoforms across multiple datasets. Bioinformatics 2026;42:btaf680. [doi:10.1093/bioinformatics/btaf680](https://doi.org/10.1093/bioinformatics/btaf680).

Mudge JM, Carbonell-Sala S, Diekhans M et al. GENCODE 2025: reference gene annotation for human and mouse. Nucleic Acids Res 2025;53:D966–D975. [doi:10.1093/nar/gkae1078](https://doi.org/10.1093/nar/gkae1078).

Pardo-Palacios FJ, Arzalluz-Luque A, Kondratova L et al. SQANTI3: curation of long-read transcriptomes for accurate identification of known and novel isoforms. Nat Methods 2024;21:793–797. [doi:10.1038/s41592-024-02229-2](https://doi.org/10.1038/s41592-024-02229-2).

Pertea G, Pertea M. GFF Utilities: GffRead and GffCompare [version 2; peer review: 3 approved]. F1000Research 2020;9:304. [doi:10.12688/f1000research.23297.2](https://doi.org/10.12688/f1000research.23297.2).

Prjibelski AD, Mikheenko A, Joglekar A et al. Accurate isoform discovery with IsoQuant using long reads. Nat Biotechnol 2023;41:915–918. [doi:10.1038/s41587-022-01565-y](https://doi.org/10.1038/s41587-022-01565-y).

Wagner AH, Babb L, Alterovitz G et al. The GA4GH Variation Representation Specification: a computational framework for variation representation and federated identification. Cell Genom 2021;1:100027. [doi:10.1016/j.xgen.2021.100027](https://doi.org/10.1016/j.xgen.2021.100027).

Wyman D, Balderrama-Gutierrez G, Reese F et al. A technology-agnostic long-read analysis pipeline for transcriptome discovery and quantification. bioRxiv 2020;672931, version 2, posted 24 March 2020. [doi:10.1101/672931](https://doi.org/10.1101/672931).
