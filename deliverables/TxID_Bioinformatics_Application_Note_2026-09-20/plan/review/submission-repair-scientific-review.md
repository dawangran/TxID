# Scientific and bilingual review of the submission repair

Reviewed 23 September 2026. This was a read-only review of the scientific sources and retained results; only this review record was written. No experimental JSON, manuscript paragraph or translation was changed.

## Conclusion

No mandatory correction to the scientific meaning, reported numbers or Chinese translation was found in the reviewed abstract, implementation revision, evaluation chapter or S9 methods. The added case supports a bounded demonstration of exact structural interoperability between two real caller outputs. It does not establish novel-transcript discovery, superiority over reference-alias matching, biological accuracy or abundance accuracy. The manuscript now states these boundaries directly.

This conclusion is not a declaration that the package is ready for immediate journal submission. Author identities, affiliations, Contact and author-supplied declarations remain to be completed; public software/data availability and persistent archive identifiers must be settled, and final four-page journal layout remains unverified. The Availability paragraph accurately describes an accompanying local source archive rather than asserting a verified public repository or inventing a DOI.

## Evidence checked

The review compared `front-matter.md`, `chapters/02_implementation.md`, `chapters/03_evaluation_discussion.md` and `supplementary-real-interop.md` with the retained case at `benchmark-work/submission-repair-2026-09-23/interop-evaluation/`. Its `summary.json` SHA-256 was `1d343f6b3336604490623efcb1ef3a4dddafa66e4a308cee2a1c0341278ca564`. The corresponding English source blocks were compared programmatically with all blocks in the three Chinese translation JSON files, then the translations were read for preservation of scientific scope.

The per-observation evidence table was independently aggregated during review. It contained 126 IsoQuant and 185 StringTie observations, each with a distinct form within its caller. Their form intersection was 110, their union was 201, and their original transcript-name intersection was zero. Direct contig/strand/exon-tuple comparison reproduced the 110 exact cross-caller pairs. Both the false-join denominator and the missed-join denominator are therefore 110; the manuscript's two `0/110` statements do not conceal a zero denominator.

All 110 shared forms were labelled `known` in both independently initialized registries. A separate simple parsing of the supplied reference GTF confirmed that all 220 observations belonging to these shared forms matched reference exon-coordinate tuples, without relying on TxID classification labels for that check. The abstract and evaluation correctly retain the reference-known qualification and explicitly exclude a claim of advantage over reference-alias matching.

The multi-exon layer had 122 and 176 distinct caller chains, 109 shared chain groups and 115 cross-caller chain-equal pairs. Seven of these pairs had different exact forms, distributed across six shared chains. The statement about six chains and seven terminal-variant pairs agrees with the evidence. Group counts and pair counts must remain distinct in tables and any subsequent condensed description.

The four staged registry audits recorded zero coordinate/identity mismatches, zero canonical/full-digest mismatches, zero existing exact-key disappearances or changes and zero existing canonical/digest changes. The final reverse-order comparisons recorded no exact-observation or canonical/digest changes. These checks support the stated incremental exact-key and form-to-chain stability; they do not establish portability of GL1 or FC1 accessions.

## Method and language boundaries retained

The protocol specified the first 100,000 complete FASTQ records and chr22 before caller execution. Whole-genome alignment preceded selection of the 2,765 chr22 primary alignments actually supplied to both callers. The fixed region, one public sample, shared alignment and annotation, and unchanged output GTFs support an interoperability case; they do not make the caller rows independent biological replicates or provide representative discovery-performance estimates.

The annotation contained 4,614 transcript and 27,854 exon records and no gene records. The initial IsoQuant log explicitly requested restarting without `--complete_genedb`; the retained correction and final run record show that change in the fresh `isoquant-inferred` directory. The final caller records use identical BAM and annotation checksums. The delivered caller GTFs match their original successful outputs byte for byte. Both the initial failed attempt and successful correction must remain in the final evidence package.

S9 now states that chain truth includes contig, strand and the complete ordered intron chain, and that a zero denominator yields an undefined rate. It correctly limits independent verification to the separately implemented GTF parser and canonical identity calculation: the assembly fingerprint is supplied from the registry, not independently derived from the FASTA. This is a within-project implementation check, not external replication.

The implementation's “model-defined transcription start and end coordinates” and its Chinese counterpart accurately refer to supplied model boundaries. Neither implies experimental confirmation of transcription initiation or termination. The revised Chinese abstract and evaluation preserve the historical/current version distinction, simulated-output qualification, shared-parser limitation of the historical benchmark, reference-known overlap, two caller rows from one sample, and absence of discovery/abundance claims. No translation reverses a conclusion, drops a denominator or strengthens an independence claim.

## Required final assembly check

At the review snapshot, the evaluation cited Table S6 but no S6 file or result table had yet appeared in the manuscript package. Before delivery, the final supplementary assembly must actually include Table S6, with both exact-form and chain pair denominators and clear distinction between 109 shared chain groups and 115 chain-equal pairs. This is a pending cross-reference/assembly check rather than a request to alter the observed results. The final builder/verifier run must also confirm the updated translations and the English word cap after assembly; this review does not substitute for those checks.

## Optional edits

None is needed to sustain the current scientific claim. Additional discovery benchmarks, wet-lab confirmation or broad literature expansion are not prerequisites for the narrowly stated interoperability result.

## Reviewed source fingerprints

| Source | SHA-256 |
|---|---|
| `front-matter.md` | `272a6e4c93947ac04d752664b6cbd307099536cb8543a79e196a464d4f80d658` |
| `chapters/02_implementation.md` | `82a97ad143a7c741d118d05801e4bb7fdbc6ad46938856e29c0b9d437ef33c97` |
| `chapters/03_evaluation_discussion.md` | `b665b0ff7d2c4fb3b81e8e621bf98907a23dfba52d40b7e4ccaf75af9c1aa1b2` |
| `supplementary-real-interop.md` | `925aa184f04c43cc644f99e2714e6f5a912165dc941d8341017c22ec291922fe` |

## Follow-up: supplementary assembly and evidence verification completed

The pending Table S6 cross-reference above was closed later on 23 September 2026. The package now contains the sixteen-row `tables/table-s6-real-interop.tsv` and the assembled S9 measured-result table and interpretation. The strengthened manuscript verifier checked every Table S6 metric against the retained case and directly recalculated the form/chain joins and every occupancy-matrix cell. Full manuscript verification passed 1,561 checks at this follow-up snapshot, with 1,987 English words including references and caption and 1,267 body words. A later source-manifest addition may change the check count without changing the scientific result.

Three deliberate errors in temporary fixtures were rejected: a changed exon coordinate despite correspondingly refreshed fixture checksums, an incorrect Table S6 overlap count despite a refreshed checksum, and removal of the main-text reference-known qualification. A separate temporary package with no original experiment directory also passed the real-case verifier: all 93 copied evidence files matched their recorded source hashes, and all 19 audited source files were verified directly inside the included source ZIP. The report explicitly distinguishes this package mode from validation against available original files; no unverified source-only hashes remained in that package test.

These follow-up checks resolve the earlier scientific assembly concern. Final root-level deterministic regeneration, archive packaging and submission metadata remain separate responsibilities. The submission-readiness limitations stated above remain in force.
