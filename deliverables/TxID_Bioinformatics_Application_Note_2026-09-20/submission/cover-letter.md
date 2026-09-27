# Cover letter — final scientific text, author fields pending

Complete the bracketed author and availability fields before use. The separate `author-information-needed.md` records the required declarations and access checks. This letter has not been submitted.

---

[SUBMISSION DATE]

Dear Editors,

Please consider our manuscript, “TxID: a reference-aware identity specification for long-read transcript models,” as an Application Note in *Bioinformatics*.

Independent long-read RNA-seq analyses produce transcript catalogs whose identifiers are difficult to reconcile across tools and annotation releases. TxID provides a versioned specification that makes structural identity explicit and independently recomputable. It defines distinct keys for complete splice chains, exact transcript forms and single-exon models within a declared reference sequence context. A registry associates these keys with annotation-dependent classifications, reference aliases and provenance. The contribution is a testable identity specification and its implementation, building on established structural comparison and content-derived naming approaches.

We evaluate this specification using synthetic structures, 538,804 observations from six released ENCODE WTC11 annotations, and a real-read IsoQuant–StringTie interoperability case. Exact identities remained stable across input permutations and GENCODE annotation contexts, while 29,435 classifications changed. Two independently initialized registries joined all 110 coordinate-identical forms shared by the two callers and preserved existing identities during reciprocal incremental imports. A separately implemented parser and identity calculator verified these results. The shared forms were reference-known, so this case demonstrates identity consistency across independent registries without establishing an advantage over reference-alias matching. Comparisons with gffcompare and isoSeQL explicitly distinguish their grouping semantics.

TxID is intended for researchers who maintain or integrate transcript catalogs across samples and independently processed datasets. The Python implementation uses the BSD 3-Clause license and includes a versioned specification, conformance examples, tests and reproducible evaluation materials. The manuscript distinguishes historical TxID 0.1.0 evaluations from the new 0.1.3 case and documents their methodological limits in the text and supplement.

[INSERT VERIFIED PUBLIC SOFTWARE URL, SOFTWARE ARCHIVE IDENTIFIER AND MANUSCRIPT-DATA ARCHIVE IDENTIFIER. Synchronize these with the manuscript Availability and Data availability statements.]

OpenAI Codex assisted with manuscript drafting and revision, translation, code revision, figure-related programming and verification workflows; this assistance is disclosed in the Acknowledgements. [COMPLETE THE TOOL/VERSION/DATE RECORD AND CONFIRM THE NAMED AUTHORS' REVIEW AND RESPONSIBILITY FOR THE RETAINED MATERIAL.]

[CONFIRM ALL AUTHORS' APPROVAL, ORIGINALITY AND EXCLUSIVE SUBMISSION; DISCLOSE RELATED MANUSCRIPTS OR PREPRINTS, FUNDING AND COMPETING INTERESTS AS APPLICABLE.]

Sincerely,

[CORRESPONDING AUTHOR'S FULL NAME]

[DEPARTMENT, INSTITUTION, POSTAL ADDRESS]

[CONTACT EMAIL]
