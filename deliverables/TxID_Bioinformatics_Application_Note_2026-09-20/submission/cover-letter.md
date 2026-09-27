# Cover letter

27 September 2026

Dear Editors,

Please consider our manuscript, “TxID: a reference-aware identity specification for long-read transcript models,” as an Application Note in *Bioinformatics*.

Independent long-read RNA-seq analyses produce transcript catalogs whose identifiers are difficult to reconcile across tools and annotation releases. TxID provides a versioned specification that makes structural identity explicit and independently recomputable. It defines distinct keys for complete splice chains, exact transcript forms and single-exon models within a declared reference sequence context. A registry associates these keys with annotation-dependent classifications, reference aliases and provenance. The contribution is a testable identity specification and its implementation, building on established structural comparison and content-derived naming approaches.

We evaluate this specification using synthetic structures, 538,804 observations from six released ENCODE WTC11 annotations, and a real-read IsoQuant–StringTie interoperability case. Exact identities remained stable across input permutations and GENCODE annotation contexts, while 29,435 classifications changed. Two independently initialized registries joined all 110 coordinate-identical forms shared by the two callers and preserved existing identities during reciprocal incremental imports. A separately implemented parser and identity calculator verified these results. The shared forms were reference-known, so this case demonstrates identity consistency across independent registries without establishing an advantage over reference-alias matching. Comparisons with gffcompare and isoSeQL explicitly distinguish their grouping semantics.

TxID is intended for researchers who maintain or integrate transcript catalogs across samples and independently processed datasets. The Python implementation uses the BSD 3-Clause license and includes a versioned specification, conformance examples, tests and reproducible evaluation materials. The manuscript distinguishes historical TxID 0.1.0 evaluations from the new 0.1.3 case and documents their methodological limits in the text and supplement.

Source code, the versioned specification, tests and installation instructions are provided at [https://github.com/dawangran/TxID](https://github.com/dawangran/TxID). The Data availability statement links the experimental evidence at the fixed repository revision used for the final scientific text.

OpenAI Codex assisted with drafting and revision, translation, code revision, figure-related programming and verification workflows, as disclosed in the Acknowledgements. The installed Codex CLI version recorded during final preparation was 0.142.2. [AUTHOR REVIEW STATEMENT]

[AUTHOR APPROVAL, ORIGINALITY AND EXCLUSIVE-SUBMISSION STATEMENT]

[AUTHOR DISCLOSURE OF RELATED WORK, FUNDING AND COMPETING INTERESTS]

Sincerely,

[CORRESPONDING AUTHOR'S FULL NAME]

[DEPARTMENT, INSTITUTION, POSTAL ADDRESS]

[CONTACT EMAIL]
