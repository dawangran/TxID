# Verification report — submission format revision, 27 September 2026

The revised submission documents passed source, numerical, translation and document-structure checks. The scientific body, figure and experimental evidence retain their previous content. This report distinguishes document preparation from completed public deposition or journal submission.

- The manuscript verifier passed 1673 checks, including exact source assembly, numerical recounts, source hashes, paragraph retention, figure embedding and DOCX package validation. Details are in `verification-report.json`.
- The four required abstract headings are present. Eight explicitly authorized author fields remain in the manuscript; their exact allowlist is recorded in `submission/author-fields.json`. Funding, contributions, competing interests and author review are not inferred.
- English and bilingual manuscripts use A4 pages, 12-point body text, double spacing, continuous line numbering, automatic page numbers and active external hyperlinks. The new cover-letter DOCX has page numbers and letter spacing, without manuscript line numbering. XML and content were checked; no Office layout renderer was available.
- Translation source alignment, quantities and inline code are checked in `translation-verification.json`. The Chinese version remains a paragraph-aligned review copy; `manuscript.en.docx` is the submission-preparation file.
- All 16 presentation artifacts reproduce byte for byte on regeneration; see `regeneration-check.json`. This count now includes the cover letter. Figure geometry, sources, DPI and font embedding are unchanged.
- The real-read case retains 311 observations, 110 shared exact forms, 109 shared chain groups and 115 chain pairs. All shared forms are reference-known. The verifier distinguishes packaged evidence from additional comparisons against original local evidence.
- Earlier software validation remains recorded in `submission/release-validation.json` and `software-test-report.json`: 99 included automated tests, four archive-boundary tests, installation, conformance vectors and an installed CLI smoke test. These historical results are not represented as new experimental runs.
- README installation and the complete quickstart were executed using TxID 0.1.3: nine mapping rows, nine catalog rows, successful validation, byte-identical retry exports and unchanged input annotations. Both SVG marks were rendered on light and dark backgrounds; all 35 local README links resolve.

Visible-text counts include headings, author fields, declarations, caption, accessibility text and references; image alt text and URL targets are excluded.

| Component | Words |
|---|---:|
| front-matter.md | 160 |
| 01_introduction.md | 249 |
| 02_implementation.md | 439 |
| 03_evaluation_discussion.md | 580 |
| end-matter.md | 78 |
| figure-caption.md | 227 |
| references.md | 241 |
| main_total | 1974 |
| body_total | 1268 |

Initial submission permits Format-Free presentation; final journal pagination is not claimed. Author details and declarations intentionally remain for the authors to complete. The 27 September anonymous repository checks still returned HTTP 404, while the authenticated connector reported private visibility, despite the user's report that the repository is public. The declared GitHub links identify repository locations; they do not demonstrate anonymous access or permanent archival deposition. No archive DOI or public deposit is invented.

Historical large-scale TxID 0.1.0 results have not been rerun or relabelled. Complete Conda solving/building and a pinned published container remain unverified. The PowerShell writing gate was unavailable; manuscript-specific Python checks and independent source review were used. Scientific claims and author-only placeholders received a separate read-only review with no necessary corrections.
