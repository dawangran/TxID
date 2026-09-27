# Verification report — final editorial revision, 23 September 2026

The contribution-focused manuscript revision passed its source, numerical, translation and artifact checks. No experimental results or software behavior changed in this editorial revision; the software installation/test results below are retained from the preceding validated repair, not newly rerun tests. This is not a statement that author declarations, public deposits or journal pagination are complete.

- The frozen software source archive passed 99 included automated tests, 4 archive-boundary tests, offline wheel installation in a fresh virtual environment, independent golden-vector checks and a separate installed CLI init/add/export/validate smoke. See `submission/release-validation.json` and `software-test-report.json`.
- The real-read case actually executed alignment, IsoQuant, StringTie and four staged imports into two independently initialized TxID 0.1.3 registries. Canonical objects, full digests, observation coordinates, exact joins and reciprocal incremental stability were checked against a separately implemented reader/calculator. All numerical results and the initial execution correction are retained in `real-interop/`.
- The manuscript verifier passed 1562 checks. It independently aggregates the 311 observation rows, distinguishes 110 exact-form matches from 109 chain groups and 115 chain pairs, verifies Table S6 and matrix cells, and confirms that all shared forms are reference-known. Its detailed checks and explicit archive/original-evidence modes are in `verification-report.json`.
- All 44 bilingual blocks match the current English components. Quantities and inline code agree, with one reviewed English-month-name/Chinese-digit conversion. Semantic cross-review retained all version, sampling and known-only limitations. See `translation-verification.json`.
- The English and bilingual DOCX files contain all current paragraphs and the current embedded figure; ZIP/XML checks pass. Office pagination was not rendered.
- Regeneration reproduced all 15 scientific presentation artifacts byte for byte; see `regeneration-check.json`. The figure's checked geometry, transitions, sources, DPI and font embedding remain intact.
- The style scan found no stock emphasis/transitions, unnecessary bold text, body lists or subjective first-person phrases. This does not establish authorship or an absence of AI assistance; actual assistance is disclosed.
- Eleven primary references and their narrow citation roles were checked. The direct transcript-hash precedent, IsoQuant and StringTie methods were added with primary-source evidence.

The English visible-text counts include headings, declarations, caption, accessibility paragraph and reference text; image alt text and URL targets are excluded.

| Component | Words |
|---|---:|
| front-matter.md | 141 |
| 01_introduction.md | 249 |
| 02_implementation.md | 439 |
| 03_evaluation_discussion.md | 580 |
| end-matter.md | 65 |
| figure-caption.md | 227 |
| references.md | 241 |
| main_total | 1942 |
| body_total | 1268 |

The PowerShell research-writing gate is unavailable; the manuscript-specific Python verifier and direct source/scientific review were used. Missing author/Contact/declaration fields remain outside the scientific draft in the author-information form; no values were fabricated. The declared GitHub web URL returned HTTP 404 anonymously, and permanent public software/data archive identifiers remain absent. Complete Conda solving/building, the final release container and the journal's four-page layout have not been verified. Local archives must not be represented as completed public deposits. Historical large-scale 0.1.0 benchmarks were not rerun or relabelled.

The two-sentence Summary and concise title follow the checked Application Note guidance. Initial submission may use Format-Free presentation; final journal pagination remains unrendered. Two separate read-only reviews checked the revised frontmatter/caption/cover letter and the contribution/result claims. Both found no necessary scientific corrections. See `plan/review/final-editorial-review-2026-09-23.md`.
