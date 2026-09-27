# TxID — Bioinformatics Application Note, prepared 27 September 2026

The manuscript presents TxID as an independently recomputable and verifiable identity specification. The English text and paragraph-aligned Chinese review copy are synchronized. The English draft contains 1,974 visible-text words, including headings, caption, author completion fields and eleven references; the three body sections contain 1,268 words. The scientific body and figure retain the verified 23 September revision.

- [English manuscript, Word](manuscript.en.docx) · [Markdown](manuscript.en.md)
- [English–Chinese review copy, Word](manuscript.docx) · [Markdown](manuscript.md)
- [Supplement, Word](supplementary.docx) · [Markdown](supplementary.md)
- [Current Chinese handoff](editorial-notes.zh-CN.md)
- [Source archive](submission/TxID-0.1.3-submission-source.zip)
- Complete local review/reproduction ZIP: generate with `python deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/submission/build_materials_archive.py`; the generated bundle is kept outside Git history.
- [Author information still required](submission/author-information-needed.md) · [Allowed author fields](submission/author-fields.json)
- [Cover letter, Word](submission/cover-letter.docx) · [Markdown](submission/cover-letter.md)

The abstract has all four Application Note headings: Summary, Availability and Implementation, Contact, and Supplementary Information. Eight designated fields remain for author names, affiliations, corresponding author, contact email, funding, contributions, human review and competing interests. These are the only remaining manuscript text placeholders. The Word files have active hyperlinks and automatic page numbers; the manuscripts and supplement use 12 pt body text, double spacing and continuous line numbering. The cover letter omits line numbering. This is a Format-Free submission layout; final journal pagination has not been rendered.

The revision credits the Isosceles transcript-hash precedent and defines TxID's contribution as a versioned, reference-aware exact-identity contract and registry. Historical six-TALON results remain explicitly assigned to TxID 0.1.0. A new prespecified real-read case used TxID 0.1.3, IsoQuant and StringTie on 2,765 chr22 primary alignments from the first 100,000 records of public accession ENCFF105WIJ. Two independent registries joined 110 coordinate-identical forms into a two-row, 201-column caller-presence matrix, with no false or missed exact joins. All 110 shared forms were reference-known, and their reference aliases also agree. The case does not demonstrate shared novel discovery, biological correctness or superiority to reference-alias joins.

Supplementary Sections S8–S9 and Tables S5–S6 add tool semantics and the new case. The [compact evidence directory](reproducibility/real-interop/) contains unchanged caller GTFs, the small reference annotation, per-observation coordinates, matrices, mapping/catalog exports, initial failure and successful retry logs, and checked manifests. Raw reads, full FASTA, BAMs and SQLite databases are not in the compact archive. The separately coded parser and identity oracle does not import TxID production modules; assembly fingerprint derivation is outside its independent scope.

Figure 1 remains the validated 178 × 122 mm three-panel composition: identity workflow, annotation-classification transitions and historical six-row matrix dimensions. [PNG](figures/figure1.png), [SVG](figures/figure1.svg), [PDF](figures/figure1.pdf), [1200 dpi TIFF](figures/figure1.tif) and [source data](figures/source-data.json) are supplied. The new two-caller matrix is reported separately in Table S6 and the evidence TSVs. Superseded draft ZIPs and agent task plans remain local and are excluded from this repository snapshot.

Rebuild within the repository or the extracted complete reproduction package:

```bash
MPLCONFIGDIR=/tmp/txid-note-mpl PYTHONDONTWRITEBYTECODE=1 \
  python deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/reproducibility/build_artifacts.py
PYTHONDONTWRITEBYTECODE=1 \
  python deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/reproducibility/verify_manuscript.py
```

The artifact environment is pinned in `reproducibility/requirements-artifacts.txt`. Verification records distinguish numerical/source consistency, translation alignment, package installation and actual experimental execution. They are not an editorial acceptance prediction.

Author completion is separate from access and archival status. On 27 September 2026, anonymous requests to both the declared GitHub web URL and repository API returned HTTP 404; the authenticated repository API identified `dawangran/TxID` as private. See the [current access record](reproducibility/public-access-check-2026-09-27.json). The manuscript supplies the declared software URL and a fixed-revision evidence link, without claiming public access or a completed archival deposit. Anonymous software/test-data access and dedicated immutable software/data archives with verified persistent identifiers remain unresolved.

The AI disclosure records Codex CLI 0.142.2 only as the installed version observed during final preparation. Authors must complete the actual human-review statement and the cover letter's approval and submission declarations. Frozen source archives and historical access reports retain their original scope; they do not certify the current revision or public availability. This package has not been submitted to the journal.
