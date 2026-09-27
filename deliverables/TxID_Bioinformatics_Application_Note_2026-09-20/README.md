# TxID — Bioinformatics Application Note, revised 23 September 2026

The final editorial revision foregrounds the independently recomputable and verifiable identity specification in the title, abstract, methods and conclusions. The English submission text and paragraph-aligned Chinese review copy are synchronized. The English draft currently contains 1,942 visible-text words, including headings, caption, declarations and eleven references. Final journal pagination has not been rendered.

- [English manuscript, Word](manuscript.en.docx) · [Markdown](manuscript.en.md)
- [English–Chinese review copy, Word](manuscript.docx) · [Markdown](manuscript.md)
- [Supplement, Word](supplementary.docx) · [Markdown](supplementary.md)
- [Current Chinese handoff](editorial-notes.zh-CN.md)
- [Source archive](submission/TxID-0.1.3-submission-source.zip)
- Complete local review/reproduction ZIP: generate with `python deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/submission/build_materials_archive.py`; the generated bundle is kept outside Git history.
- [Author information still required](submission/author-information-needed.md) · [Cover-letter draft](submission/cover-letter.md)

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

Authors, Contact, funding, competing interests and author approval are still required. The declared GitHub address returned HTTP 404 to an anonymous check on 23 September; the cause is not inferred. Public software/test-data access and permanent archive links must be completed before submission. The journal permits Format-Free initial submission; journal-style pagination remains to be checked when preparing the formatted version. This synchronization does not submit the manuscript to the journal or make a public release.

GitHub synchronization, 27 September 2026: the authenticated repository API confirms that `dawangran/TxID` is private. Pushing these materials does not make the software publicly accessible or satisfy the journal public-availability requirement. The dated anonymous-access report remains an unchanged historical record. Frozen source archives record the earlier tested snapshot; the current Git commit is the source of truth for this repository revision.
