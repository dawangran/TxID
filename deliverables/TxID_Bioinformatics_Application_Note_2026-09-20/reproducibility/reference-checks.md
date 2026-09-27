# Reference and venue checks — 20 September 2026

All eight cited works were checked against the primary article, its author-deposited full text, an authoritative publication record, or the originating project. These checks establish bibliographic identity and the narrow cited facts; they do not verify TxID's own results. No citation count or journal-impact claim is used.

| Citation | Checked record | Checked fact |
|---|---|---|
| ENCODE Project Consortium et al., 2020 | https://www.nature.com/articles/s41586-020-2493-4 | Title, consortium authorship, Nature 583:699–710, public ENCODE data context |
| Jackman et al., 2015 | https://pubmed.ncbi.nlm.nih.gov/26020645/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC4447347/ | Three authors, PLoS ONE 10(5):e0128026, representative sequence-derived k-mer identifiers |
| Liu and Chun, 2026 | https://academic.oup.com/bioinformatics/article/42/1/btaf680/8405383 | Authors, issue year, title, SQANTI3 inputs and extensible isoSeQL database; abundance and metadata tracking |
| Mudge et al., 2025 | https://www.gencodegenes.org/pages/publications.html and https://academic.oup.com/nar/issue/53/D1 | Authors, title, NAR 53:D966–D975 and DOI; reference resource context |
| Pardo-Palacios et al., 2024 | https://www.nature.com/articles/s41592-024-02229-2 and https://pmc.ncbi.nlm.nih.gov/articles/PMC11093726/ | Authors, title, Nature Methods 21:793–797, transcript QC/curation role |
| Pertea and Pertea, 2020 | https://f1000research.com/articles/9-304 | Geo Pertea then Mihaela Pertea, version 2 dated 9 September 2020, DOI ending .2; GFF comparison/tracking role |
| Wagner et al., 2021 | https://pubmed.ncbi.nlm.nih.gov/35311178/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC8929418/ | Authors, title, Cell Genomics 1(2):100027, normalized variation representation and computed-identification precedent |
| Wyman et al., 2020 | https://www.biorxiv.org/content/10.1101/672931v2 | Preprint version 2, posted 24 March 2020, authors, title and TALON's transcript/abundance tracking role |

Publisher instructions checked: https://academic.oup.com/bioinformatics/pages/author-guidelines. The Application Note limit is four journal pages, approximately 2,600 words or 2,000 words plus one figure. The abstract is short and structured (Summary, Availability and Implementation, Contact, Supplementary Information). Contact is deliberately omitted from the scientific draft until the author provides it. Software should be freely available to non-commercial users and maintained for two years following publication.

Direct retrieval of several publisher pages initially failed; indexed primary-page contents and author-deposited articles supplied the relevant material. A direct Crossref API check was attempted from the shell but DNS was unavailable. The report therefore does not claim all DOIs were re-resolved through Crossref in this session.

The TxID GitHub URL is taken from the repository's own pyproject.toml/README. Public access to that URL was not confirmed by this environment: web retrieval failed and shell DNS was unavailable. This is an outstanding submission check, not evidence that the repository is absent or private.

The BibTeX file uses the first three authors plus `and others` for longer author lists; full publisher metadata should be imported if the final journal style requests a longer list.

## 23 September repair — three additional primary references

The current manuscript cites eleven works. The original eight records above remain applicable. Added records were checked against these primary sources:

| Citation | Checked record | Checked fact |
|---|---|---|
| Kabza et al., 2024 | https://www.nature.com/articles/s41467-024-51584-3.pdf and publisher Crossmark | Nature Communications 15:7316; stable transcript-path hashes on the same genome build; 50 bp default annotation end-bin policy; complete hash serialization is not established by this paper |
| Prjibelski et al., 2023 | https://www.nature.com/articles/s41587-022-01565-y | Accurate isoform discovery with IsoQuant using long reads; Nature Biotechnology 41:915–918; DOI year suffix does not replace the 2023 issue year |
| Kovaka et al., 2019 | https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1910-1 and https://ccb.jhu.edu/software/stringtie/ | StringTie2 long-read assembly method; Genome Biology 20:278; the official StringTie page still directs long-read use to this citation; actual executed version 3.0.3 is separately recorded |

The software-access check above is historical. A new anonymous HTTP request on 23 September returned 404 for the declared GitHub web URL (page title “Page not found · GitHub · GitHub”); the API returned 403. See `public-access-check.json`. This establishes unavailable anonymous access in that check, without determining whether the repository is private, renamed or absent. The draft no longer asserts that its source code is publicly accessible there. Local archive validation does not replace the public software/test-data deposition required by the journal.

## Final editorial revision — journal-format recheck

On 23 September 2026, the indexed primary author-guidelines page was retrieved again after direct retrieval failed because the response exceeded the browsing limit. It explicitly permits Format-Free initial submission, with Application Notes guided by approximately 2,000 words plus one figure, and specifies a one- or two-sentence Summary. The final Summary follows that sentence guidance. Missing Contact/author declarations and public access remain separate author-supplied items; unrendered journal-template pagination is not presented as a distinct initial-submission blocker.
