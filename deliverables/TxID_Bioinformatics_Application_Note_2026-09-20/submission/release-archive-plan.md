# Local release archive and public-deposit plan

Status: updated 27 September 2026. The software archive remains the frozen TxID 0.1.3 snapshot validated on 23 September; the current manuscript and formatting changes are a separate revision. Repository synchronization does not establish anonymous public access, create an archival identifier or submit a manuscript. The real-read evaluation retains its recorded version and source hashes; packaging checks do not assess discovery accuracy. Any later software release must identify its reviewed commit and version.

The current manuscript contains 1,974 visible-text words, including eight author completion fields, with 1,268 words in the three body sections. Its abstract has four required headings. Word artifacts include active hyperlinks and automatic page numbers; manuscripts and supplement have double-spaced 12 pt text and continuous line numbers. The [cover letter](cover-letter.docx) is also supplied in Word without line numbers. Author fields and the final-preparation Codex CLI 0.142.2 qualification are recorded in [author-information-needed.md](author-information-needed.md); actual human review remains for the authors to confirm.

## Materials already present

| Material | Local evidence | Release status |
| --- | --- | --- |
| Software and license | `src/txid/`, `pyproject.toml`, `LICENSE`, `README.md`, `CITATION.cff` | Python package 0.1.3; BSD 3-Clause license; local source available |
| Normative rules and schemas | `docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md`, `schemas/` | Draft specification 1.0.1 records the fuzzy-bridge clarification; identity families and existing clustering behavior unchanged |
| Installation and examples | `docs/tutorial.md`, `docs/cli.md`, `docs/reproducibility.md`, `workflows/` | Local documentation, synthetic smoke-data generator and benchmark workflow present |
| Small test inputs and expected results | `tests/helpers.py`, test modules, `tests/conformance/exact-v1.json` | Synthetic inputs are embedded in tests/helpers; conformance vectors include input models and expected canonical objects, full digests and identifiers |
| Independent identity verifier | `benchmarks/independent_identity.py`, `tests/test_independent_identity.py`, `docs/independent-identity.md` | Local standard-library implementation included; independently parses models and calculates exact identities under a supplied assembly fingerprint, not an independently derived FASTA fingerprint |
| Real-read interoperability runner and evaluator | `benchmarks/run_real_interop_callers.py`, `benchmarks/evaluate_real_interop.py`, `tests/test_real_interop.py`, `tests/test_caller_annotation_policy.py` | Frozen source and small synthetic regression tests included in the final local archive; real-case results belong to the separate evidence deposit |
| Conda recipe | `conda-recipe/meta.yaml`, `tests/test_packaging.py` | Recipe and installed-version assertion corrected to 0.1.3; metadata rendering and wheel-based recipe commands checked; full Conda solve/build remains unverified |
| Container | `Dockerfile`, WDL `docker_image` input | Existing base image pinned by digest and current wheel installed; final 0.1.3 image build and published image digest remain unverified |
| Historical evaluation | Retained compact results in `benchmarks/results/`, input manifests in `workflows/publication_benchmark/config/`, manuscript supplementary files | Historical TxID 0.1.0 labels remain intact; a separate manuscript-data deposit is still needed |

The retained real-read run used a Python 3.12 benchmark environment. Its runner uses `hashlib.file_digest`, available in Python 3.11 and later. This benchmark requirement is separate from the TxID runtime requirement of Python 3.10 and later.

The source archive contains the selected code, schemas, specification, tests, compact input manifests and workflow files. It deliberately excludes historical benchmark result directories, active databases, manuscript files, figures, presentation history and large primary data. This keeps the software snapshot independently installable and testable without mistaking it for the complete manuscript-data deposit.

## Rebuild and inspect the local archive

From the repository root:

```bash
python -B deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/submission/build_release_archive.py
python -B deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/submission/test_release_archive.py
python -B deliverables/TxID_Bioinformatics_Application_Note_2026-09-20/submission/verify_release_archive.py
```

The builder writes `TxID-0.1.3-submission-source.zip`, `TxID-0.1.3-submission-source.manifest.json` and `TxID-0.1.3-submission-source.zip.sha256` in this directory. It uses an explicit allowlist, fixed ZIP timestamps, sorted members, stored compression, normalized permissions and SHA-256 for every payload. The archive includes `RELEASE-MANIFEST.json` and `SHA256SUMS`. It refuses symlinks, path traversal, sensitive path components, unexpected biological-data files, files over 2 MiB and total input over 20 MiB. Only tiny explicitly selected text fixtures may contain FASTA/GTF/GFF3 input. It writes no external service and reads no credentials. Rebuilding replaces only these three named local outputs.

The tests check byte-for-byte archive determinism and payload integrity, plus rejection of unapproved paths, symlinks, oversized files and secret-key material. The extracted snapshot can run its included tests with `PYTHONPATH=src python -B -m unittest discover -s tests -v`. Figure/site tests requiring excluded authoring files are not in this source snapshot. Build outputs are not claimed to be published packages.

Local verification succeeded with Python 3.13.13: the extracted snapshot rebuilt an identical ZIP; a wheel was built offline using the host's available build tools, then installed in a fresh virtual environment; 99 included tests and standalone independent golden-vector verification passed. The CLI ran outside the extracted source directory to initialize a registry, import two synthetic samples, export a catalog and validate the database; both sample mappings contained the same nine exact-form keys. Four archive-boundary tests also passed. See `release-validation.json` for commands and the verified archive digest. Some test subprocesses intentionally use the extracted `src/` tree, so the test suite exercises both archived source and installed code. The external-directory CLI smoke supplies the direct installed-wheel check. This is not a fully isolated build-toolchain or a Conda/container test.

The frozen snapshot was rebuilt after the real-case code freeze and verified from its extracted contents on 23 September. Its internal per-file hashes and retained validation report describe that snapshot, not later working-tree changes. A new software archive requires renewed validation and a distinct release record; preserve the frozen evidence when preparing it. The current repository commit identifies the synchronized source revision and does not turn the earlier ZIP into a tagged release or immutable public deposit.

## Public records and remaining actions

| Record | Current state | Action before submission |
| --- | --- | --- |
| Declared repository | https://github.com/dawangran/TxID; anonymous web and repository API checks on 27 September 2026 returned HTTP 404; the authenticated API identified the repository as private | Establish and verify anonymous software/test-data access before submission. Supplying a URL or synchronizing a private repository does not establish public availability |
| Versioned software archive | Local ZIP prepared; no verified public archive URL or persistent DOI in this package | Deposit the approved software snapshot and small test data in a dedicated archive; record its persistent identifier, release commit and SHA-256 |
| Manuscript-data archive | The manuscript links the evidence at fixed Git commit `35c124be9b0e76cd4d71b39c3ec394d035e5be8d`; a dedicated immutable deposit and persistent identifier remain unverified | Deposit final figure source values, compact reports, command/provenance manifests, reproducibility scripts and checksums; record the verified identifier and relationship to the software release |
| Primary public reads/reference data | External public inputs identified by accession/URL/checksum manifests | Cite source accessions, versions and checksums; do not duplicate large FASTQ/BAM/FASTA files in this source ZIP |
| Historical recovery limits | Sixth v49 input classified without persisted observations; incremental snapshots include reconstruction | Retain these caveats; a newly completed smaller experiment does not retroactively complete historical imports |
| New real-read case | Evaluation owner confirmed successful completion and froze the oracle, evaluator and tests before this final source snapshot | Include the validated final reports and provenance in the separate manuscript-data deposit; the source ZIP includes the evaluator and its small synthetic regression tests, not the real reads or registry databases |
| Immutable runtime | Base-image digest exists; final-image digest and full Conda build not verified | Test Conda installation and final container, record the immutable final digest/platform, and update archival execution inputs |
| Authorship and declarations | Eight explicit author-only manuscript fields and the cover letter's author confirmations remain unfilled | Authors complete names, affiliations, Contact, funding, contributions, competing interests and actual human review; finalize approval and submission declarations in the cover letter |

The 23 September access report remains unchanged: its anonymous web check returned HTTP 404 and its API check returned HTTP 403. The [27 September report](../reproducibility/public-access-check-2026-09-27.json) records the newer anonymous results; authenticated repository visibility is a separate check. No public-deposit or DOI claim is made. Resolve access and dedicated archival records separately from author text completion, then update Availability and Data availability with verified links. The [Bioinformatics author guidelines](https://academic.oup.com/bioinformatics/pages/author-guidelines) remain the submission reference; a local ZIP, fixed Git commit or repository homepage is not itself a verified dedicated archival deposit.
