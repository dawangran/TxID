# Local release archive and public-deposit plan

Status: 23 September 2026. This directory prepares a reviewable local software snapshot. It does not publish a release, create a DOI or submit a manuscript. The real-read evaluation is recorded separately; packaging and installation checks are not evidence of discovery accuracy. The snapshot retains package version 0.1.3; a final release must identify the exact reviewed commit and any subsequent version decision.

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

The real-read benchmark runner currently uses `hashlib.file_digest`, available in Python 3.11 and later, and is being run in a Python 3.12 benchmark environment. This benchmark requirement is separate from the TxID runtime requirement of Python 3.10 and later; it must be rechecked if the runner is revised before the final snapshot.

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

The final local snapshot was rebuilt after the real-case code freeze and verified from its extracted contents. Repeat both steps after any further software edits. The internal per-file hashes define this snapshot even when the working tree has uncommitted edits; the archive must not be described as a clean tagged release until a reviewed release commit/tag exists.

## Public records and remaining actions

| Record | Current state | Action before submission |
| --- | --- | --- |
| Declared repository | https://github.com/dawangran/TxID is recorded in README/CITATION; anonymous `curl` GET on 23 September 2026 returned HTTP 404 (an earlier API check returned HTTP 403) | Anonymous access was unavailable in this check. Establish and verify public repository/release access before submission; these responses do not determine whether the repository is private, deleted or otherwise inaccessible |
| Versioned software archive | Local ZIP prepared; no verified public archive URL or persistent DOI in this package | Deposit the approved software snapshot and small test data in a dedicated archive; record its persistent identifier, release commit and SHA-256 |
| Manuscript-data archive | Not yet a verified public deposit | Deposit final figure source values, compact result tables/JSON, exact command/provenance manifests, reproducibility scripts and checksums; record its DOI and relationship to the software release |
| Primary public reads/reference data | External public inputs identified by accession/URL/checksum manifests | Cite source accessions, versions and checksums; do not duplicate large FASTQ/BAM/FASTA files in this source ZIP |
| Historical recovery limits | Sixth v49 input classified without persisted observations; incremental snapshots include reconstruction | Retain these caveats; a newly completed smaller experiment does not retroactively complete historical imports |
| New real-read case | Evaluation owner confirmed successful completion and froze the oracle, evaluator and tests before this final source snapshot | Include the validated final reports and provenance in the separate manuscript-data deposit; the source ZIP includes the evaluator and its small synthetic regression tests, not the real reads or registry databases |
| Immutable runtime | Base-image digest exists; final-image digest and full Conda build not verified | Test Conda installation and final container, record the immutable final digest/platform, and update archival execution inputs |
| Authorship and declarations | No verified author block or author-approved declarations | Complete `author-information-needed.md`, finalize the cover letter and synchronize the manuscript's availability/AI/data statements |

The [Bioinformatics author guidelines](https://academic.oup.com/bioinformatics/pages/author-guidelines) call for the submitted software version and test data to be archived in a dedicated repository, with the archive URL in Availability and implementation, and sufficient information to reproduce results and figures. Persistent software and manuscript-data identifiers in this plan remain pending actual public deposits. Do not substitute the local ZIP path or the mutable repository homepage for a completed archival record.
