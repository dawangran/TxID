This directory is a compact, unchanged copy of selected files from
benchmark-work/submission-repair-2026-09-23, preserving paths relative to that
case directory. The original JSON/log paths and checksums remain unmodified;
evidence-manifest.json maps each repository source path to its package copy.
Absolute execution paths in historical manifests are provenance, not portable
installation requirements. Source code is supplied separately in the submission
source archive; the evaluation summary records the exact source hashes used.

The original selection.json preserves the first IsoQuant attempt's failure.
execution-correction.json and callers/isoquant.run.json document the successful
retry without --complete_genedb; both the initial and retry logs are retained.
inputs/gencode-v29.identity.chr22.gtf is the exact shared annotation input.

Excluded large or reconstructible artifacts: reference FASTA, raw/selected
FASTQ, BAM and indexes, SQLite registries, tool binary/archive and unrelated
caller abundance products. Their available input hashes remain in the original
manifests. Copied source snapshots under provenance/first-attempt preserve the
initial wrapper behavior. No original case files or frozen code were modified.

New supplemental files are supplement-classification-audit.json and
results-for-s9.md; the companion compact table is ../../tables/table-s6-real-interop.tsv.
These additions do not replace or edit interop-evaluation/summary.json.
