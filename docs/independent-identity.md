# Independently coded exact-identity verifier

`benchmarks/independent_identity.py` provides a second implementation of the
published TxID SC1/TF1/SE1 rules and a separate GTF reader. Both use only the
Python standard library. They do not import TxID's parser, model, reference,
identity or registry modules. The calculator was implemented from the normative
specification, canonical schema and existing golden vectors; production
`identity.py` was not copied. This is implementation separation within the same
project, not a claim of an external replication or independent authorship.

The purpose is to test whether coordinates obtained through a separate reader
produce the same exact canonical bytes and identifiers. This supplements the
older fixed-model scorer, which reuses the production TxID parser. It does not
evaluate transcript discovery, biological truth, gene assignment or abundance.
The original golden vectors and historical benchmark outputs remain unchanged.

## Importable API

From the repository root:

```python
from benchmarks.independent_identity import read_gtf_models, compute_identity

models = read_gtf_models(
    "models.gtf.gz",
    contig_aliases={"1": "chr1"},
    contig_lengths={"chr1": 248956422},
)

# models is keyed by the original, independently parsed transcript_id.
model = models["upstream_transcript_id"]
keys = model.identity(assembly_fingerprint)
assert keys == compute_identity(
    assembly_fingerprint=assembly_fingerprint,
    contig=model.contig,
    strand=model.strand,
    exons=model.exons,
)
form_id = keys["form"]["public_id"]
splice_id = None if keys["splice_chain"] is None else keys["splice_chain"]["public_id"]
```

`IndependentTranscript` is immutable and exposes `transcript_id`, optional
`gene_id`, resolved primary `contig`, `strand` and a tuple of sorted `(start,end)`
exon intervals. `compute_identity` accepts those coordinates directly and
returns `form` plus `splice_chain`. Each non-null identity contains
`canonical_json`, `full_digest` and `public_id`; a single-exon model has a SE1
form and a null splice chain. The public suffix is fixed at the normative first
24 hexadecimal characters of the SHA-256 digest.

`sha256_file(path)` returns the SHA-256 digest of exact file bytes, including
gzip headers. `verify_conformance(path)` returns a report comparing all three
identity fields against the unmodified golden JSON. `IndependentIdentityError`
is a `ValueError` subclass; parser failures identify the source file and line
where possible.

## Reference and parser scope

The assembly fingerprint is a required, syntax-checked input to identity
calculation. This verifier does not independently calculate it from a FASTA or
establish that a supplied fingerprint belongs to a supplied length map. The
calling experiment must record the selected FASTA/context and its checksums.
Consequently, agreement tests the parser and structural identity calculation
under the supplied context, not independent validation of that context's origin.

An optional `contig_lengths` mapping restricts records to declared primary
names and rejects coordinates beyond those lengths. If it is absent, only
positive coordinate and model consistency checks are possible. Explicit aliases
require a length map so their targets can be checked. Aliases must target an
existing primary contig and cannot override another primary name. No `chr`
prefix stripping, sequence inference or alias-chain traversal is performed.

The reader accepts plain or gzip GTF, detecting gzip by its magic bytes rather
than its filename. Transcript meta-features are optional; exons can be shuffled
or interleaved among transcripts. Coordinates remain 1-based closed. On the
negative strand, introns are ordered from high to low genomic position and TSS
and TES use the corresponding outer boundaries. No coordinate is rounded,
shifted or merged.

This is a strict identity-reading subset, not a second full annotation writer:

- Every non-comment record must have nine tab-separated columns, valid positive
  coordinates, a finite score or `.`, a valid strand and a valid phase.
- Attribute strings support quoted strings with escaped quotes, backslashes,
  newlines, carriage returns and tabs, as well as bare single-token values.
  Unknown escapes, unterminated strings and missing delimiters are errors.
- Transcript/exon records require one nonempty `transcript_id` and strand `+`
  or `-`. `gene_id` is optional, but all provided values within a transcript
  must agree. Repeated identity attributes are rejected even when identical.
- Duplicate transcript meta-features, conflicting contigs/strands, missing
  exons, duplicate/overlapping exons and transcript bounds that do not equal
  the exon span are errors. Adjacent exons are also rejected because their
  implied intron is empty; they are never silently merged. These strict checks
  delimit the accepted domain rather than change production parser behavior.
- Repeated nonidentity attributes are accepted. After syntax validation the
  reader intentionally excludes nonidentity attributes, scores and other
  feature types from its returned structural model. Their round-trip
  preservation is outside this verifier's scope.
- GFF3 is outside this reader's scope. Unsupported or malformed input must not
  be converted silently. Experiments should disclose any subset selection or
  rejected models before asserting cross-format agreement.

## CLI and report provenance

Verify the published vectors without adding the production package to the
import path:

```bash
python -I benchmarks/independent_identity.py conformance \
  --vectors tests/conformance/exact-v1.json
```

Calculate keys for a GTF using an explicitly supplied context:

```bash
python -I benchmarks/independent_identity.py gtf \
  --input models.gtf.gz \
  --assembly-fingerprint sha256:0000000000000000000000000000000000000000000000000000000000000000 \
  --contig-lengths contigs.json \
  --contig-aliases aliases.tsv \
  --output independent-model-keys.json
```

The all-zero fingerprint above is an illustration, not a reference identity for
real data. `contigs.json` is a JSON object mapping primary names to positive
integer lengths. `aliases.tsv` has the exact header `alias<TAB>primary` and one
unique alias per row. Omit the alias option if names already match.

Reports record the calculator version, Python version, supplied assembly
fingerprint and SHA-256/byte size/path of the exact input files and calculator.
Conformance reports additionally record the specification and schema checksums.
GTF reports explicitly state whether reference lengths were validated and retain
canonical JSON plus full digests for every model. Reports contain no timestamps,
and repeated execution on the same files produces the same bytes. `--output`
creates an atomic new file and refuses to replace an existing output or input.
Without that option JSON goes to stdout; diagnostics go to stderr. Exit codes
are 0 for success, 1 for a conformance mismatch, and 2 for invalid input or I/O
failure.

## Recorded verification

The new compact conformance result is
`benchmarks/results/independent-identity-2026-09-23.json`, included in the explicit
allowlist of compact verification reports. It covers the existing
three vectors and five non-null identities: positive-strand SC1/TF1,
negative-strand SC1/TF1, and single-exon SE1. Exact canonical JSON, full SHA-256
and public identifiers all match. It does not amend prior benchmark results.

The separate automated tests run with:

```bash
python -m unittest tests.test_independent_identity -v
```

They cover both strands and exon permutations, transcript-end and one-base splice
changes, single-exon intervals, Unicode canonical bytes, context changes,
interleaved GTF parsing, quoted escaping, gzip, explicit aliases and lengths,
malformed/conflicting input rejection, mismatch reporting, output preservation
and deterministic repeated CLI output. An isolated `python -I` subprocess and
an import audit check that execution does not depend on the production package.
These are synthetic conformance/engineering checks; real-data agreement must be
reported separately by the calling experiment.
