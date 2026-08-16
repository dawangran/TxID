# TxID v1 identity and registry specification

Status: normative draft 1.0.0  
Identity algorithm: SC1 / TF1 / SE1  
Registry schema: 1

## Purpose and non-goals

TxID assigns stable structural identifiers to transcript models on one validated
reference sequence collection. It joins equivalent outputs from independent
long-read RNA-seq callers without making transcript discovery, alignment, or
abundance claims. Exact identity, annotation-relative classification, locus
assignment, and fuzzy similarity are stored as separate facts.

## Coordinates and common transcript model

All external and internal exon coordinates are 1-based, closed intervals.
Coordinates must be positive, `start <= end`, and use strand `+` or `-`. Exons
are normalized into increasing genomic order for validation and serialization.
Overlapping exons in one transcript are invalid. Input record order has no
meaning.

For a multi-exon model, genomic introns are `[left_exon.end + 1,
right_exon.start - 1]`. The canonical intron array is in transcription order:
increasing genomic order on `+` and decreasing genomic order on `-`. TSS and TES
are exact outer exon boundaries: on `+`, TSS is the minimum start and TES the
maximum end; on `-`, TSS is the maximum end and TES the minimum start.

## Reference sequence collection

FASTA sequence is upper-cased and all ASCII whitespace is removed. Each contig
record stores its primary FASTA name, length, and SHA-256 sequence digest. The
assembly fingerprint is SHA-256 over UTF-8 records
`name<TAB>length<TAB>sequence_digest<LF>` sorted bytewise by primary name. The
fingerprint is written as `sha256:<64 lowercase hex characters>`.

An input contig resolves only by an exact primary name or an explicitly supplied
alias table recorded with the registry. Alias targets must exist in the loaded
FASTA and an alias cannot target more than one primary contig. No prefix removal
or other string heuristic is permitted. Every exon end must be within the
resolved contig length.

## Canonical serialization and exact identifiers

Canonical objects are UTF-8 JSON with keys sorted lexicographically, no
insignificant whitespace, JSON Unicode preserved, and integer coordinates. The
following members are normative. Member values shown in angle brackets are
placeholders.

Multi-exon splice chain (`SC1`):

```json
{"algorithm":"SC1","assembly":"<fingerprint>","contig":"<primary name>","introns":[[101,199]],"strand":"+"}
```

Multi-exon transcript form (`TF1`):

```json
{"algorithm":"TF1","assembly":"<fingerprint>","contig":"<primary name>","introns":[[101,199]],"strand":"+","tes":300,"tss":1}
```

Single-exon form (`SE1`):

```json
{"algorithm":"SE1","assembly":"<fingerprint>","contig":"<primary name>","end":300,"start":1,"strand":"+"}
```

The full digest is lowercase SHA-256 of the canonical JSON bytes. The public
digest is the first 24 hexadecimal characters (96 bits). Public identifiers are
`txid:SC1.<public_digest>`, `txid:TF1.<public_digest>`, and
`txid:SE1.<public_digest>`. Registries must store the canonical JSON and full
digest. If one public identifier is presented with a different full digest or
canonical object, the transaction fails; suffix repair is forbidden.

Changing any canonical member, coordinate rule, serialization rule, digest
algorithm, public truncation rule, or identifier meaning requires a new family
version.

## Annotation contexts and classification

An annotation fingerprint is computed from the canonical, sorted list of parsed
reference transcripts containing primary contig, strand, exon intervals,
reference `gene_id`, and reference `transcript_id`. It therefore ignores record
order, comments, and formatting but changes when identity-relevant annotation
content changes. Annotation name and fingerprint are separate.

Classification is relative to one annotation context and does not affect exact
TxIDs:

1. `known`: exact form matches one reference transcript. Rewritten `gene_id` and
   `transcript_id` are the reference values; exact TxIDs are emitted as aliases.
2. `novel_in_known_gene`: no exact form match and exactly one same-strand
   reference gene span overlaps the transcript. The known `gene_id` is retained
   and the form TxID (or SE TxID) becomes `transcript_id`.
3. `ambiguous_gene`: more than one same-strand reference gene span overlaps.
   Candidate genes are emitted explicitly and a registry-managed locus is used;
   no arbitrary known gene is selected.
4. `new_locus`: no same-strand reference gene span overlaps. A
   registry-managed `txid:GL1.<six-digit accession>` is used.

Antisense-only overlap is not evidence for assignment to a known gene. A
read-through overlapping multiple genes is ambiguous. Retained-intron and
single-exon candidates use the same conservative span rule and retain their
separate classification labels in provenance.

GL accessions are persistent registry entities, not exact structural identities.
They are allocated in deterministic transcript sort order within an atomic
import, never reused or renumbered, and may differ between independently created
registries. Exact TxIDs do not.

## Parsing, provenance, and rewriting

GTF and GFF3 exon models enter one shared model. Transcript and gene
meta-features are optional. A malformed coordinate, missing transcript parent,
mixed contig/strand transcript, or overlapping exon is an error. GFF3 multiple
`Parent` values create an exon observation for each parent.

Rewritten GTF is a new file. It contains deterministic gene/transcript/exon order,
preserves safe upstream transcript and exon attributes, and namespaces replaced
values as `txid_original_gene_id` and `txid_original_transcript_id`. It adds
`txid_classification`, `txid_annotation`, `txid_sc` (multi-exon only),
`txid_form`, and optional `txid_gene_candidates`/`txid_fuzzy_cluster`. Exact TxIDs
are always present when a fuzzy cluster is emitted.

Repeated non-identity GTF attributes are retained as ordered key/value
occurrences rather than rejected or silently collapsed. This includes common
GENCODE multi-valued attributes such as `tag` and `ont`. Repeated `gene_id` or
`transcript_id` keys on one GTF feature remain an error because the feature's
identity or parent would be ambiguous. GFF3 continues to express multiple values
inside one attribute value and rejects duplicate keys. `ID`, `Parent`, and
output-reserved keys that cannot retain their original role are preserved with a
`txid_upstream_` prefix. Exon score, phase, source, and all attribute occurrences
are retained in the registry so a retry reproduces the same rewritten bytes.

## Registry and transactions

Schema 1 separates structural objects, observations, annotation classifications,
aliases, loci, fuzzy clusters, and import manifests. Foreign keys are enabled.
An import is completely parsed, normalized, and reference-validated before a
write transaction. Manifest uniqueness makes an identical sample/tool/input
retry idempotent. Failed imports leave no manifest or observations.

The `batch` workflow accepts multiple GTF/GFF3 inputs through a manifest and
processes rows in deterministic key order. Each row remains a separate atomic
import, while immutable assembly metadata and annotation lookup indexes may be
reused across rows and structural objects/observations may be written with
bounded bulk SQL operations. These execution optimizations must not change
canonical objects, exact identifiers, locus assignment, output ordering, or
retry behavior. When optional fuzzy mode is requested for a batch, clustering is
computed once after all rows have imported successfully, and every rewritten
output is generated from that same final fuzzy run.

The `multi-add` workflow is a manifest-free convenience interface for one or
more GTF/GFF3 files that share an upstream tool, annotation context, and optional
format override. Sample names are either supplied one-to-one with the input
paths or derived deterministically from input basenames. It constructs the same
per-input provenance rows and executes the same deterministic batch workflow as
`batch`; therefore it must produce the same registry facts and rewritten bytes
as an equivalent manifest. Heterogeneous tools or annotation contexts require
the manifest interface so their provenance remains explicit.

Input checksum, software version, command options, reference and annotation
fingerprints, sample, upstream tool, path label, original IDs and attributes are
recorded. Database row numbers never appear in exact IDs.

## Optional fuzzy layer

Fuzzy mode requires non-negative splice and end tolerances. Candidate forms must
share assembly, primary contig, strand, exon count, and family. Corresponding
intron boundaries must each be within splice tolerance and TSS/TES (or SE
start/end) within end tolerance.

Clusters use deterministic complete linkage: every pair in a cluster must satisfy
the tolerances. A candidate compatible with more than one existing cluster is an
ambiguous bridge and remains in a singleton cluster with that status recorded.
Thus a chain of pairwise similarities cannot merge incompatible endpoints.
Cluster accessions use `txid:FC1.<six-digit accession>` and record algorithm,
tolerances, membership, and creation software. They never alias or replace exact
IDs.

## Deterministic output and error behavior

Transcript sort key is primary contig, minimum exon start, maximum exon end,
strand, exact form ID, original transcript ID. UTF-8 TSV, JSON, and GTF outputs
use LF line endings and a final newline. Diagnostics go to stderr, data to the
selected file or stdout. Conflicting assembly metadata, unknown contigs,
out-of-range coordinates, digest collisions, malformed input, and registry
integrity failures return non-zero status.

## Cross-assembly relations

TxID v1 does not implement liftover. Registries with different assembly
fingerprints cannot be merged or imported into one another. Future mappings must
retain distinct source and target identifiers and record their validation;
liftover is never exact identity.
