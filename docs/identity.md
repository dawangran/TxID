# Identity and reference context

TxID separates exact structure, annotation-relative classification, locus
assignment and optional fuzzy similarity. This guide explains the output; the
[identity v1 specification](spec/identity-v1.md) defines the normative rules.

## Identifier families

| Family | Object represented | Assignment and portability |
| --- | --- | --- |
| `txid:SC1.<digest>` | Complete intron chain of a multi-exon model; outer transcript ends are excluded. | Computed exact key; reproducible across registries sharing the reference context. |
| `txid:TF1.<digest>` | Multi-exon form: intron chain plus model-defined transcription start (TSS) and end (TES). | Computed exact key; changing an end changes TF1 even when SC1 is unchanged. |
| `txid:SE1.<digest>` | Single-exon interval and strand. | Computed exact key; single-exon models are classified separately. |
| `txid:FC1.<accession>` | Fuzzy cluster with recorded splice/end tolerances and membership. | Registry-managed; not a portable exact structural key. |
| `txid:GL1.<accession>` | New or unresolved gene/locus assignment. | Registry-managed; not a portable exact structural key. |

Exact objects include the algorithm family, reference fingerprint, primary contig
name and strand. Coordinates are **1-based closed**. Canonical JSON fixes key
ordering and serialization; SHA-256 supplies the full digest and its first 24
hexadecimal characters form the public suffix. Registries retain the complete
canonical object and digest. A conflicting public identifier is an error and
cannot be repaired by adding a suffix.

Sample labels, tool names, annotation releases and import order are excluded from
exact hash inputs. Changing identity semantics requires a new algorithm-family
version. The [schemas](../schemas/) and [conformance vectors](../tests/conformance/exact-v1.json)
make the contract testable independently of the CLI and registry.

## Reference and annotation context

The reference fingerprint includes normalized sequence digests, contig lengths
and **primary FASTA names**. Matching assembly labels alone are insufficient:
renaming primary contigs or changing the sequence collection changes the context.
Input contigs must match primary names or use an explicit alias table whose
targets exist in the selected FASTA. TxID does not infer aliases by stripping
`chr`. Registries with different reference fingerprints cannot be combined as
exact identity contexts; cross-assembly liftover is not implemented.

Annotation versions have separate normalized fingerprints. On the same reference
context, changing annotation release can change classification and reference
aliases while leaving exact structural IDs unchanged:

| Relationship to the selected annotation | Rewritten identifiers |
| --- | --- |
| Exact reference-form match | Preserve the selected reference `gene_id` and `transcript_id`; record TxIDs separately. |
| No exact match; overlap with one same-strand reference gene span | Retain that reference `gene_id`; use the TF1 or SE1 transcript ID. |
| Multiple candidate genes or no same-strand gene overlap | Keep ambiguous candidates explicit where present; assign a GL1 locus and a TF1 or SE1 transcript ID. |

Gene-span overlap is an assignment rule, not evidence of biological gene
membership. Novelty is relative to the annotation and is never part of a
permanent exact ID.

## Exact joins and fuzzy grouping

Use **`txid_form`** as the exact-form join key in mapping tables: TF1 for multi-exon
models and SE1 for single-exon models. Use `txid_sc` to compare complete multi-exon
splice chains; single-exon models have no SC1. Mapping and catalog counts describe
model observations, not read counts or TPM. Abundance measurements must come from
the upstream analysis.

Exact mode is the default and needs no tolerance. Optional fuzzy grouping requires
both `--fuzzy-splice-tolerance` and `--fuzzy-end-tolerance` in base pairs. Every
pair in a cluster must meet the tolerances, and ambiguous bridges are marked.
Exact IDs remain available alongside FC1 assignments. Fuzzy grouping does not
change exact identifiers or establish biological equivalence.

## Import and export guarantees

TxID reads GTF and GFF3, including gzip-compressed inputs and exon models lacking
transcript meta-features. All callers enter the same internal model. Unknown
contigs, inconsistent models, overlapping exons and out-of-range coordinates are
rejected. GFF3 inputs are rewritten as GTF; this is not a byte-preserving format
conversion. Upstream attributes are retained, with reserved attributes namespaced
where necessary.

Each input is staged and validated before its atomic registry transaction. An
identical sample/tool/input retry is idempotent. A batch consists of separate
per-input transactions, so a later failed input does not undo earlier completed
imports. Exports from a fixed registry are sorted deterministically; exact-ID
invariance does not imply identical locus assignments or complete catalog bytes
across different import histories.

Exact IDs preserve differences in caller boundaries. They do not establish
transcript discovery accuracy, biological correctness, abundance accuracy or the
validity of locus assignments. See the [CLI reference](cli.md) for commands and
fields and the [benchmark guide](benchmark.md) for evaluated claims.
