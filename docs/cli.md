# TxID command-line reference

All commands are non-interactive. JSON summaries or requested data are written to
stdout/files; diagnostics are written to stderr. Expected input or registry
errors return a non-zero exit status.

## `txid init`

Creates a new schema-v1 registry and refuses to overwrite an existing path.
Required inputs are `--db`, `--fasta`, `--assembly`, `--annotation`, and
`--annotation-name`. `--contig-aliases` accepts a two-column tab-delimited
`alias`/`primary` table. Alias targets must be exact FASTA primary names.

Use `--format gtf` or `--format gff3` only when suffix/content detection is not
sufficient.

## `txid annotation-add`

Registers another annotation release on the registry assembly. The command
revalidates every contig and coordinate and computes a normalized annotation
fingerprint. Annotation context affects known/novel classification, not exact
structural TxIDs.

## `txid add`

Imports one GTF or GFF3 annotation atomically. Required provenance fields are
`--sample`, `--tool`, and `--annotation-name`. `--output-gtf` and `--mapping`
must differ from each other and the input.

Exact mode is the default. Optional fuzzy grouping requires both parameters:

```bash
--fuzzy-splice-tolerance 5 --fuzzy-end-tolerance 25
```

Fuzzy clusters use complete linkage. Exact identifiers remain in every mapping
and rewritten record. An ambiguous bridge becomes its own cluster with
`fuzzy_bridge_status=ambiguous_bridge`.

## `txid batch`

The manifest is tab-delimited and must contain:

```text
input	sample	tool	annotation_name	format
/data/s1.gtf	s1	IsoQuant	GENCODE-v49	gtf
```

`format` is optional. Rows are processed in a deterministic key order; each input
is its own atomic, retry-safe import. This means a late failure does not corrupt
an earlier completed registry import, but the whole manifest is not one global
transaction.

## `txid export`

`--catalog` writes a deterministic cohort TSV containing form/splice identifiers,
gene assignment, classifications, and observation/sample/tool counts.

## `txid validate`

Runs SQLite integrity and foreign-key checks, then independently recomputes every
stored full and public digest. The JSON result includes registry summary counts.

## `txid inspect`

Accepts an `SC1`, `TF1`, `SE1`, `FC1`, `GL1`, or observed output transcript
identifier and returns JSON with canonical structure, observations, locus, or
fuzzy membership as applicable.

## `txid plot`

Writes a dependency-free accessible SVG overview of annotation classifications.

