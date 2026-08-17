# Minimal end-to-end tutorial

This example assumes one FASTA assembly, one reference annotation, and two caller
annotations produced independently.

## 1. Validate the identity context

Use the exact FASTA used for alignment. FASTA sequence is upper-cased and
whitespace-normalized before contig digests are calculated. TxID does not infer
that `1` means `chr1`; provide a reviewed alias table when required:

```text
alias	primary
1	chr1
MT	chrM
```

Initialize once:

```bash
txid init --db cohort.sqlite --fasta genome.fa --assembly GRCh38 \
  --contig-aliases aliases.tsv --annotation reference.gtf \
  --annotation-name release-1
```

## 2. Import independent samples

```bash
txid add --db cohort.sqlite --input sample-a.gtf --sample A --tool FLAIR \
  --annotation-name release-1 --output-gtf sample-a.txid.gtf \
  --mapping sample-a.mapping.tsv

txid add --db cohort.sqlite --input sample-b.gff3 --format gff3 \
  --sample B --tool StringTie --annotation-name release-1 \
  --output-gtf sample-b.txid.gtf --mapping sample-b.mapping.tsv
```

For an exact reference match, rewritten `gene_id` and `transcript_id` are the
reference identifiers. For a novel form overlapping exactly one same-strand
reference gene, `gene_id` is retained and `transcript_id` is the exact `TF1` or
`SE1` identifier. New and ambiguous loci receive persistent `GL1` accessions.

## 3. Build a cohort catalog

```bash
txid export --db cohort.sqlite --catalog cohort.tsv
txid validate --db cohort.sqlite
txid plot --db cohort.sqlite --output cohort.svg
txid plot --db cohort.sqlite --gene MY_GENE_ID --output MY_GENE_ID.svg
```

The gene view labels each exact TxID beside its upstream transcript ID and colors
tracks by fuzzy cluster when fuzzy grouping is available. Same-color structures
that differ substantially can reveal possible over-grouping; near-identical
structures with different colors can reveal possible under-grouping. These are
diagnostic cues, not changes to exact identity.

Use `txid_form` from the mapping files as the join key for a sample-by-transcript
matrix. Do not substitute `FC1` for exact columns unless a tolerance-dependent
analysis is explicitly intended and recorded.

## 4. Add a new annotation release

```bash
txid annotation-add --db cohort.sqlite --annotation reference-v2.gtf \
  --annotation-name release-2
```

Reimporting a structure against release 2 can change `known` versus
`novel_in_known_gene`, and can change preserved reference aliases, while its exact
TxID remains unchanged because annotation version is not part of its canonical
structural object.
