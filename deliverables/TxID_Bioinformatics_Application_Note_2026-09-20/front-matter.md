# TxID: a reference-aware identity specification for long-read transcript models

## Abstract

### Summary

TxID provides a versioned, reference-aware identity specification and registry for GTF/GFF3 transcript models, enabling independent analyses to recompute common structural keys for exact splice chains, transcript forms and single-exon models while recording annotation-dependent classifications and provenance separately. Evaluation across 538,804 ENCODE WTC11 observations showed stable exact identifiers despite input permutation and 29,435 annotation-classification changes; a real-read IsoQuant–StringTie case joined 110 shared reference-known forms across independently initialized registries, with a separately implemented parser and identity calculator confirming the joins.

### Availability and implementation

TxID requires Python 3.10 or later, has no runtime dependencies and uses a BSD 3-Clause license. Source code, the versioned specification, tests and installation instructions are supplied in the accompanying TxID 0.1.3 source archive.

### Supplementary information

Supplementary methods, input accessions, comparison tables and figure source data accompany this manuscript.
