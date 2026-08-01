# TxID

TxID is a reference-aware identity and interoperability layer for transcript
models produced by long-read RNA-seq tools. It gives the same exact structure the
same stable identifier across samples, runs, input order, and supported GTF/GFF3
producers. It does **not** discover transcripts, align reads, or estimate
abundance.

Current release: `0.1.0` (research alpha). The normative identity definition is
the [TxID v1 specification](docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md).

## Identity families

- `txid:SC1.<digest>` — exact multi-exon splice chain.
- `txid:TF1.<digest>` — exact multi-exon form including TSS and TES.
- `txid:SE1.<digest>` — exact single-exon form.
- `txid:FC1.<accession>` — optional tolerance-dependent fuzzy cluster.
- `txid:GL1.<accession>` — registry-managed new or unresolved locus.

Known reference `gene_id` and `transcript_id` values remain unchanged in
rewritten annotations. Exact TxIDs are added as aliases. Novelty is stored as an
annotation-relative classification and never appears in a permanent structural
identifier.

## Install

TxID requires Python 3.10 or newer and has no runtime dependencies.

```bash
python -m pip install .
txid --version
```

For development without installation:

```bash
PYTHONPATH=src python -m txid --help
PYTHONPATH=src python -m unittest discover -s tests -v
```

## Minimal workflow

```bash
txid init \
  --db cohort.sqlite \
  --fasta genome.fa \
  --assembly GRCh38 \
  --annotation gencode.v49.gtf \
  --annotation-name GENCODE-v49

txid add \
  --db cohort.sqlite \
  --input sample01.isoquant.gtf \
  --sample sample01 \
  --tool IsoQuant \
  --annotation-name GENCODE-v49 \
  --output-gtf sample01.txid.gtf \
  --mapping sample01.mapping.tsv

txid export --db cohort.sqlite --catalog cohort.catalog.tsv
txid validate --db cohort.sqlite
```

Inputs are never modified in place. `add` parses and reference-validates the
entire annotation before opening a transaction. A failed import leaves no
manifest or observations, and an identical retry is idempotent.

See the [tutorial](docs/tutorial.md), [CLI reference](docs/cli.md), and
[benchmark guide](docs/benchmark.md).

## Reproduce the synthetic evaluation

```bash
PYTHONPATH=src python benchmarks/run_simulation.py \
  --output-dir benchmarks/results
```

The documented evaluation uses seed `20260730`, 120 truth transcripts, six simulated
samples, four explicitly labelled `-like` callers, and a 5% caller-boundary error
rate. It is a functional evaluation, not a biological dataset. The same fixed
truth manifest was used for real executions of gffcompare 0.12.10, isoSeQL 1.0.1,
and TALON 6.0.1; see the benchmark guide for stage scope, commands, results, and
the explicit note that SQANTI3 preprocessing was not run.

For the newer end-to-end publication workflow, including real IsoQuant, FLAIR,
StringTie, TALON, gffcompare, SQANTI3, and isoSeQL execution, see the
[workflow guide](workflows/publication_benchmark/README.md) and the
[Tier-0 execution report](docs/benchmark/publication-smoke-2026-07-31.md).
Tier 0 verifies software integration. The separate
[ENCODE full-genome GTF report](docs/benchmark/encode-gtf-full-2026-07-31.md)
evaluates released fixed models. The current Technical Note is scoped to
deterministic identity and interoperability of those models, not discovery
accuracy or full-genome caller superiority.

## Scope and important limitations

- A registry contains one reference assembly fingerprint. Cross-assembly
  liftover is not exact identity and is not implemented in v1.
- Exact transcript ends are intentional. Boundary errors yield distinct `TF1` or
  `SE1` identifiers; optional fuzzy clusters are reported separately.
- Gene assignment is deliberately conservative: one same-strand overlapping
  reference gene is assigned; multiple genes remain explicit and use a `GL1`
  locus.
- The current alpha includes synthetic validation, version-pinned integration
  tests, and a documented six-sample public full-genome fixed-model benchmark.
  Generated inputs and result payloads are deliberately excluded from this
  source-only repository. TxID does not evaluate discovery sensitivity or
  abundance accuracy from raw reads.

## License

BSD 3-Clause, an OSI-approved license. See [LICENSE](LICENSE).
