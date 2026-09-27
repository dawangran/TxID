# Workflow integrations

Choose the integration that matches your input organization:

| Workflow | Use | Documentation |
| --- | --- | --- |
| [txid_multi_add.wdl](txid_multi_add.wdl) | Multiple annotation files sharing one tool and annotation context | [Example inputs](txid_multi_add.inputs.example.json), [Chinese guide](../docs/txid-multi-add-wdl-guide.zh-CN.md) |
| [txid_batch.wdl](txid_batch.wdl) | Cohort imports with per-file sample, tool and annotation provenance | [Example inputs](txid_batch.inputs.example.json), [Chinese guide](../docs/txid-batch-wdl-guide.zh-CN.md) |
| [publication_benchmark/](publication_benchmark/README.md) | Snakemake workflow for caller and identity evaluation | [Benchmark guide](../docs/benchmark.md), [protocol](../docs/design/benchmark-protocol.md) |

The WDL integrations import inputs within one task because isolated scatter tasks
cannot safely mutate a shared SQLite registry. The image input can be overridden;
use the final published image digest for archival runs. See
[installation](../docs/installation.md) for container details.

Benchmark tools have their own environments and are not TxID runtime dependencies.
Generated workflow data belongs in the ignored `benchmark-work/` directory.
