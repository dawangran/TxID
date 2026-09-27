# Machine-readable schemas

| File | Contract |
| --- | --- |
| [canonical-object.schema.json](canonical-object.schema.json) | Exact SC1, TF1 and SE1 canonical objects |
| [mapping-table.schema.json](mapping-table.schema.json) | Mapping-table fields and identifier representations |
| [registry-v1.sql](registry-v1.sql) | SQLite registry schema, version 1 |

The [identity specification](../docs/spec/identity-v1.md) defines semantics and
canonical serialization. The [conformance vectors](../tests/conformance/exact-v1.json)
provide expected canonical objects, full digests and public IDs. Schema validation
alone does not establish that a digest matches its object or that an annotation
matches its reference context.

Changes to identity meaning require a new algorithm family. Registry changes need
an explicit schema version and migration; see [Contributing](../CONTRIBUTING.md).
