"""Versioned SQLite registry with atomic, retry-safe imports."""

from __future__ import annotations

import hashlib
import json
import sqlite3
import sys
from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator

from .classify import ReferenceEntry
from .errors import IdentityCollisionError, RegistryError
from .identity import PUBLIC_DIGEST_LENGTH
from .fuzzy import cluster_forms
from .models import (
    Assignment,
    ContigRecord,
    Exon,
    IdentityBundle,
    SequenceCollection,
    StructuralIdentity,
    TranscriptModel,
    freeze_attributes,
)

SCHEMA_VERSION = 1


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _schema_sql() -> str:
    root = Path(__file__).resolve().parents[2]
    schema_path = root / "schemas" / "registry-v1.sql"
    if schema_path.exists():
        return schema_path.read_text(encoding="utf-8")
    # Installed wheels place the SQL beside this module through package-data.
    packaged = Path(__file__).with_name("registry-v1.sql")
    if packaged.exists():
        return packaged.read_text(encoding="utf-8")
    data_file = Path(sys.prefix) / "share" / "txid" / "registry-v1.sql"
    if data_file.exists():
        return data_file.read_text(encoding="utf-8")
    raise RegistryError("registry schema SQL is missing from the installation")


class Registry:
    def __init__(self, path: str | Path, *, read_only: bool = False):
        self.path = Path(path)
        if read_only:
            uri = f"file:{self.path.resolve()}?mode=ro"
            self.connection = sqlite3.connect(uri, uri=True, isolation_level=None)
        else:
            self.connection = sqlite3.connect(self.path, isolation_level=None)
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("PRAGMA foreign_keys = ON")
        self.connection.execute("PRAGMA busy_timeout = 30000")

    def close(self) -> None:
        self.connection.close()

    def __enter__(self) -> "Registry":
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    @classmethod
    def create(
        cls,
        path: str | Path,
        reference: SequenceCollection,
        *,
        fasta_checksum: str,
        alias_checksum: str | None = None,
        software_version: str,
    ) -> "Registry":
        target = Path(path)
        if target.exists():
            raise RegistryError(f"refusing to overwrite existing registry {target}")
        target.parent.mkdir(parents=True, exist_ok=True)
        registry = cls(target)
        try:
            registry.connection.executescript(_schema_sql())
            with registry.transaction():
                registry.connection.executemany(
                    "INSERT INTO schema_meta(key, value) VALUES (?, ?)",
                    [
                        ("schema_version", str(SCHEMA_VERSION)),
                        ("identity_version", "SC1/TF1/SE1"),
                        ("software_version_created", software_version),
                    ],
                )
                registry.connection.execute(
                    """
                    INSERT INTO reference_context(
                        id, assembly_name, fingerprint, fasta_checksum, alias_checksum, created_utc
                    ) VALUES (1, ?, ?, ?, ?, ?)
                    """,
                    (
                        reference.name,
                        reference.fingerprint,
                        fasta_checksum,
                        alias_checksum,
                        utc_now(),
                    ),
                )
                registry.connection.executemany(
                    "INSERT INTO contig(name, length, sequence_digest) VALUES (?, ?, ?)",
                    [
                        (contig.name, contig.length, contig.sequence_digest)
                        for contig in reference.contigs
                    ],
                )
                registry.connection.executemany(
                    "INSERT INTO contig_alias(alias, primary_name) VALUES (?, ?)",
                    list(reference.aliases),
                )
        except Exception:
            registry.close()
            if target.exists():
                target.unlink()
            raise
        return registry

    @contextmanager
    def transaction(self) -> Iterator[None]:
        self.connection.execute("BEGIN IMMEDIATE")
        try:
            yield
        except Exception:
            # SQLite can roll a transaction back automatically after fatal I/O
            # errors such as a full filesystem.  Do not mask the actionable
            # original exception with "cannot rollback - no transaction is
            # active".
            if self.connection.in_transaction:
                self.connection.execute("ROLLBACK")
            raise
        else:
            self.connection.execute("COMMIT")

    def check_schema(self) -> None:
        try:
            row = self.connection.execute(
                "SELECT value FROM schema_meta WHERE key = 'schema_version'"
            ).fetchone()
        except sqlite3.Error as error:
            raise RegistryError(f"{self.path} is not a TxID registry: {error}") from error
        if row is None or int(row["value"]) != SCHEMA_VERSION:
            found = None if row is None else row["value"]
            raise RegistryError(
                f"unsupported registry schema {found!r}; this release requires {SCHEMA_VERSION}"
            )

    def reference(self) -> SequenceCollection:
        self.check_schema()
        context = self.connection.execute(
            "SELECT assembly_name, fingerprint FROM reference_context WHERE id = 1"
        ).fetchone()
        if context is None:
            raise RegistryError("registry has no reference context")
        contigs = tuple(
            ContigRecord(row["name"], row["length"], row["sequence_digest"])
            for row in self.connection.execute(
                "SELECT name, length, sequence_digest FROM contig ORDER BY name"
            )
        )
        aliases = tuple(
            (row["alias"], row["primary_name"])
            for row in self.connection.execute(
                "SELECT alias, primary_name FROM contig_alias ORDER BY alias"
            )
        )
        return SequenceCollection(
            name=context["assembly_name"],
            fingerprint=context["fingerprint"],
            contigs=contigs,
            aliases=aliases,
        )

    def register_structural(self, identity: StructuralIdentity) -> None:
        row = self.connection.execute(
            "SELECT full_digest, canonical_json, family FROM structural_object WHERE public_id = ?",
            (identity.public_id,),
        ).fetchone()
        if row is not None:
            if (
                row["full_digest"] != identity.full_digest
                or row["canonical_json"] != identity.canonical_json
                or row["family"] != identity.family
            ):
                raise IdentityCollisionError(
                    f"public identifier collision for {identity.public_id}: stored and incoming "
                    "canonical objects have different full SHA-256 digests"
                )
            return
        self.connection.execute(
            """
            INSERT INTO structural_object(
                public_id, family, public_digest, full_digest, canonical_json
            ) VALUES (?, ?, ?, ?, ?)
            """,
            (
                identity.public_id,
                identity.family,
                identity.public_digest,
                identity.full_digest,
                identity.canonical_json,
            ),
        )

    def add_annotation(
        self,
        *,
        name: str,
        fingerprint: str,
        input_checksum: str,
        transcripts: Iterable[tuple[TranscriptModel, IdentityBundle]],
    ) -> int:
        existing = self.connection.execute(
            "SELECT id, fingerprint FROM annotation_context WHERE name = ?", (name,)
        ).fetchone()
        if existing is not None:
            if existing["fingerprint"] != fingerprint:
                raise RegistryError(
                    f"annotation name {name!r} already identifies a different fingerprint"
                )
            return int(existing["id"])
        cursor = self.connection.execute(
            """
            INSERT INTO annotation_context(name, fingerprint, input_checksum, created_utc)
            VALUES (?, ?, ?, ?)
            """,
            (name, fingerprint, input_checksum, utc_now()),
        )
        annotation_id = int(cursor.lastrowid)
        for transcript, identities in transcripts:
            if identities.splice_chain is not None:
                self.register_structural(identities.splice_chain)
            self.register_structural(identities.form)
            self.connection.execute(
                """
                INSERT INTO reference_transcript(
                    annotation_id, reference_gene_id, reference_transcript_id,
                    contig, strand, start, end, exons_json, attributes_json,
                    splice_chain_id, form_id
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    annotation_id,
                    transcript.original_gene_id,
                    transcript.original_transcript_id,
                    transcript.contig,
                    transcript.strand,
                    transcript.start,
                    transcript.end,
                    json.dumps([[e.start, e.end] for e in transcript.exons], separators=(",", ":")),
                    json.dumps(dict(transcript.attributes), sort_keys=True, separators=(",", ":")),
                    identities.splice_chain.public_id if identities.splice_chain else None,
                    identities.form.public_id,
                ),
            )
        return annotation_id

    def annotation(self, name: str) -> sqlite3.Row:
        row = self.connection.execute(
            "SELECT * FROM annotation_context WHERE name = ?", (name,)
        ).fetchone()
        if row is None:
            available = ", ".join(
                item["name"]
                for item in self.connection.execute(
                    "SELECT name FROM annotation_context ORDER BY name"
                )
            )
            raise RegistryError(
                f"annotation context {name!r} is not registered"
                + (f"; available: {available}" if available else "")
            )
        return row

    def references(self, annotation_id: int) -> list[ReferenceEntry]:
        return [
            ReferenceEntry(
                gene_id=row["reference_gene_id"],
                transcript_id=row["reference_transcript_id"],
                contig=row["contig"],
                strand=row["strand"],
                start=row["start"],
                end=row["end"],
                form_id=row["form_id"],
                splice_chain_id=row["splice_chain_id"],
            )
            for row in self.connection.execute(
                """
                SELECT reference_gene_id, reference_transcript_id, contig, strand,
                       start, end, form_id, splice_chain_id
                FROM reference_transcript WHERE annotation_id = ?
                ORDER BY reference_transcript_id
                """,
                (annotation_id,),
            )
        ]

    def find_import(
        self, annotation_id: int, sample: str, tool: str, input_checksum: str
    ) -> int | None:
        row = self.connection.execute(
            """
            SELECT id FROM import_manifest
            WHERE annotation_id = ? AND sample = ? AND upstream_tool = ? AND input_checksum = ?
            """,
            (annotation_id, sample, tool, input_checksum),
        ).fetchone()
        return None if row is None else int(row["id"])

    def allocate_locus(
        self,
        *,
        contig: str,
        strand: str,
        start: int,
        end: int,
        status: str,
        gene_candidates: tuple[str, ...],
    ) -> tuple[int, str]:
        candidates_json = json.dumps(gene_candidates, separators=(",", ":"))
        rows = self.connection.execute(
            """
            SELECT accession, public_id, status, gene_candidates_json
            FROM locus
            WHERE contig = ? AND strand = ? AND start <= ? AND end >= ?
            ORDER BY accession
            """,
            (contig, strand, end, start),
        ).fetchall()
        compatible = [
            row
            for row in rows
            if row["status"] == status
            and (status != "ambiguous_gene" or row["gene_candidates_json"] == candidates_json)
        ]
        if len(compatible) == 1:
            return int(compatible[0]["accession"]), compatible[0]["public_id"]
        cursor = self.connection.execute(
            """
            INSERT INTO locus(
                public_id, contig, strand, start, end, status,
                gene_candidates_json, created_utc
            ) VALUES ('pending', ?, ?, ?, ?, ?, ?, ?)
            """,
            (contig, strand, start, end, status, candidates_json, utc_now()),
        )
        accession = int(cursor.lastrowid)
        public_id = f"txid:GL1.{accession:06d}"
        self.connection.execute(
            "UPDATE locus SET public_id = ? WHERE accession = ?",
            (public_id, accession),
        )
        return accession, public_id

    def store_import(
        self,
        *,
        annotation_id: int,
        sample: str,
        tool: str,
        source_label: str,
        input_checksum: str,
        options: dict[str, Any],
        software_version: str,
        provisional: list[tuple[TranscriptModel, IdentityBundle, str, str | None, str, tuple[str, ...]]],
        annotation_name: str,
    ) -> tuple[int, list[Assignment], bool]:
        existing = self.find_import(annotation_id, sample, tool, input_checksum)
        if existing is not None:
            return existing, self.load_assignments(existing, annotation_name), False
        assignments: list[Assignment] = []
        with self.transaction():
            cursor = self.connection.execute(
                """
                INSERT INTO import_manifest(
                    annotation_id, sample, upstream_tool, source_label, input_checksum,
                    options_json, software_version, created_utc
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    annotation_id,
                    sample,
                    tool,
                    source_label,
                    input_checksum,
                    json.dumps(options, sort_keys=True, separators=(",", ":")),
                    software_version,
                    utc_now(),
                ),
            )
            import_id = int(cursor.lastrowid)
            for ordinal, item in enumerate(provisional, 1):
                transcript, identities, label, output_gene, output_transcript, gene_candidates = item
                if identities.splice_chain is not None:
                    self.register_structural(identities.splice_chain)
                self.register_structural(identities.form)
                locus_accession: int | None = None
                if output_gene is None:
                    locus_accession, output_gene = self.allocate_locus(
                        contig=transcript.contig,
                        strand=transcript.strand,
                        start=transcript.start,
                        end=transcript.end,
                        status=label if label in {"new_locus", "ambiguous_gene"} else "new_locus",
                        gene_candidates=gene_candidates,
                    )
                assignment = Assignment(
                    transcript=transcript,
                    identities=identities,
                    classification=label,
                    output_gene_id=output_gene,
                    output_transcript_id=output_transcript,
                    annotation_name=annotation_name,
                    gene_candidates=gene_candidates,
                )
                assignments.append(assignment)
                self.connection.execute(
                    """
                    INSERT INTO observation(
                        import_id, ordinal, original_gene_id, original_transcript_id,
                        output_gene_id, output_transcript_id, classification,
                        gene_candidates_json, contig, strand, start, end,
                        exons_json, attributes_json, exon_attributes_json,
                        source, source_format, exon_metadata_json,
                        splice_chain_id, form_id, locus_accession, fuzzy_cluster_id
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, NULL)
                    """,
                    (
                        import_id,
                        ordinal,
                        transcript.original_gene_id,
                        transcript.original_transcript_id,
                        output_gene,
                        output_transcript,
                        label,
                        json.dumps(gene_candidates, separators=(",", ":")),
                        transcript.contig,
                        transcript.strand,
                        transcript.start,
                        transcript.end,
                        json.dumps([[e.start, e.end] for e in transcript.exons], separators=(",", ":")),
                        json.dumps(dict(transcript.attributes), sort_keys=True, separators=(",", ":")),
                        json.dumps([dict(e.attributes) for e in transcript.exons], sort_keys=True, separators=(",", ":")),
                        transcript.source,
                        transcript.source_format,
                        json.dumps(
                            [[e.source, e.score, e.phase] for e in transcript.exons],
                            separators=(",", ":"),
                        ),
                        identities.splice_chain.public_id if identities.splice_chain else None,
                        identities.form.public_id,
                        locus_accession,
                    ),
                )
        return import_id, assignments, True

    def _identity(self, public_id: str | None) -> StructuralIdentity | None:
        if public_id is None:
            return None
        row = self.connection.execute(
            "SELECT * FROM structural_object WHERE public_id = ?", (public_id,)
        ).fetchone()
        if row is None:
            raise RegistryError(f"observation refers to missing structural object {public_id}")
        return StructuralIdentity(
            family=row["family"],
            public_id=row["public_id"],
            public_digest=row["public_digest"],
            full_digest=row["full_digest"],
            canonical_json=row["canonical_json"],
        )

    def load_assignments(self, import_id: int, annotation_name: str) -> list[Assignment]:
        rows = self.connection.execute(
            "SELECT * FROM observation WHERE import_id = ? ORDER BY ordinal", (import_id,)
        ).fetchall()
        assignments: list[Assignment] = []
        for row in rows:
            exon_pairs = json.loads(row["exons_json"])
            exon_attributes = json.loads(row["exon_attributes_json"])
            exon_metadata = json.loads(row["exon_metadata_json"])
            exons = tuple(
                Exon(
                    start,
                    end,
                    source=metadata[0],
                    score=metadata[1],
                    phase=metadata[2],
                    attributes=freeze_attributes(attrs),
                )
                for (start, end), attrs, metadata in zip(
                    exon_pairs, exon_attributes, exon_metadata
                )
            )
            transcript = TranscriptModel(
                contig=row["contig"],
                strand=row["strand"],
                exons=exons,
                original_transcript_id=row["original_transcript_id"],
                original_gene_id=row["original_gene_id"],
                source=row["source"],
                source_format=row["source_format"],
                attributes=freeze_attributes(json.loads(row["attributes_json"])),
            )
            form = self._identity(row["form_id"])
            assert form is not None
            assignments.append(
                Assignment(
                    transcript=transcript,
                    identities=IdentityBundle(
                        form=form, splice_chain=self._identity(row["splice_chain_id"])
                    ),
                    classification=row["classification"],
                    output_gene_id=row["output_gene_id"],
                    output_transcript_id=row["output_transcript_id"],
                    annotation_name=annotation_name,
                    gene_candidates=tuple(json.loads(row["gene_candidates_json"])),
                    fuzzy_cluster=row["fuzzy_cluster_id"],
                    fuzzy_bridge_status=self._fuzzy_status(row["fuzzy_cluster_id"]),
                )
            )
        return assignments

    def _fuzzy_status(self, public_id: str | None) -> str | None:
        if public_id is None:
            return None
        row = self.connection.execute(
            "SELECT bridge_status FROM fuzzy_cluster WHERE public_id = ?", (public_id,)
        ).fetchone()
        return None if row is None else str(row["bridge_status"])

    def run_fuzzy(
        self,
        *,
        splice_tolerance: int,
        end_tolerance: int,
        software_version: str,
    ) -> dict[str, str]:
        if splice_tolerance < 0 or end_tolerance < 0:
            raise RegistryError("fuzzy tolerances must be non-negative")
        rows = self.connection.execute(
            """
            SELECT DISTINCT s.*
            FROM structural_object s
            JOIN observation o ON o.form_id = s.public_id
            WHERE s.family IN ('TF1', 'SE1')
            ORDER BY s.public_id
            """
        ).fetchall()
        identities = [
            StructuralIdentity(
                family=row["family"],
                public_id=row["public_id"],
                public_digest=row["public_digest"],
                full_digest=row["full_digest"],
                canonical_json=row["canonical_json"],
            )
            for row in rows
        ]
        input_fingerprint = "sha256:" + hashlib.sha256(
            "".join(f"{item.public_id}\t{item.full_digest}\n" for item in identities).encode("utf-8")
        ).hexdigest()
        clusters = cluster_forms(
            identities,
            splice_tolerance=splice_tolerance,
            end_tolerance=end_tolerance,
        )
        mapping: dict[str, str] = {}
        with self.transaction():
            run = self.connection.execute(
                """
                SELECT id FROM fuzzy_run
                WHERE algorithm = 'complete-linkage-v1'
                  AND splice_tolerance = ? AND end_tolerance = ?
                  AND input_fingerprint = ?
                """,
                (splice_tolerance, end_tolerance, input_fingerprint),
            ).fetchone()
            if run is None:
                cursor = self.connection.execute(
                    """
                    INSERT INTO fuzzy_run(
                        algorithm, splice_tolerance, end_tolerance,
                        input_fingerprint, software_version, created_utc
                    ) VALUES ('complete-linkage-v1', ?, ?, ?, ?, ?)
                    """,
                    (
                        splice_tolerance,
                        end_tolerance,
                        input_fingerprint,
                        software_version,
                        utc_now(),
                    ),
                )
                run_id = int(cursor.lastrowid)
                for cluster in clusters:
                    signature = hashlib.sha256("\n".join(cluster.members).encode("utf-8")).hexdigest()
                    cursor = self.connection.execute(
                        """
                        INSERT INTO fuzzy_cluster(
                            public_id, fuzzy_run_id, signature, bridge_status
                        ) VALUES ('pending', ?, ?, ?)
                        """,
                        (run_id, signature, cluster.bridge_status),
                    )
                    accession = int(cursor.lastrowid)
                    public_id = f"txid:FC1.{accession:06d}"
                    self.connection.execute(
                        "UPDATE fuzzy_cluster SET public_id = ? WHERE accession = ?",
                        (public_id, accession),
                    )
                    self.connection.executemany(
                        "INSERT INTO fuzzy_member(fuzzy_cluster_accession, form_id) VALUES (?, ?)",
                        [(accession, form_id) for form_id in cluster.members],
                    )
                    mapping.update({form_id: public_id for form_id in cluster.members})
            else:
                run_id = int(run["id"])
                for row in self.connection.execute(
                    """
                    SELECT c.public_id, m.form_id
                    FROM fuzzy_cluster c
                    JOIN fuzzy_member m ON m.fuzzy_cluster_accession = c.accession
                    WHERE c.fuzzy_run_id = ?
                    """,
                    (run_id,),
                ):
                    mapping[row["form_id"]] = row["public_id"]
            for form_id, cluster_id in mapping.items():
                self.connection.execute(
                    "UPDATE observation SET fuzzy_cluster_id = ? WHERE form_id = ?",
                    (cluster_id, form_id),
                )
        return mapping

    def validate(self) -> list[str]:
        issues: list[str] = []
        self.check_schema()
        integrity = self.connection.execute("PRAGMA integrity_check").fetchall()
        issues.extend(row[0] for row in integrity if row[0] != "ok")
        for row in self.connection.execute("PRAGMA foreign_key_check"):
            issues.append(
                f"foreign-key violation table={row[0]} rowid={row[1]} parent={row[2]}"
            )
        for row in self.connection.execute(
            "SELECT public_id, family, public_digest, full_digest, canonical_json FROM structural_object"
        ):
            computed = hashlib.sha256(row["canonical_json"].encode("utf-8")).hexdigest()
            expected_id = f"txid:{row['family']}.{computed[:PUBLIC_DIGEST_LENGTH]}"
            if computed != row["full_digest"]:
                issues.append(f"full digest mismatch for {row['public_id']}")
            if row["public_digest"] != computed[:PUBLIC_DIGEST_LENGTH] or row["public_id"] != expected_id:
                issues.append(f"public digest mismatch for {row['public_id']}")
        return issues

    def inspect(self, public_id: str) -> dict[str, Any] | None:
        structural = self.connection.execute(
            "SELECT * FROM structural_object WHERE public_id = ?", (public_id,)
        ).fetchone()
        observations = self.connection.execute(
            """
            SELECT o.*, i.sample, i.upstream_tool, a.name AS annotation_name
            FROM observation o
            JOIN import_manifest i ON i.id = o.import_id
            JOIN annotation_context a ON a.id = i.annotation_id
            WHERE o.form_id = ? OR o.splice_chain_id = ? OR o.output_transcript_id = ?
               OR o.output_gene_id = ? OR o.fuzzy_cluster_id = ?
            ORDER BY i.sample, i.upstream_tool, o.original_transcript_id
            """,
            (public_id, public_id, public_id, public_id, public_id),
        ).fetchall()
        locus = self.connection.execute(
            "SELECT * FROM locus WHERE public_id = ?", (public_id,)
        ).fetchone()
        fuzzy = self.connection.execute(
            """
            SELECT c.*, r.algorithm, r.splice_tolerance, r.end_tolerance,
                   r.input_fingerprint
            FROM fuzzy_cluster c
            JOIN fuzzy_run r ON r.id = c.fuzzy_run_id
            WHERE c.public_id = ?
            """,
            (public_id,),
        ).fetchone()
        fuzzy_members = [] if fuzzy is None else [
            row["form_id"]
            for row in self.connection.execute(
                """
                SELECT m.form_id FROM fuzzy_member m
                JOIN fuzzy_cluster c ON c.accession = m.fuzzy_cluster_accession
                WHERE c.public_id = ? ORDER BY m.form_id
                """,
                (public_id,),
            )
        ]
        if structural is None and not observations and locus is None and fuzzy is None:
            return None
        return {
            "query": public_id,
            "structural_object": None if structural is None else dict(structural),
            "locus": None if locus is None else dict(locus),
            "fuzzy_cluster": None
            if fuzzy is None
            else {**dict(fuzzy), "members": fuzzy_members},
            "observations": [dict(row) for row in observations],
        }

    def catalog_rows(self) -> list[sqlite3.Row]:
        return self.connection.execute(
            """
            WITH grouped AS (
                SELECT
                    o.form_id, o.splice_chain_id, o.output_gene_id,
                    COUNT(*) AS observation_count,
                    COUNT(DISTINCT i.sample) AS sample_count,
                    COUNT(DISTINCT i.upstream_tool) AS tool_count
                FROM observation o
                JOIN import_manifest i ON i.id = o.import_id
                GROUP BY o.form_id, o.splice_chain_id, o.output_gene_id
            ), classes AS (
                SELECT form_id, splice_chain_id, output_gene_id,
                       GROUP_CONCAT(classification) AS classifications
                FROM (
                    SELECT DISTINCT form_id, splice_chain_id, output_gene_id, classification
                    FROM observation
                    ORDER BY form_id, splice_chain_id, output_gene_id, classification
                )
                GROUP BY form_id, splice_chain_id, output_gene_id
            )
            SELECT grouped.*, classes.classifications
            FROM grouped JOIN classes
              ON classes.form_id = grouped.form_id
             AND classes.splice_chain_id IS grouped.splice_chain_id
             AND classes.output_gene_id = grouped.output_gene_id
            ORDER BY grouped.output_gene_id, grouped.form_id
            """
        ).fetchall()

    def summary(self) -> dict[str, Any]:
        classification_counts = {
            row["classification"]: row["count"]
            for row in self.connection.execute(
                "SELECT classification, COUNT(*) AS count FROM observation GROUP BY classification"
            )
        }
        return {
            "assembly": self.reference().name,
            "assembly_fingerprint": self.reference().fingerprint,
            "annotations": self.connection.execute("SELECT COUNT(*) FROM annotation_context").fetchone()[0],
            "imports": self.connection.execute("SELECT COUNT(*) FROM import_manifest").fetchone()[0],
            "observations": self.connection.execute("SELECT COUNT(*) FROM observation").fetchone()[0],
            "exact_forms": self.connection.execute(
                "SELECT COUNT(DISTINCT form_id) FROM observation"
            ).fetchone()[0],
            "splice_chains": self.connection.execute(
                "SELECT COUNT(DISTINCT splice_chain_id) FROM observation WHERE splice_chain_id IS NOT NULL"
            ).fetchone()[0],
            "loci": self.connection.execute("SELECT COUNT(*) FROM locus").fetchone()[0],
            "classifications": dict(sorted(classification_counts.items())),
        }
