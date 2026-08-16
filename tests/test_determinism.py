from __future__ import annotations

import concurrent.futures
import csv
import tempfile
import unittest
from pathlib import Path

from txid.identity import identify
from txid.models import Exon, TranscriptModel
from txid.parser import annotation_fingerprint, parse_annotation_lines
from txid.reference import read_fasta, resolve_and_validate
from txid.registry import Registry
from txid.workflow import batch_import, import_file
from txid.writers import write_catalog

from tests.helpers import CALLER_GTF, make_registry, write_inputs


class DeterminismTests(unittest.TestCase):
    def test_annotation_record_order_and_thread_count_invariance(self):
        lines = [line + "\n" for line in CALLER_GTF.splitlines()]
        forward = parse_annotation_lines(lines, fmt="gtf")
        reverse = parse_annotation_lines(reversed(lines), fmt="gtf")
        self.assertEqual(annotation_fingerprint(forward), annotation_fingerprint(reverse))
        assembly = "sha256:" + "f" * 64

        def compute(model):
            return identify(model, assembly).form.public_id

        serial = [compute(model) for model in forward]
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
            threaded = list(pool.map(compute, reversed(forward)))
        self.assertEqual(sorted(serial), sorted(threaded))

    def test_explicit_contig_alias_and_fasta_normalization(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            first = root / "a.fa"
            second = root / "b.fa"
            aliases = root / "aliases.tsv"
            first.write_text(">chr1\nacgt\nACGT\n", encoding="utf-8")
            second.write_text(">chr1 description\nACGTACGT\n", encoding="utf-8")
            aliases.write_text("alias\tprimary\n1\tchr1\n", encoding="utf-8")
            a = read_fasta(first, assembly_name="a", alias_path=aliases)
            b = read_fasta(second, assembly_name="b")
            self.assertEqual(a.fingerprint, b.fingerprint)
            self.assertEqual(
                a.fingerprint,
                "sha256:6dd937f58b7b2aaf135aa2ee50611d8a83a8c290089da09c9ec416a18d6d3114",
            )
            self.assertEqual(
                a.contigs[0].sequence_digest,
                "b28b7e7e6b70661dfee15d5290c4bca097ca145f721c4fbc4de73ad1d1660b8b",
            )
            model = TranscriptModel("1", "+", (Exon(1, 4),), "t")
            self.assertEqual(resolve_and_validate(model, a).contig, "chr1")

    def test_batch_manifest_order_and_catalog_serialization(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            db1_root = root / "one"
            db2_root = root / "two"
            db1_root.mkdir()
            db2_root.mkdir()
            db1, caller1 = make_registry(db1_root)
            db2, caller2 = make_registry(db2_root)
            caller1b = db1_root / "caller_b.gtf"
            caller2b = db2_root / "caller_b.gtf"
            caller1b.write_text(CALLER_GTF.replace("caller", "other"), encoding="utf-8")
            caller2b.write_text(CALLER_GTF.replace("caller", "other"), encoding="utf-8")
            manifest1 = db1_root / "manifest.tsv"
            manifest2 = db2_root / "manifest.tsv"
            header = "input\tsample\ttool\tannotation_name\n"
            rows1 = [
                f"{caller1}\ts2\ttoolB\tref-v1\n",
                f"{caller1b}\ts1\ttoolA\tref-v1\n",
            ]
            rows2 = [
                f"{caller2b}\ts1\ttoolA\tref-v1\n",
                f"{caller2}\ts2\ttoolB\tref-v1\n",
            ]
            manifest1.write_text(header + "".join(rows1), encoding="utf-8")
            manifest2.write_text(header + "".join(rows2), encoding="utf-8")
            results1 = batch_import(
                db1, manifest=manifest1, output_dir=db1_root / "out"
            )
            results2 = batch_import(
                db2,
                manifest=manifest2,
                output_dir=db2_root / "out",
            )

            outputs1 = {
                (row["sample"], row["tool"]): (
                    Path(row["output_gtf"]).read_bytes(),
                    Path(row["mapping"]).read_bytes(),
                )
                for row in results1
            }
            outputs2 = {
                (row["sample"], row["tool"]): (
                    Path(row["output_gtf"]).read_bytes(),
                    Path(row["mapping"]).read_bytes(),
                )
                for row in results2
            }
            self.assertEqual(outputs1, outputs2)

            catalogs = []
            for database, directory in [(db1, db1_root), (db2, db2_root)]:
                target = directory / "catalog.tsv"
                with Registry(database, read_only=True) as registry:
                    rows = [
                        {
                            "form_id": row["form_id"],
                            "splice_chain_id": row["splice_chain_id"] or "",
                            "gene_id": row["output_gene_id"],
                            "classifications": row["classifications"],
                            "observation_count": row["observation_count"],
                            "sample_count": row["sample_count"],
                            "tool_count": row["tool_count"],
                        }
                        for row in registry.catalog_rows()
                    ]
                write_catalog(target, rows)
                catalogs.append(target.read_bytes())
            self.assertEqual(catalogs[0], catalogs[1])


if __name__ == "__main__":
    unittest.main()
