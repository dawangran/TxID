"""Non-interactive TxID command-line interface."""

from __future__ import annotations

import argparse
import json
import sqlite3
import sys
from pathlib import Path

from . import __version__
from .errors import TxIDError
from .plotting import write_gene_svg, write_registry_svg
from .registry import Registry
from .workflow import (
    add_annotation_context,
    batch_import,
    import_file,
    initialize_registry,
    multi_import,
)
from .writers import write_catalog


def _format_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--format", choices=["gtf", "gff3"], help="override format detection")


def _fuzzy_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--fuzzy-splice-tolerance", type=int, metavar="BP")
    parser.add_argument("--fuzzy-end-tolerance", type=int, metavar="BP")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="txid",
        description="Reference-aware deterministic identity registry for transcript models",
    )
    parser.add_argument("--version", action="version", version=f"txid {__version__}")
    commands = parser.add_subparsers(dest="command", required=True)

    init = commands.add_parser("init", help="initialize a registry from FASTA and reference annotation")
    init.add_argument("--db", required=True, type=Path)
    init.add_argument("--fasta", required=True, type=Path)
    init.add_argument("--assembly", required=True)
    init.add_argument("--annotation", required=True, type=Path)
    init.add_argument("--annotation-name", required=True)
    init.add_argument("--contig-aliases", type=Path)
    _format_option(init)

    annotation_add = commands.add_parser(
        "annotation-add", help="register another annotation release on the same assembly"
    )
    annotation_add.add_argument("--db", required=True, type=Path)
    annotation_add.add_argument("--annotation", required=True, type=Path)
    annotation_add.add_argument("--annotation-name", required=True)
    _format_option(annotation_add)

    add = commands.add_parser("add", help="atomically import one caller annotation")
    add.add_argument("--db", required=True, type=Path)
    add.add_argument("--input", required=True, type=Path)
    add.add_argument("--sample", required=True)
    add.add_argument("--tool", required=True)
    add.add_argument("--annotation-name", required=True)
    add.add_argument("--output-gtf", required=True, type=Path)
    add.add_argument("--mapping", required=True, type=Path)
    _format_option(add)
    _fuzzy_options(add)

    batch = commands.add_parser("batch", help="import a tab-delimited manifest")
    batch.add_argument("--db", required=True, type=Path)
    batch.add_argument("--manifest", required=True, type=Path)
    batch.add_argument("--output-dir", required=True, type=Path)
    _fuzzy_options(batch)

    multi_add = commands.add_parser(
        "multi-add",
        help="import multiple annotations with shared provenance without a manifest",
    )
    multi_add.add_argument("--db", required=True, type=Path)
    multi_add.add_argument("--input", required=True, nargs="+", type=Path, metavar="PATH")
    multi_add.add_argument(
        "--samples",
        nargs="+",
        metavar="SAMPLE",
        help="one sample name per input; defaults to input basenames",
    )
    multi_add.add_argument("--tool", required=True)
    multi_add.add_argument("--annotation-name", required=True)
    multi_add.add_argument("--output-dir", required=True, type=Path)
    _format_option(multi_add)
    _fuzzy_options(multi_add)

    export = commands.add_parser("export", help="write a deterministic cohort catalog")
    export.add_argument("--db", required=True, type=Path)
    export.add_argument("--catalog", required=True, type=Path)

    validate = commands.add_parser("validate", help="validate schema, foreign keys, and stored digests")
    validate.add_argument("--db", required=True, type=Path)

    inspect = commands.add_parser("inspect", help="inspect an exact, locus, or observed identifier")
    inspect.add_argument("--db", required=True, type=Path)
    inspect.add_argument("identifier")

    plot = commands.add_parser(
        "plot", help="render an SVG registry overview or one gene's transcript structures"
    )
    plot.add_argument("--db", required=True, type=Path)
    plot.add_argument("--output", required=True, type=Path)
    plot.add_argument(
        "--gene",
        help="gene/locus ID, unique reference gene_name, or upstream gene ID",
    )
    return parser


def _json_stdout(value: object) -> None:
    json.dump(value, sys.stdout, sort_keys=True, indent=2, ensure_ascii=False)
    sys.stdout.write("\n")


def run(args: argparse.Namespace) -> int:
    if args.command == "init":
        _json_stdout(
            initialize_registry(
                args.db,
                fasta=args.fasta,
                assembly_name=args.assembly,
                annotation=args.annotation,
                annotation_name=args.annotation_name,
                alias_path=args.contig_aliases,
                annotation_format=args.format,
            )
        )
    elif args.command == "annotation-add":
        _json_stdout(
            add_annotation_context(
                args.db,
                annotation=args.annotation,
                annotation_name=args.annotation_name,
                annotation_format=args.format,
            )
        )
    elif args.command == "add":
        result = import_file(
            args.db,
            input_path=args.input,
            sample=args.sample,
            tool=args.tool,
            annotation_name=args.annotation_name,
            output_gtf=args.output_gtf,
            mapping_path=args.mapping,
            annotation_format=args.format,
            fuzzy_splice_tolerance=args.fuzzy_splice_tolerance,
            fuzzy_end_tolerance=args.fuzzy_end_tolerance,
        )
        _json_stdout(
            {
                "import_id": result.import_id,
                "created": result.created,
                "input_checksum": result.input_checksum,
                "transcripts": len(result.assignments),
                "output_gtf": str(args.output_gtf),
                "mapping": str(args.mapping),
            }
        )
    elif args.command == "batch":
        _json_stdout(
            {
                "imports": batch_import(
                    args.db,
                    manifest=args.manifest,
                    output_dir=args.output_dir,
                    fuzzy_splice_tolerance=args.fuzzy_splice_tolerance,
                    fuzzy_end_tolerance=args.fuzzy_end_tolerance,
                )
            }
        )
    elif args.command == "multi-add":
        _json_stdout(
            {
                "imports": multi_import(
                    args.db,
                    inputs=args.input,
                    samples=args.samples,
                    tool=args.tool,
                    annotation_name=args.annotation_name,
                    output_dir=args.output_dir,
                    annotation_format=args.format,
                    fuzzy_splice_tolerance=args.fuzzy_splice_tolerance,
                    fuzzy_end_tolerance=args.fuzzy_end_tolerance,
                )
            }
        )
    elif args.command == "export":
        with Registry(args.db, read_only=True) as registry:
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
        write_catalog(args.catalog, rows)
        _json_stdout({"catalog": str(args.catalog), "rows": len(rows)})
    elif args.command == "validate":
        with Registry(args.db, read_only=True) as registry:
            issues = registry.validate()
            summary = registry.summary()
        _json_stdout({"valid": not issues, "issues": issues, "summary": summary})
        return 0 if not issues else 1
    elif args.command == "inspect":
        with Registry(args.db, read_only=True) as registry:
            result = registry.inspect(args.identifier)
        if result is None:
            raise TxIDError(f"identifier not found: {args.identifier}")
        _json_stdout(result)
    elif args.command == "plot":
        with Registry(args.db, read_only=True) as registry:
            if args.gene is None:
                write_registry_svg(args.output, registry)
                result = {"plot": str(args.output)}
            else:
                plot_summary = write_gene_svg(args.output, registry, args.gene)
                result = {
                    "plot": str(args.output),
                    **plot_summary,
                }
        _json_stdout(result)
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    try:
        return run(parser.parse_args(argv))
    except (TxIDError, OSError, sqlite3.Error, ValueError) as error:
        print(f"txid: error: {error}", file=sys.stderr)
        return 2
