#!/usr/bin/env python3
"""Build presentation artifacts from retained TxID reports; never rerun benchmarks."""
from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import re
import sys
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


PACKAGE = Path(__file__).resolve().parents[1]
REPO = PACKAGE.parents[1]
RESULTS = REPO / "benchmarks/results/encode-gtf-full"


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def tsv(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def md_table(rows):
    keys = list(rows[0])
    return "\n".join(["| " + " | ".join(keys) + " |", "| " + " | ".join(["---"] * len(keys)) + " |"] + ["| " + " | ".join(str(row[k]) for k in keys) + " |" for row in rows])


def generate_tables():
    partitions = read_json(RESULTS / "identity-partitions.json")
    annotations = read_json(RESULTS / "annotation-release-v29-v49.json")
    extraction = read_json(RESULTS / "talon-model-extraction-audit.json")
    manifest = list(csv.DictReader((RESULTS / "encode-files-checked.tsv").open(), delimiter="\t"))
    by_sample = {row["dataset"]: row for row in extraction["datasets"]}
    inputs = []
    for row in manifest:
        entry = by_sample[row["dataset_id"]]
        inputs.append({"Accession": row["public_accession"], "Platform": row["platform"], "Replicate": row["replicate"], "Transcripts": entry["selected_transcripts"], "Exons": entry["selected_exon_lines"], "SHA256": row["sha256"], "URL": row["source_url"]})
    assert sum(row["Transcripts"] for row in inputs) == 538804
    assert sum(row["Exons"] for row in inputs) == 4248994
    tsv(PACKAGE / "tables/table-s1-inputs.tsv", inputs)

    layers = [("TxID exact form", "txid-sorted", "Exact form"), ("gffcompare tracking", "gffcompare-sorted", "Exact form"), ("gffcompare tracking", "gffcompare-splice-sorted", "Complete splice chain"), ("isoSeQL exact ends", "isoseql-exact-sorted", "Exact form"), ("isoSeQL common junction", "isoseql-junction-sorted", "Complete splice chain")]
    comparison = []
    matrix = []
    for label, key, target in layers:
        v = partitions["evaluations"][key]
        for num, den, rate in [("observation_pairs_false_merged", "eligible_same_cluster_pairs", "false_merge_rate"), ("observation_pairs_false_split", "eligible_same_target_pairs", "false_split_rate")]:
            assert abs(v[num] / v[den] - v[rate]) < 1e-12
        comparison.append({"Layer": label, "Target": target, "Emitted": v["observations_emitted"], "Missing": v["observations_missing"], "Groups": v["identity_clusters"], "FM pairs": v["observation_pairs_false_merged"], "FM denominator": v["eligible_same_cluster_pairs"], "FM rate": f'{v["false_merge_rate"]:.7f}', "FS pairs": v["observation_pairs_false_split"], "FS denominator": v["eligible_same_target_pairs"], "FS rate": f'{v["false_split_rate"]:.7f}'})
        if key != "gffcompare-splice-sorted":
            assert abs(1 - v["matrix_occupied_cells"] / (v["matrix_rows"] * v["matrix_columns"]) - v["matrix_sparsity"]) < 1e-12
            matrix.append({"Layer": label, "Rows": v["matrix_rows"], "Columns": v["matrix_columns"], "Occupied cells": v["matrix_occupied_cells"], "Sparsity (%)": f'{100*v["matrix_sparsity"]:.3f}'})
    tsv(PACKAGE / "tables/table-s2-structural-comparison.tsv", comparison)
    tsv(PACKAGE / "tables/table-s4-matrix-dimensions.tsv", matrix)

    order_rows = []
    order_specs = [("TxID exact form", partitions["order_comparisons"]["txid"]["exact_form_partition"]), ("TxID gene/locus", partitions["order_comparisons"]["txid"]["gene_locus_partition"]), ("gffcompare tracking", partitions["order_comparisons"]["gffcompare_exact"]), ("isoSeQL exact ends", partitions["order_comparisons"]["isoseql_exact_ends"]), ("isoSeQL common junction", partitions["order_comparisons"]["isoseql_common_junction"])]
    for label, v in order_specs:
        order_rows.append({"Layer": label, "Shared observations": v["members_shared"], "Changed observation labels": v["cluster_identifier_changes"], "Changed pair relations": v["pair_relation_changes"]})
    tsv(PACKAGE / "tables/table-s3-order-comparison.tsv", order_rows)
    totals = annotations["totals"]
    assert totals["classification_changes"] == sum(x["classification_changes"] for x in annotations["per_sample"])
    source_data = {"empirical_source": "retained fixed-model evaluation, TxID 0.1.0", "annotation": totals, "matrix": matrix, "order": order_rows, "schematic": {"type": "illustrative, not empirical", "context": "same reference collection, contig and positive strand", "A": [[100,200],[300,400]], "B": [[90,200],[300,420]], "C": [[100,201],[300,400]]}}
    write_json(PACKAGE / "figures/source-data.json", source_data)
    fig_rows = [{"Panel": "b", "Metric": name, "Value": totals[key], "Unit": "observations", "Source": "annotation-release-v29-v49.json/totals/" + key} for name,key in [("Exact form ID changes", "txid_form_changes"), ("Classification changes", "classification_changes"), ("Compared observations", "shared_observations")]]
    for transition, count in sorted(totals["classification_transitions"].items()):
        fig_rows.append({"Panel": "b", "Metric": transition, "Value": count,
                         "Unit": "changed observations", "Source": "annotation-release-v29-v49.json/totals/classification_transitions/" + transition})
    for row in order_rows:
        for key in ["Changed observation labels", "Changed pair relations"]:
            fig_rows.append({"Panel": "Table S3", "Metric": row["Layer"] + ": " + key, "Value": row[key], "Unit": "observations" if "labels" in key else "pairs", "Source": "identity-partitions.json/order_comparisons"})
    for row in matrix:
        fig_rows.append({"Panel": "c", "Metric": row["Layer"], "Value": row["Columns"], "Unit": "columns", "Source": "identity-partitions.json/evaluations"})
    tsv(PACKAGE / "figures/source-data.tsv", fig_rows)
    categories = ["known", "novel_in_known_gene", "new_locus", "ambiguous_gene"]
    transition_rows = []
    for old in categories:
        for new in categories:
            transition_rows.append({"GENCODE v29": old, "GENCODE v49": new,
                "Changed observations": "" if old == new else totals["classification_transitions"].get(f"{old} -> {new}", 0),
                "Display": "diagonal not shown" if old == new else "off-diagonal count"})
    tsv(PACKAGE / "figures/figure1-annotation-transitions.tsv", transition_rows)

    inputs_display = [{k:v for k,v in row.items() if k not in {"SHA256","URL"}} for row in inputs]
    compact = [{"Layer / target": row["Layer"] + " / " + row["Target"], "Emitted / input": str(row["Emitted"]) + " / 538804", "Groups": row["Groups"], "False-merge pairs / eligible pairs": f'{row["FM pairs"]} / {row["FM denominator"]}', "False-split pairs / eligible pairs": f'{row["FS pairs"]} / {row["FS denominator"]}'} for row in comparison]
    tables_text = "\n\n".join(["# Supplementary tables", "## Table S1. Public fixed-model inputs", md_table(inputs_display), "All inputs are WTC11 TALON annotations. PacBio replicates belong to ENCSR309IKK; ONT replicates to ENCSR392BGY. Full download URLs and SHA-256 checksums are in table-s1-inputs.tsv.", "## Table S2. Structural pair comparisons", md_table(compact), "Rates and separate denominators are included in table-s2-structural-comparison.tsv. Targets are coordinate-based comparisons using the shared parser. Missing observations are excluded from eligible pairs and explicitly reported. The isoSeQL comparison used the compatibility patch described in Supplementary Section S4. Different equivalence relations are not a general performance ranking.", "## Table S3. Input-order effects", md_table(order_rows), "Changed observation labels and changed pair relations have different units. The TxID locus layer is separate from exact identity. isoSeQL label changes coexist with unchanged partitions.", "## Table S4. Sample-by-group dimensions", md_table(matrix), "Each matrix has six sample rows. Occupancy represents model observations, not RNA abundance. These layers define different kinds of columns."])
    (PACKAGE / "supplementary-tables.md").write_text(tables_text + "\n", encoding="utf-8")
    return source_data


def figure(data):
    source = PACKAGE / "reproducibility/render_figure1_editorial.py"
    spec = importlib.util.spec_from_file_location("txid_publication_figure", source)
    assert spec and spec.loader
    renderer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(renderer)
    renderer.render(data, PACKAGE / "figures")


def load_docx_builder():
    source=REPO/"deliverables/TxID_manuscript_cohort_revision_2026-08-08/reproducibility/create_manuscript_docx.py"
    spec=importlib.util.spec_from_file_location("txid_note_docx_builder",source)
    assert spec and spec.loader
    module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module)
    module.ZIP_TIME=(2026,9,20,0,0,0)
    return module


def docx(source):
    builder=load_docx_builder()
    lines=source.read_text(encoding="utf-8").splitlines()
    # Convert Markdown links to readable bibliographic text while retaining images.
    lines=[re.sub(r"(?<!!)\[([^\]]+)\]\(([^)]+)\)",r"\1 (\2)",line) for line in lines]
    assets,payloads=builder.collect_assets(lines,source_path=source)
    document=builder.document_xml(lines,assets)
    document=document.replace('w:line="270"','w:line="480"')
    document=document.replace('<w:docGrid ', '<w:lnNumType w:countBy="1" w:restart="continuous"/><w:docGrid ')
    styles=builder.STYLES.replace('w:val="20"','w:val="24"').replace('w:val="23"','w:val="24"')
    core=builder.CORE.replace('TxID GigaScience Technical Note — cohort-scale identity revision','TxID Bioinformatics Application Note draft').replace('Technical Note submission manuscript','Application Note review manuscript').replace('Bilingual-ready Technical Note manuscript with three code-generated scientific figures and three result tables.','English scientific draft; author details and submission declarations are maintained separately.').replace('2026-08-08','2026-09-20')
    if source.name == "manuscript.md" and (PACKAGE / "translations").is_dir():
        core = core.replace('English scientific draft;', 'English–Chinese paragraph-aligned review draft;')
    output=source.with_suffix(".docx")
    with zipfile.ZipFile(output,"w") as archive:
        entries={"[Content_Types].xml":builder.CONTENT_TYPES,"_rels/.rels":builder.ROOT_RELS,"docProps/core.xml":core,"docProps/app.xml":builder.APP,"word/document.xml":document,"word/styles.xml":styles,"word/settings.xml":builder.SETTINGS,"word/_rels/document.xml.rels":builder.document_relationships(assets)}
        for name,payload in entries.items():
            ET.fromstring(payload)
            builder.write_member(archive,name,payload)
        for name,payload in payloads.items(): builder.write_member(archive,"word/media/"+name,payload)
    with zipfile.ZipFile(output) as archive: assert archive.testzip() is None


def bilingual_text(english):
    """Insert reviewed translations without altering the English source blocks."""
    paths = sorted((PACKAGE / "translations").glob("*.zh-CN.json"))
    if not paths:
        return english
    translations = {}
    for path in paths:
        for row in read_json(path):
            source, translated = row["source"], row["translation"]
            if source in translations or not translated.strip() or "\n\n" in translated:
                raise ValueError(f"Duplicate or invalid translation block in {path}")
            translations[source] = translated
    blocks = english.strip().split("\n\n")
    expected = {block for block in blocks if not block.startswith("![")}
    if set(translations) != expected:
        raise ValueError("Translations do not match the English components; update the corresponding translation blocks before rebuilding.")
    paired = []
    for block in blocks:
        paired.append(block)
        if block in translations:
            paired.append(translations[block])
    return "\n\n".join(paired) + "\n"


def supplementary_components():
    """Return reviewed supplement sources, including optional new evidence."""
    paths = [PACKAGE / name for name in [
        "supplementary-methods.md",
        "supplementary-tables.md",
        "supplementary-tool-semantics.md",
    ]]
    for path in paths:
        if not path.is_file():
            raise ValueError(f"Supplementary component is missing: {path.name}")
    interop = PACKAGE / "supplementary-real-interop.md"
    if interop.exists():
        if not interop.is_file():
            raise ValueError("Supplementary interoperation component must be a file")
        paths.append(interop)
    return paths


def manuscript_components():
    """Keep optional declarations after the three body chapters."""
    chapter_names = ["01_introduction.md", "02_implementation.md", "03_evaluation_discussion.md"]
    if sorted(path.name for path in (PACKAGE / "chapters").glob("*.md")) != chapter_names:
        raise ValueError("Manuscript must contain exactly the three planned body chapters")
    paths = [PACKAGE / "front-matter.md"] + [PACKAGE / "chapters" / name for name in chapter_names]
    end_matter = PACKAGE / "end-matter.md"
    if end_matter.exists():
        paths.append(end_matter)
    paths += [PACKAGE / "figure-caption.md", PACKAGE / "references.md"]
    if not all(path.is_file() for path in paths):
        raise ValueError("Manuscript components are incomplete")
    return paths


def assemble():
    paths = manuscript_components()
    english = "\n\n".join(path.read_text().strip() for path in paths)+"\n"
    main = bilingual_text(english)
    if main != english:
        (PACKAGE/"manuscript.en.md").write_text(english,encoding="utf-8")
        docx(PACKAGE/"manuscript.en.md")
    (PACKAGE/"manuscript.md").write_text(main,encoding="utf-8")
    supplement="\n\n".join(path.read_text(encoding="utf-8").strip() for path in supplementary_components())+"\n"
    (PACKAGE/"supplementary.md").write_text(supplement,encoding="utf-8")
    for name in ["manuscript.md","supplementary.md"]: docx(PACKAGE/name)


def source_manifest():
    paths=[RESULTS/name for name in ["identity-partitions.json","invariance-summary.json","annotation-release-v29-v49.json","talon-model-extraction-audit.json","identity-view-audit.json","gencode-v49-identity-view-audit.json","encode-files-checked.tsv","isoseql-compat-audit.json","txid-sorted-recovery-audit.json","txid-shuffled-recovery-audit.json","v49-pacbio-r3-classification-only-audit.json"]]
    paths += [REPO/name for name in ["benchmarks/results/simulation-summary.json","docs/benchmark/encode-gtf-full-2026-07-31.md","docs/benchmark/publication-smoke-2026-07-31.md","docs/superpowers/specs/2026-07-19-novel-transcript-identity-registry-design.md","workflows/publication_benchmark/scripts/evaluate_identity_partitions.py","deliverables/TxID_manuscript_cohort_revision_2026-08-08/reproducibility/create_manuscript_docx.py"]]
    paths += [PACKAGE/name for name in ["supplementary-tool-semantics.md","tables/table-s5-identity-semantics.tsv"]]
    for name in ["supplementary-real-interop.md", "end-matter.md", "translations/end-matter.zh-CN.json", "reproducibility/software-test-report.json"]:
        path = PACKAGE / name
        if path.is_file():
            paths.append(path)
    for name in ["tables/table-s6-real-interop.tsv", "reproducibility/public-access-check.json"]:
        path = PACKAGE / name
        if path.is_file():
            paths.append(path)
    compact = PACKAGE / "reproducibility/real-interop"
    if compact.is_dir():
        paths.extend(sorted(path for path in compact.rglob('*') if path.is_file()))
    independent_report = REPO / 'benchmarks/results/independent-identity-2026-09-23.json'
    if independent_report.is_file():
        paths.append(independent_report)
    write_json(PACKAGE/"reproducibility/source-manifest.json",[{"path":str(path.relative_to(REPO)),"bytes":path.stat().st_size,"sha256":hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths])


if __name__=="__main__":
    source_manifest()
    data=generate_tables()
    figure(data)
    assemble()
    print("Built four numeric source tables, included identity-semantics Table S5 and any interoperation supplement, figure PNG/SVG/PDF/TIFF, main manuscript and supplement DOCX.")
