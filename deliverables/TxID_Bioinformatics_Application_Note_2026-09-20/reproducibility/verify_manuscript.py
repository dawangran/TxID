#!/usr/bin/env python3
"""Validate the draft package against retained evidence and assembled source text."""
from __future__ import annotations
import csv
import hashlib
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path
from xml.etree import ElementTree as ET
from PIL import Image

PACKAGE=Path(__file__).resolve().parents[1]
REPO=PACKAGE.parents[1]
checks=[]
EXPECTED_REFERENCES = {
    "ENCODE2020": "10.1038/s41586-020-2493-4",
    "Jackman2015": "10.1371/journal.pone.0128026",
    "Kabza2024": "10.1038/s41467-024-51584-3",
    "Kovaka2019": "10.1186/s13059-019-1910-1",
    "Liu2026": "10.1093/bioinformatics/btaf680",
    "Mudge2025": "10.1093/nar/gkae1078",
    "PardoPalacios2024": "10.1038/s41592-024-02229-2",
    "Pertea2020": "10.12688/f1000research.23297.2",
    "Prjibelski2023": "10.1038/s41587-022-01565-y",
    "Wagner2021": "10.1016/j.xgen.2021.100027",
    "Wyman2020": "10.1101/672931",
}
AUTHOR_FIELDS = {
    "front-matter.md": ["[AUTHOR NAMES AND AFFILIATION INDICES]",
                        "[AFFILIATIONS: DEPARTMENT, INSTITUTION, CITY, POSTCODE, COUNTRY]",
                        "[CORRESPONDING AUTHOR NAME]", "[CONTACT EMAIL]"],
    "end-matter.md": ["[AUTHOR FUNDING STATEMENT]", "[AUTHOR CONTRIBUTIONS]",
                      "[AUTHOR REVIEW STATEMENT]", "[AUTHOR COMPETING-INTEREST STATEMENT]"],
}

def require(condition, message):
    if not condition: raise AssertionError(message)
    checks.append(message)

def plain(text):
    text=re.sub(r"!\[[^\]]*\]\([^)]*\)","",text)
    text=re.sub(r"\[([^\]]+)\]\([^)]*\)",r"\1",text)
    return text.replace("`","").replace("#","")

def words(text): return len(plain(text).split())

def software_test_report():
    path = PACKAGE / "reproducibility/software-test-report.json"
    if not path.is_file():
        return {"status": "NOT RUN", "reason": "No retained software-test-report.json is available"}
    report = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(report, dict) or not report:
        raise ValueError("Software test report must contain a nonempty JSON object")
    return report

def public_access_report():
    path = PACKAGE / "reproducibility/public-access-check-2026-09-27.json"
    if not path.is_file():
        path = PACKAGE / "reproducibility/public-access-check.json"
    if not path.is_file():
        return {"status": "NOT CHECKED"}, "Public repository availability has no retained access check"
    report = json.loads(path.read_text())
    require(report.get("authentication") == "none" and isinstance(report.get("web_http_status"), int),
            "Public access report records an anonymous browser HTTP check")
    status = report["web_http_status"]
    if status != 200:
        limitation = (f"Anonymous repository access returned HTTP {status} on {report.get('date', 'the recorded date')}; "
                      "a public software link remains unresolved (the response does not establish why)")
    else:
        limitation = "Public repository access is recorded, but persistent software/data archive records require confirmation"
    return report, limitation


def verify_real_interop(manuscript, supplement):
    """Recount the bounded case from retained observation coordinates, without rerunning it."""
    package = PACKAGE / "reproducibility/real-interop"
    evidence = package / "interop-evaluation"
    def read_json(path): return json.loads(path.read_text())
    def read_rows(path):
        with path.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()

    manifest = read_json(package / "evidence-manifest.json")
    require(manifest["schema"] == "txid.real-interop-evidence-package.v1" and manifest["valid"] is True,
            "Real-case compact evidence manifest has the expected schema and valid status")
    copied = manifest["copied_files"]
    original_copies_checked = 0
    source_files_checked = 0
    archive_sources_checked = 0
    source_archives = sorted((PACKAGE / "submission").glob("TxID-*-submission-source.zip"))
    require(len(copied) == manifest["copied_file_count"] == len({r["destination_package_path"] for r in copied}),
            "Real-case evidence copy list has unique paths and the reported count")
    for item in copied:
        destination = package / item["destination_package_path"]
        source = REPO / item["source_repository_path"]
        require(destination.resolve().is_relative_to(package.resolve()) and source.resolve().is_relative_to(REPO.resolve()),
                "Real-case evidence paths stay inside declared roots: " + item["destination_package_path"])
        require(item["verified_equal"] is True and destination.stat().st_size == item["bytes"] and
                sha(destination) == item["copy_sha256"] == item["source_sha256"],
                "Real-case evidence copy matches recorded original SHA-256: " + item["destination_package_path"])
        if source.is_file():
            require(sha(source) == item["source_sha256"], "Available original evidence checksum agrees: " + item["source_repository_path"])
            original_copies_checked += 1
    for item in manifest["source_hash_checks"]:
        source = REPO / item["source_repository_path"]
        require(item["verified_against_original"] is True, "Real-case source verification is recorded: " + item["source_repository_path"])
        if source.is_file():
            require(sha(source) == item["sha256"], "Available real-case audited source checksum agrees: " + item["source_repository_path"])
            source_files_checked += 1
        else:
            found = False
            for archive_path in source_archives:
                with zipfile.ZipFile(archive_path) as archive:
                    matches = [name for name in archive.namelist() if name.endswith("/" + item["source_repository_path"])]
                    if matches:
                        require(len(matches) == 1 and hashlib.sha256(archive.read(matches[0])).hexdigest() == item["sha256"],
                                "Included source archive checksum agrees: " + item["source_repository_path"])
                        found = True
                        break
            archive_sources_checked += int(found)
    for item in manifest["generated_supplemental_files"]:
        path = PACKAGE / item["manuscript_relative_path"]
        require(path.stat().st_size == item["bytes"] and sha(path) == item["sha256"],
                "Real-case derived supplement checksum agrees: " + path.name)

    summary = read_json(evidence / "summary.json")
    for item in manifest["source_hash_checks"]:
        records = [record for record in summary["sources"] if record["path"].endswith("/" + item["source_repository_path"])]
        require(len(records) == 1 and records[0]["sha256"] == item["sha256"],
                "Source hash provenance agrees between immutable summary and package manifest: " + item["source_repository_path"])
    require(summary["schema"] == "txid.real-interop-evaluation.v1" and summary["valid"] is True and
            summary["software_version"] == "txid 0.1.3" and summary["evidence_kind"] == "real-read-case",
            "Real-case result is valid, uses TxID 0.1.3 and is labelled as actual caller evidence")
    require(manifest["original_summary_modified"] is False and sha(evidence / "summary.json") == manifest["original_summary_sha256"],
            "Original real-case summary remains unchanged")
    rows = read_rows(evidence / "per-observation-evidence.tsv")
    require(len(rows) == 311 and len({(r["caller"], r["original_transcript_id"]) for r in rows}) == 311,
            "Real-case evidence contains 311 unique caller observations")
    require(all(r["initial_mapping_matches_oracle"] == r["both_final_registries_match_oracle"] == "True" for r in rows),
            "Every retained observation agrees with the oracle in initial and final registries")
    by_caller = {caller: [r for r in rows if r["caller"] == caller] for caller in ("isoquant", "stringtie")}
    require(sum(map(len, by_caller.values())) == 311, "Real-case evidence contains exactly the two declared callers")
    for row in rows:
        exons = tuple(tuple(x) for x in json.loads(row["exons"]))
        require(exons == tuple(sorted(exons)) and row["strand"] in {"+", "-"} and
                all(0 < start <= end for start, end in exons), "Real-case exon coordinates are ordered and valid: " + row["original_transcript_id"])
        row["form_coordinate"] = (row["contig"], row["strand"], exons)
        introns = tuple((left[1]+1, right[0]-1) for left, right in zip(exons, exons[1:]))
        if row["strand"] == "-": introns = introns[::-1]
        row["chain_coordinate"] = (row["contig"], row["strand"], introns) if introns else None
        require(bool(row["splice_chain_id"]) == bool(introns), "Only multi-exon observations receive SC IDs: " + row["original_transcript_id"])
    sets = {}
    for caller, expected in [("isoquant", (126, 126, 122)), ("stringtie", (185, 185, 176))]:
        observations = by_caller[caller]
        forms = {r["form_id"] for r in observations}
        chains = {r["splice_chain_id"] for r in observations if r["splice_chain_id"]}
        sets[caller] = {"form": forms, "chain": chains}
        counts = summary["caller_counts"][caller]
        require((len(observations), len(forms), len(chains)) == expected ==
                (counts["observations"], counts["distinct_forms"], counts["distinct_multi_exon_chains"]),
                "Caller observation/form/chain totals reproduce the result: " + caller)
    shared = sets["isoquant"]["form"] & sets["stringtie"]["form"]
    require(len(shared) == 110 and len(sets["isoquant"]["form"] | sets["stringtie"]["form"]) == 201,
            "Real-case exact-form intersection and union reproduce 110 and 201")
    paired = [(i, j, a, b) for i, a in enumerate(by_caller["isoquant"]) for j, b in enumerate(by_caller["stringtie"])]
    for layer, field, coordinate, key, count in [
            ("form", "form_id", "form_coordinate", "exact_form_join", 110),
            ("chain", "splice_chain_id", "chain_coordinate", "multi_exon_chain_join", 115)]:
        target = {(i,j) for i,j,a,b in paired if a[coordinate] is not None and a[coordinate] == b[coordinate]}
        joined = {(i,j) for i,j,a,b in paired if a[field] and a[field] == b[field]}
        measured = summary[key]
        require(target == joined and len(target) == count and
                all(measured[k] == count for k in ["correctly_joined_pairs", "equal_coordinate_pairs", "joined_pairs", "false_merge_denominator", "false_split_denominator"]) and
                measured["false_merge_pairs"] == measured["false_split_pairs"] == 0,
                "Direct coordinate comparison reproduces join counts and both error denominators: " + layer)
        require(measured["false_merge_rate"] == measured["false_split_rate"] == 0.0,
                "Nonzero real-case denominators support zero observed join-error rates: " + layer)
        matrix = read_rows(evidence / f"caller-by-{layer}.tsv")
        identifiers = sorted(sets["isoquant"][layer] | sets["stringtie"][layer])
        expected_matrix = [{"caller": caller, **{key: str(int(key in sets[caller][layer])) for key in identifiers}} for caller in ("isoquant", "stringtie")]
        require(matrix == expected_matrix, "Every binary matrix cell matches observation membership: " + layer)
        occupied = sum(len(sets[c][layer]) for c in by_caller)
        overlap = len(sets["isoquant"][layer] & sets["stringtie"][layer])
        m = summary["matrices"][layer]
        require((m["rows"], m["columns"], m["occupied_cells"], m["shared_columns"]) == (2, len(identifiers), occupied, overlap) and
                abs(m["sparsity"] - (1-occupied/(2*len(identifiers)))) < 1e-12,
                "Real-case matrix dimensions, occupancy and sparsity agree: " + layer)
    variants = [(a,b) for _,_,a,b in paired if a["chain_coordinate"] is not None and a["chain_coordinate"] == b["chain_coordinate"] and a["form_coordinate"] != b["form_coordinate"]]
    variant_chains = {a["splice_chain_id"] for a,b in variants}
    require(len(variants) == summary["same_chain_different_form_cross_caller_pairs"] == 7 and
            len(variant_chains) == summary["shared_chains_with_distinct_forms"] == 6 and
            variant_chains == set(summary["end_variant_chain_ids"]),
            "Seven coordinate-defined terminal-variant pairs occur in six shared chains")
    require(len(sets["isoquant"]["chain"] & sets["stringtie"]["chain"]) == 109,
            "109 shared chain groups remain distinct from 115 matched observation pairs")
    native = {r["original_transcript_id"] for r in by_caller["isoquant"]} & {r["original_transcript_id"] for r in by_caller["stringtie"]}
    require(not native and summary["original_id_intersection"] == 0 and summary["coordinate_join_of_original_ids"]["false_merge_rate"] is None,
            "Native caller-name intersection is zero and its zero-denominator false-merge rate is undefined")
    require(summary["independent_initial_registries"] == 2 and summary["import_orders"] == [["isoquant","stringtie"],["stringtie","isoquant"]] and
            summary["reverse_order_exact_observation_changes"] == summary["reverse_order_canonical_or_digest_changes"] == 0,
            "Real-case summary records two initial registries and unchanged opposite-order exact results")
    zero_fields = ["canonical_or_full_digest_mismatches", "existing_canonical_or_digest_changes", "existing_exact_observation_changes", "existing_exact_observation_disappearances", "observation_coordinate_or_identity_mismatches"]
    require(len(summary["stage_audits"]) == 4 and all(stage[k] == 0 for stage in summary["stage_audits"] for k in zero_fields),
            "All four registry stages retain zero recorded structural/key mismatches")

    audit = read_json(package / "supplement-classification-audit.json")
    require(audit["schema"] == "txid.real-interop-classification-audit.v1" and audit["original_summary_modified"] is False and
            audit["original_summary_sha256"] == sha(evidence / "summary.json"), "Classification supplement identifies the unchanged source summary")
    mappings = {}
    for caller in by_caller:
        mapping = read_rows(evidence / f"{caller}-first/01-{caller}.mapping.tsv")
        mappings[caller] = {r["txid_form"]: r for r in mapping}
        require(set(mappings[caller]) == sets[caller]["form"] and dict(Counter(r["classification"] for r in mapping)) == audit["caller_classification_counts"][caller],
                "Classification counts and initial mapping forms agree: " + caller)
    require(audit["shared_exact_form_count"] == audit["shared_forms_with_equal_preserved_reference_transcript_id"] == 110 and
            audit["shared_classification_pairs"] == [{"forms":110,"isoquant":"known","stringtie":"known"}] and
            audit["shared_form_family_counts"] == {"txid:TF1":108,"txid:SE1":2} and
            all(mappings[c][key]["classification"] == "known" for c in mappings for key in shared) and
            all(mappings["isoquant"][key]["txid_transcript_id"] == mappings["stringtie"][key]["txid_transcript_id"] for key in shared),
            "All 110 shared forms are reference-known and preserve the same reference transcript IDs")
    details = audit["shared_form_details"]
    require(len(details) == 110 and {r["form_id"] for r in details} == shared and
            all(r[f"{c}_classification"] == mappings[c][r["form_id"]]["classification"] and
                r[f"{c}_preserved_reference_transcript_id"] == mappings[c][r["form_id"]]["txid_transcript_id"] for r in details for c in mappings),
            "Every classification-audit shared-form detail agrees with the initial mappings")
    for caller, other, expected in [("isoquant","stringtie",16),("stringtie","isoquant",75)]:
        unique = sets[caller]["form"] - sets[other]["form"]
        unmatched = audit["unmatched_exact_forms"][caller]
        other_coordinates = {r["form_coordinate"] for r in by_caller[other]}
        require(len(unique) == unmatched["count"] == expected and set(unmatched["form_ids"]) == unique and
                unmatched["coordinate_counterparts_in_other_caller"] == 0 and
                all(r["form_coordinate"] not in other_coordinates for r in by_caller[caller] if r["form_id"] in unique),
                "Unmatched forms have unequal coordinates rather than missed identity joins: " + caller)

    table = read_rows(PACKAGE / "tables/table-s6-real-interop.tsv")
    expected = [
        ("Input observations","126","185","311 total"),
        ("Exact transcript forms (TF1 + SE1)","126","185","110 shared; 201 union"),
        ("Multi-exon splice chains (SC1)","122","176","109 shared; 189 union"),
        ("Forms without an exact counterpart in the other caller","16","75","91 total"),
        ("Shared chains with distinct transcript forms","NA","NA","6 chains; 7 cross-caller pairs"),
        ("Exact-form join: correctly matched pairs","NA","NA","110/110"),
        ("Exact-form join: false merges; false splits","NA","NA","0/110; 0/110"),
        ("Splice-chain join: correctly matched pairs","NA","NA","115/115"),
        ("Splice-chain join: false merges; false splits","NA","NA","0/115; 0/115"),
        ("Original transcript-ID string matches (descriptive)","126","185","0 shared strings; 0 joined pairs"),
        ("Original transcript-ID join: missed equal forms (descriptive)","NA","NA","110/110"),
        ("Shared exact forms classified known in both callers","NA","NA","110/110"),
        ("Incremental import: changes to existing exact observation keys","0/126","0/185","0 changes"),
        ("Incremental import: changes to existing canonical objects or full digests","0/8673","0/8669","0 changes"),
        ("Reverse-order final comparison: exact observation keys","NA","NA","0/311 changed"),
        ("Reverse-order final comparison: canonical objects or full digests","NA","NA","0/8684 changed"),
    ]
    require([(r["Metric"],r["IsoQuant"],r["StringTie"],r["Cross-caller result"]) for r in table] == expected and
            all(r["Unit and interpretation"].strip() for r in table), "All sixteen Table S6 metric values and denominators agree with the checked source case")
    require([s["structural_objects_checked"] for s in summary["stage_audits"]] == [8673,8684,8669,8684],
            "Table S6 canonical-object denominators agree with the four source audits")
    require("Table S6" in supplement and "110/110 equal-coordinate pairs recovered" in supplement and
            "115/115 equal-chain pairs recovered" in supplement and "201 columns and 311 occupied cells" in supplement and
            "189 columns and 298 occupied cells" in supplement,
            "Assembled S9 includes the checked exact/chain pair results and both matrix occupancies")
    require("reference-known" in manuscript and "shared novel" in manuscript and "reference-alias" in manuscript,
            "Main manuscript retains reference-known overlap and the novel/reference-alias limitations")
    require("All 110 shared exact forms" in supplement and "no example of a shared novel exact form" in supplement and
            "does not show an advantage over joining shared reference aliases" in supplement and
            "does not independently derive the FASTA fingerprint" in supplement,
            "S9 retains known-only evidence, reference-alias limits and the supplied FASTA fingerprint boundary")
    return {"observations":311,"shared_forms":110,"form_union":201,"shared_chains":109,"chain_pairs":115,"end_variant_pairs":7,"shared_forms_all_reference_known":True,
            "evidence_verification_mode": "package-and-original-files" if original_copies_checked == len(copied) else "package-with-recorded-original-hashes",
            "packaged_evidence_files_checked":len(copied), "available_original_evidence_files_checked":original_copies_checked,
            "source_files_checked":source_files_checked,"source_archive_members_checked":archive_sources_checked,
            "source_hashes_only":len(manifest["source_hash_checks"])-source_files_checked-archive_sources_checked}

def main():
    chapter_names=["01_introduction.md","02_implementation.md","03_evaluation_discussion.md"]
    require(sorted(p.name for p in (PACKAGE/"chapters").glob("*.md"))==chapter_names,"Exactly the three planned body chapters exist")
    body_paths=[PACKAGE/"chapters"/name for name in chapter_names]
    components=[PACKAGE/"front-matter.md"]+body_paths
    if (PACKAGE/"end-matter.md").is_file():
        components.append(PACKAGE/"end-matter.md")
    components += [PACKAGE/"figure-caption.md",PACKAGE/"references.md"]
    assembled="\n\n".join(path.read_text().strip() for path in components)+"\n"
    display_manuscript=(PACKAGE/"manuscript.md").read_text()
    manuscript=assembled
    bilingual=bool(list((PACKAGE/"translations").glob("*.zh-CN.json")))
    if bilingual:
        from build_artifacts import bilingual_text
        require(display_manuscript==bilingual_text(assembled),"Every English block is followed by its reviewed Chinese translation")
        require((PACKAGE/"manuscript.en.md").read_text()==assembled,"English manuscript is exactly assembled from reviewed components")
    else:
        require(display_manuscript==assembled,"Main manuscript is exactly assembled from reviewed components")
    for path,minimum in zip(body_paths,[220,380,550]):
        require(words(path.read_text())>=minimum,path.name+" meets the revised internal drafting budget")
    require(words(manuscript)<2200,"Main draft including references/caption remains below 2,200 whitespace-counted words")
    require(manuscript.count("![")==1,"Exactly one composite figure is embedded")
    require(not re.search(r"(?im)^\s*(?:[-*+]\s|\d+\.\s)","\n".join(p.read_text() for p in body_paths)),"Main body is continuous prose without lists")
    require(not re.search(r"(?i)TO BE ADDED|\bTODO\b|\bTBD\b|\[INSERT|待补|待真实",manuscript+(PACKAGE/"supplementary.md").read_text()),"Scientific manuscript and supplement contain no fabricated replacement placeholders")
    author_manifest = json.loads((PACKAGE/"submission/author-fields.json").read_text())
    require(author_manifest["allowed_placeholders"] == AUTHOR_FIELDS,
            "Only the explicitly authorized author information and declaration fields are allowed")
    remaining_author_fields = []
    for path in components:
        fields = re.findall(r"\[[A-Z][A-Z0-9 ,:/-]+\](?!\()", path.read_text())
        require(set(fields) <= set(AUTHOR_FIELDS.get(path.name, [])),
                "No unapproved scientific or availability placeholder appears in " + path.name)
        remaining_author_fields.extend(fields)
    front = (PACKAGE/"front-matter.md").read_text()
    require(re.findall(r"(?m)^### (.+)$", front) ==
            ["Summary", "Availability and Implementation", "Contact", "Supplementary Information"],
            "Application Note abstract has all four required headings in order")
    require("https://github.com/dawangran/TxID" in front and
            "35c124be9b0e76cd4d71b39c3ec394d035e5be8d" in (PACKAGE/"end-matter.md").read_text(),
            "Software location and fixed-revision reproducibility link are supplied without a fabricated archive identifier")
    require("### Funding" in (PACKAGE/"end-matter.md").read_text() and
            "## Conflict of interest" in (PACKAGE/"end-matter.md").read_text(),
            "Funding and competing-interest declarations have designated author completion fields")
    require(not re.search(r"(?i)first and foremost|moreover|furthermore|it is worth noting|groundbreaking|transformative|seamlessly",manuscript),"Stock emphasis/transition scan passed")

    reference_text=(PACKAGE/"references.md").read_text()
    ref_doi_list=re.findall(r"https://doi\.org/([^\s)]+)",reference_text)
    bib_text=(PACKAGE/"references.bib").read_text()
    bib_entries=re.findall(r"(?ms)^@article\{([^,\s]+),\s*(.*?)(?=^@|\Z)",bib_text)
    bib_dois={}
    for key,entry in bib_entries:
        matches=re.findall(r"doi\s*=\s*\{([^}]+)\}",entry)
        require(len(matches)==1 and key not in bib_dois,"BibTeX key is unique and has one DOI: "+key)
        bib_dois[key]=matches[0]
    require(bib_dois==EXPECTED_REFERENCES and len(ref_doi_list)==11 and set(ref_doi_list)==set(EXPECTED_REFERENCES.values()),"Eleven exact reference keys and DOI entries agree, including Kabza2024 and the two caller methods")
    citation_keys=["ENCODE Project Consortium et al., 2020","Jackman et al., 2015","Kabza et al., 2024","Kovaka et al., 2019","Liu and Chun, 2026","Mudge et al., 2025","Pardo-Palacios et al., 2024","Pertea and Pertea, 2020","Prjibelski et al., 2023","Wagner et al., 2021","Wyman et al., 2020"]
    for key in citation_keys: require(key in manuscript,"In-text citation present: "+key)
    require((PACKAGE/"reproducibility/reference-checks.md").is_file(),"Primary-record verification and unresolved access limits are recorded")

    source_manifest=json.loads((PACKAGE/"reproducibility/source-manifest.json").read_text())
    additional_sources=[PACKAGE/name for name in ["supplementary-tool-semantics.md","tables/table-s5-identity-semantics.tsv","tables/table-s6-real-interop.tsv","reproducibility/public-access-check.json","submission/author-fields.json"]]
    if (PACKAGE/"reproducibility/public-access-check-2026-09-27.json").is_file():
        additional_sources.append(PACKAGE/"reproducibility/public-access-check-2026-09-27.json")
    additional_sources += sorted(path for path in (PACKAGE/"reproducibility/real-interop").rglob("*") if path.is_file())
    for name in ["supplementary-real-interop.md", "end-matter.md", "translations/end-matter.zh-CN.json", "reproducibility/software-test-report.json"]:
        if (PACKAGE/name).is_file():
            additional_sources.append(PACKAGE/name)
    manifest_paths=[item["path"] for item in source_manifest]
    require(len(manifest_paths)==len(set(manifest_paths)),"Source manifest has no duplicate paths")
    for path in additional_sources:
        require(str(path.relative_to(REPO)) in manifest_paths,"Source manifest includes current supplementary/declaration/test source: "+path.name)
    for item in source_manifest:
        path=REPO/item["path"]
        require(hashlib.sha256(path.read_bytes()).hexdigest()==item["sha256"],"Source checksum unchanged: "+item["path"])
    reports=REPO/"benchmarks/results/encode-gtf-full"
    partitions=json.loads((reports/"identity-partitions.json").read_text())
    data=json.loads((PACKAGE/"figures/source-data.json").read_text())
    annotations=json.loads((reports/"annotation-release-v29-v49.json").read_text())
    require(data["annotation"]==annotations["totals"],"Figure annotation values equal source JSON")
    svg_tree=ET.parse(PACKAGE/"figures/figure1.svg")
    svg_text=" ".join(svg_tree.getroot().itertext())
    order_sources=[partitions["order_comparisons"]["txid"]["exact_form_partition"],partitions["order_comparisons"]["txid"]["gene_locus_partition"],partitions["order_comparisons"]["gffcompare_exact"],partitions["order_comparisons"]["isoseql_exact_ends"],partitions["order_comparisons"]["isoseql_common_junction"]]
    for row,source in zip(data["order"],order_sources):
        require((row["Shared observations"],row["Changed observation labels"],row["Changed pair relations"])==(source["members_shared"],source["cluster_identifier_changes"],source["pair_relation_changes"]),"Figure order values equal source JSON: "+row["Layer"])
    plotted=json.loads((PACKAGE/"figures/figure1-plot-data.json").read_text())
    svg_groups={node.get("id"):node for node in svg_tree.getroot().iter() if node.get("id")}
    def svg_label(gid): return "".join(svg_groups[gid].itertext()).strip()
    def normalized_label(gid): return " ".join(svg_label(gid).split())
    for i,key in enumerate(["txid_form_changes","classification_changes"]):
        count=annotations["totals"][key]
        total=annotations["totals"]["shared_observations"]
        expected=100*count/total
        actual=plotted["annotation"][i]
        require(actual["changed"]==count and actual["total"]==total and abs(actual["percent"]-expected)<1e-12,"Annotation plotted total and denominator agree with raw report: "+key)
    require(svg_label("b-classification-summary")=="Changed observations only: 29,435 / 538,804 (5.463%)","Alluvial is explicitly limited to changed observations with the correct denominator")
    require(svg_label("b-exact-id-changes")=="Exact form IDs changed: 0","Exact form ID change count is shown separately")
    raw_flows=annotations["totals"]["classification_transitions"]
    require(len(plotted["flows"])==len(raw_flows)==8,"All eight nonzero annotation transitions are drawn")
    left_totals={}; right_totals={}
    for flow in plotted["flows"]:
        old,new=flow["from"],flow["to"]
        count=raw_flows[f"{old} -> {new}"]
        require(flow["count"]==count,"Alluvial transition count matches the original report: "+old+" -> "+new)
        expected_height=count*plotted["flow_scale_mm_per_observation"]
        require(abs(flow["height_mm"]-expected_height)<1e-12 and all(abs(interval[1]-interval[0]-expected_height)<1e-10 for interval in [flow["left_interval"],flow["right_interval"]]),"Alluvial ribbon width is proportional at both ends: "+old+" -> "+new)
        node=svg_groups[f"flow-{old}-to-{new}"]
        path=next(node.iter("{http://www.w3.org/2000/svg}path"))
        numbers=[float(value) for value in re.findall(r"[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?",path.get("d"))]
        require(len(numbers)==16 and abs((numbers[15]-numbers[1])-expected_height*72/25.4)<1e-5 and abs((numbers[9]-numbers[7])-expected_height*72/25.4)<1e-5,"SVG ribbon endpoints reproduce the correct count scale: "+old+" -> "+new)
        left_totals[old]=left_totals.get(old,0)+count
        right_totals[new]=right_totals.get(new,0)+count
    require(left_totals==plotted["left_totals"] and right_totals==plotted["right_totals"] and sum(left_totals.values())==sum(right_totals.values())==29435,"Alluvial conserves all changed observations on both sides")
    names={"novel_in_known_gene":"Novel in gene","new_locus":"New locus","ambiguous_gene":"Ambiguous","known":"Known"}
    for side,totals in [("left",left_totals),("right",right_totals)]:
        for category,count in totals.items():
            require(normalized_label(f"b-{side}-{category}")==f"{names[category]} {count:,}","Alluvial node label and count agree: "+side+" "+category)
    require(plotted["order_comparison_location"]=="Supplementary Table S3","Input-order results remain linked to the existing supplementary table")
    transition_rows=list(csv.DictReader((PACKAGE/"figures/figure1-annotation-transitions.tsv").open(),delimiter="\t"))
    require(len(transition_rows)==16,"Detailed annotation transitions remain available as sixteen source cells")
    transition_sum=0
    for row in transition_rows:
        old,new=row["GENCODE v29"],row["GENCODE v49"]
        if old==new:
            require(row["Changed observations"]=="","Retained transition table does not encode unchanged observations as zero: "+old)
        else:
            count=annotations["totals"]["classification_transitions"].get(f"{old} -> {new}",0)
            require(int(row["Changed observations"])==count,"Detailed transition source agrees with raw report: "+old+" -> "+new)
            transition_sum+=count
    require(transition_sum==annotations["totals"]["classification_changes"],"Detailed source transitions sum to the total classification changes")
    sys.path.insert(0,str(REPO/"src"))
    from txid.identity import identify
    from txid.models import Exon,TranscriptModel
    schematic={key:TranscriptModel(contig="schematic",strand="+",exons=tuple(Exon(*interval) for interval in data["schematic"][key]),original_transcript_id=key) for key in "ABC"}
    schematic_ids={key:identify(model,"same_schematic_reference") for key,model in schematic.items()}
    require(schematic_ids["A"].splice_chain==schematic_ids["B"].splice_chain and schematic_ids["A"].splice_chain!=schematic_ids["C"].splice_chain,"Schematic SC equality classes agree with identity implementation")
    require(len({item.form.public_id for item in schematic_ids.values()})==3,"Schematic has three distinct exact transcript forms")
    require(schematic["A"].introns==((201,299),) and schematic["C"].introns==((202,299),),"Annotated donor coordinates agree with the underlying closed intron intervals")
    layout=json.loads((PACKAGE/"figures/figure1-layout-check.json").read_text())
    require(layout["width_mm"]==178 and not layout["out_of_bounds_text"] and not layout["overlapping_text"],"Publication geometry is 178 mm wide with no detected text clipping or overlaps")
    require(layout["height_mm"]==122 and layout["minimum_font_pt"]>=5.8 and layout["maximum_body_font_pt"]<=7 and layout["panel_font_pt"]==8,"Composition uses 122 mm height, 5.8–7 pt body text and 8 pt panel labels")
    require(layout["all_text_black"] and not layout["background_gridlines"],"Figure uses black text without background grids")
    require(all(limits[0]==0 for limits in layout["quantitative_y_limits"]),"All quantitative axes start at zero")
    for suffix,dpi in [("png",600),("tif",1200)]:
        with Image.open(PACKAGE/f"figures/figure1.{suffix}") as raster:
            require(all(abs(actual-mm/25.4*dpi)<1.1 for actual,mm in zip(raster.size,[178,122])),"Raster dimensions agree with final physical size: "+suffix)
            require(all(abs(value-dpi)<.1 for value in raster.info["dpi"]),"Raster DPI metadata agrees: "+suffix)
    require(b"/FontFile2" in (PACKAGE/"figures/figure1.pdf").read_bytes(),"Figure PDF contains embedded TrueType font data")
    matrix=list(csv.DictReader((PACKAGE/"tables/table-s4-matrix-dimensions.tsv").open(),delimiter="\t"))
    source_keys=["txid-sorted","gffcompare-sorted","isoseql-exact-sorted","isoseql-junction-sorted"]
    plotted_matrix={row["layer"]:(i,row) for i,row in enumerate(plotted["matrix"])}
    for row,key in zip(matrix,source_keys):
        v=partitions["evaluations"][key]
        require(int(row["Columns"])==v["matrix_columns"],"Matrix columns agree: "+key)
        require(int(row["Occupied cells"])==v["matrix_occupied_cells"],"Matrix occupancy agrees: "+key)
        require(f'{v["matrix_columns"]:,}' in svg_text,"Rendered matrix column count present: "+key)
        i,actual=plotted_matrix[row["Layer"]]
        require(actual["columns"]==v["matrix_columns"] and actual["rows"]==v["matrix_rows"] and svg_label(f"c-columns-{i}")==f'{v["matrix_columns"]:,}',"Matrix plotted value is attached to the correct grouping relation: "+key)
    caption=(PACKAGE/"figure-caption.md").read_text()
    require("Ribbon widths encode transition counts" in caption and "unchanged classifications are omitted" in caption and "sorted-input runs" in caption,"Caption defines the transition subset and matrix run")
    require("Alt text:" in caption and "538,804" in caption and "without database persistence" in caption,"Caption supplies alt text, total denominator and persistence qualification")
    for panel in "abc":
        require("Fig. 1"+panel in "\n".join(path.read_text() for path in body_paths),"Body includes a specific reference to panel "+panel)
    require(round(100*29435/538804,3)==5.463,"Annotation-change percentage reproduces the stated rounding")
    require("reused the TxID parser" in manuscript,"Main text discloses shared parsing")
    require("not byte-identical" in manuscript,"Full-catalog order difference is retained")
    require("without persisting observations" in manuscript,"Partial v49 persistence is disclosed")
    require("form-to-splice-chain associations" in manuscript,"Incremental metadata claim is scoped to the actual check")
    supplement=(PACKAGE/"supplementary.md").read_text()
    supplement_paths=[PACKAGE/name for name in ["supplementary-methods.md","supplementary-tables.md","supplementary-tool-semantics.md"]]
    if (PACKAGE/"supplementary-real-interop.md").is_file():
        supplement_paths.append(PACKAGE/"supplementary-real-interop.md")
    expected_supplement="\n\n".join(path.read_text().strip() for path in supplement_paths)+"\n"
    require(supplement==expected_supplement,"Supplement exactly includes methods, numeric tables, identity semantics and any interoperation evidence")
    with (PACKAGE/"tables/table-s5-identity-semantics.tsv").open(encoding="utf-8",newline="") as handle:
        semantics_rows=list(csv.reader(handle,delimiter="\t"))
    semantics_text=(PACKAGE/"supplementary-tool-semantics.md").read_text()
    markdown_rows=[line for line in semantics_text.splitlines() if line.startswith("|")]
    semantics_md=[[cell.strip() for cell in line.strip("|").split("|")] for line in markdown_rows[:1]+markdown_rows[2:]]
    require(semantics_rows==semantics_md and len(semantics_rows)==6 and all(len(row)==5 for row in semantics_rows),"Table S5 Markdown and TSV agree for all five tools and five columns")
    require({row[0] for row in semantics_rows[1:]}=={"TxID","Isosceles","isoSeQL","TALON","GffCompare"},"Table S5 includes the direct transcript-hash precedent and all comparison workflows")
    require("seed 20260732" in supplement and "seed 20260731" in supplement,"Bootstrap seeds are reported separately")
    require("projected the reference and selected sample annotations" in supplement,"Historical attribute projection is disclosed")
    require("later mutually compatible models can join" in supplement,"Fuzzy description follows observed implementation")
    require("INSERT OR IGNORE" in manuscript,"Comparator compatibility patch is disclosed in main text")
    real_case = verify_real_interop(manuscript, supplement)
    public_access, public_access_limitation = public_access_report()

    docx_names=["manuscript.docx","supplementary.docx", "submission/cover-letter.docx"] + (["manuscript.en.docx"] if bilingual else [])
    for name in docx_names:
        with zipfile.ZipFile(PACKAGE/name) as archive:
            require(archive.testzip() is None,"DOCX ZIP CRC check: "+name)
            for entry in archive.namelist():
                if entry.endswith((".xml",".rels")): ET.fromstring(archive.read(entry))
            require(True,"DOCX XML parts parse: "+name)
            document_tree=ET.fromstring(archive.read("word/document.xml"))
            ns={"w":"http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
            require("word/footer1.xml" in archive.namelist() and
                    b'w:instr="PAGE"' in archive.read("word/footer1.xml"),
                    "DOCX has an automatic page-number footer: " + name)
            relations=ET.fromstring(archive.read("word/_rels/document.xml.rels"))
            relation_ids={node.get("Id") for node in relations}
            require(len(relation_ids)==len(relations), "DOCX relationship identifiers are unique: " + name)
            for node in document_tree.findall(".//w:hyperlink",ns):
                identifier=node.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
                matches=[r for r in relations if r.get("Id")==identifier]
                require(len(matches)==1 and matches[0].get("TargetMode")=="External" and
                        matches[0].get("Target","").startswith(("https://","http://","mailto:")),
                        "Word hyperlink resolves to a supported external target: "+identifier+" in "+name)
            style_tree=ET.fromstring(archive.read("word/styles.xml"))
            default_size=style_tree.find(".//w:rPrDefault/w:rPr/w:sz",ns)
            require(default_size is not None and default_size.get("{"+ns["w"]+"}val")=="24",
                    "DOCX default body font is 12 pt: "+name)
            line_numbers=document_tree.find(".//w:lnNumType",ns)
            require((line_numbers is None)==(name=="submission/cover-letter.docx"),
                    "Line numbering applies to manuscripts and supplement, not the cover letter: "+name)
            if name in {"manuscript.docx","manuscript.en.docx"}:
                require(len([x for x in archive.namelist() if x.startswith("word/media/")])==1,"Main DOCX contains one embedded figure")
                figure_entry=next(x for x in archive.namelist() if x.startswith("word/media/"))
                require(archive.read(figure_entry)==(PACKAGE/"figures/figure1.png").read_bytes(),"DOCX embeds the current figure PNG byte for byte")
                tree=ET.fromstring(archive.read("word/document.xml"))
                text="".join(n.text or "" for n in tree.iter("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t"))
                for line in (display_manuscript if name=="manuscript.docx" else manuscript).splitlines():
                    if not line.strip() or line.startswith(("#","![")): continue
                    line=re.sub(r"\[([^\]]+)\]\(([^)]+)\)",r"\1",line).replace("`","")
                    require(line in text,"DOCX preserves paragraph: "+line[:65])
            elif name=="supplementary.docx":
                tree=ET.fromstring(archive.read("word/document.xml"))
                text="".join(n.text or "" for n in tree.iter("{http://schemas.openxmlformats.org/wordprocessingml/2006/main}t"))
                for path in supplement_paths[2:]:
                    for line in path.read_text().splitlines():
                        if not line.strip(): continue
                        if line.startswith("|"):
                            fragments=[cell.strip() for cell in line.strip("|").split("|")]
                            if all(re.fullmatch(r":?-+:?",fragment) for fragment in fragments): continue
                        else:
                            fragments=[re.sub(r"^#+\s*","",line)]
                        for fragment in fragments:
                            rendered=re.sub(r"\[([^\]]+)\]\(([^)]+)\)",r"\1",fragment).replace("`","")
                            require(rendered in text,"Supplement DOCX preserves reviewed content: "+rendered[:65])
            elif name=="submission/cover-letter.docx":
                text="".join(n.text or "" for n in document_tree.iter("{"+ns["w"]+"}t"))
                for block in (PACKAGE/"submission/cover-letter.md").read_text().strip().split("\n\n"):
                    if block.startswith("#"): continue
                    rendered=re.sub(r"\[([^\]]+)\]\(([^)]+)\)",r"\1",block).replace("*","")
                    require(rendered in text,"Cover-letter Word preserves paragraph: "+rendered[:65])
    ET.parse(PACKAGE/"figures/figure1.svg")
    require((PACKAGE/"figures/figure1.pdf").read_bytes().startswith(b"%PDF"),"Figure PDF and SVG formats are valid")
    counts={path.name:words(path.read_text()) for path in components}
    counts["main_total"]=words(manuscript)
    counts["body_total"]=sum(words(path.read_text()) for path in body_paths)
    summary={"status":"PASS for scientific-draft consistency; not a submission-readiness declaration","checks_passed":len(checks),"word_count_rule":"whitespace-separated visible Markdown text; includes headings, excludes image alt text and URL targets","word_counts":counts,"checks":checks,"current_checkout_tests":software_test_report(),"real_case_recount":real_case,"public_access_check":public_access,"limitations":["Author identities, affiliations, Contact and author-supplied declarations remain incomplete",public_access_limitation,"Persistent software and manuscript-data archive identifiers remain incomplete","Four-page journal pagination not rendered","Full historical benchmarks not rerun","Complete Conda package build and final container image were not tested","PowerShell quality gate unavailable; equivalent Python document/data checks used","DOCX XML/text validated; no office layout renderer available"]}
    summary["author_fields"] = {"authorized_to_remain": True, "remaining": remaining_author_fields}
    summary["document_format"] = "Format-Free: A4, 12 pt body, manuscript line numbering and double spacing, automatic page numbers, active hyperlinks; author completion fields explicit"
    (PACKAGE/"reproducibility/verification-report.json").write_text(json.dumps(summary,ensure_ascii=False,indent=2)+"\n")
    outputs=[PACKAGE/name for name in ["manuscript.md","manuscript.docx","supplementary.md","supplementary.docx","figures/figure1.png","figures/figure1.svg","figures/figure1.pdf","figures/figure1.tif","figures/figure1-layout-check.json","figures/figure1-plot-data.json","figures/figure1-annotation-transitions.tsv","figures/source-data.json","figures/source-data.tsv"]]
    outputs.append(PACKAGE/"submission/cover-letter.docx")
    if bilingual:
        outputs += [PACKAGE/"manuscript.en.md",PACKAGE/"manuscript.en.docx"]
    (PACKAGE/"reproducibility/artifact-sha256.json").write_text(json.dumps({str(path.relative_to(PACKAGE)):hashlib.sha256(path.read_bytes()).hexdigest() for path in outputs},indent=2)+"\n")
    print(json.dumps({"checks_passed":len(checks),"word_counts":counts,"status":"PASS"},indent=2))

if __name__=="__main__": main()
