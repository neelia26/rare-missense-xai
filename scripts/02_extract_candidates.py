#!/usr/bin/env python3
"""Extract unambiguous labelled biallelic missense SNVs from the frozen ClinVar VCF."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    load_config,
    variant_key,
    write_csv,
    write_json,
)


OUTPUT_COLUMNS = [
    "variant_key",
    "genome_build",
    "chromosome",
    "position",
    "reference_allele",
    "alternate_allele",
    "clinvar_variation_id",
    "clinvar_allele_id",
    "dbsnp_id",
    "gene_symbol",
    "gene_id",
    "all_gene_info",
    "genomic_hgvs",
    "molecular_consequence",
    "clinvar_significance",
    "binary_label",
    "review_status",
    "review_stars",
    "clinical_conditions",
    "candidate_source_release",
]


REVIEW_STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_single_submitter": 1,
    "no_assertion_criteria_provided": 0,
    "no_assertion_provided": 0,
    "no_classification_provided": 0,
}


def parse_info(raw: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in raw.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = value
        else:
            parsed[item] = ""
    return parsed


def first_gene(raw: str) -> tuple[str, str]:
    first = raw.split("|", 1)[0]
    if ":" not in first:
        return first, ""
    return tuple(first.split(":", 1))  # type: ignore[return-value]


def row_from_vcf(
    fields: list[str], labels: dict[str, int], minimum_review_stars: int, release: str
) -> tuple[dict[str, object] | None, str]:
    if len(fields) < 8:
        return None, "malformed_vcf_row"
    chromosome, position, variation_id, reference, alternate, _, _, raw_info = fields[:8]
    reference = reference.upper()
    alternate = alternate.upper()
    if len(reference) != 1 or len(alternate) != 1 or "," in alternate:
        return None, "not_biallelic_snv"
    if reference not in "ACGT" or alternate not in "ACGT" or reference == alternate:
        return None, "invalid_allele"
    info = parse_info(raw_info)
    if info.get("CLNVC") != "single_nucleotide_variant":
        return None, "not_clinvar_snv"
    molecular = info.get("MC", "")
    if "missense_variant" not in molecular:
        return None, "not_clinvar_missense"
    significance = info.get("CLNSIG", "")
    if significance not in labels:
        return None, "ineligible_or_nonbinary_label"
    if info.get("CLNSIGCONF"):
        return None, "conflicting_label"
    review = info.get("CLNREVSTAT", "")
    stars = REVIEW_STARS.get(review, 0)
    if stars < minimum_review_stars:
        return None, "review_status_below_minimum"
    gene_symbol, gene_id = first_gene(info.get("GENEINFO", ""))
    if not gene_symbol:
        return None, "missing_gene"
    return {
        "variant_key": variant_key(chromosome, position, reference, alternate),
        "genome_build": "GRCh38",
        "chromosome": chromosome,
        "position": int(position),
        "reference_allele": reference,
        "alternate_allele": alternate,
        "clinvar_variation_id": variation_id,
        "clinvar_allele_id": info.get("ALLELEID", ""),
        "dbsnp_id": f"rs{info['RS']}" if info.get("RS") else "",
        "gene_symbol": gene_symbol,
        "gene_id": gene_id,
        "all_gene_info": info.get("GENEINFO", ""),
        "genomic_hgvs": info.get("CLNHGVS", ""),
        "molecular_consequence": molecular,
        "clinvar_significance": significance.replace("_", " "),
        "binary_label": labels[significance],
        "review_status": review.replace("_", " "),
        "review_stars": stars,
        "clinical_conditions": info.get("CLNDN", "").replace("_", " "),
        "candidate_source_release": release,
    }, "eligible"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/raw/clinvar_grch38.vcf.gz"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/02_clinvar_candidates.csv"))
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    parser.add_argument("--minimum-review-stars", type=int)
    args = parser.parse_args()
    config = load_config(args.config)
    labels = config["clinvar"]["accepted_significance"]
    minimum_stars = (
        args.minimum_review_stars
        if args.minimum_review_stars is not None
        else int(config["clinvar"]["minimum_review_stars"])
    )
    input_path = Path(args.input)
    output_path = Path(args.output)
    exclusions: Counter[str] = Counter()
    release = ""
    deduplicated: dict[str, dict[str, object]] = {}
    duplicate_label_conflicts: set[str] = set()
    scanned = 0

    with gzip.open(input_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##fileDate="):
                release = line.strip().split("=", 1)[1]
            if line.startswith("#"):
                continue
            scanned += 1
            row, status = row_from_vcf(
                line.rstrip("\n").split("\t"), labels, minimum_stars, release
            )
            if row is None:
                exclusions[status] += 1
                continue
            key = str(row["variant_key"])
            prior = deduplicated.get(key)
            if prior is None:
                deduplicated[key] = row
                continue
            if prior["binary_label"] != row["binary_label"]:
                duplicate_label_conflicts.add(key)
                continue
            prior_rank = (int(prior["review_stars"]), -int(prior["clinvar_variation_id"]))
            row_rank = (int(row["review_stars"]), -int(row["clinvar_variation_id"]))
            if row_rank > prior_rank:
                deduplicated[key] = row

    for key in duplicate_label_conflicts:
        deduplicated.pop(key, None)
    exclusions["duplicate_label_conflict"] += len(duplicate_label_conflicts)
    rows = sorted(
        deduplicated.values(),
        key=lambda row: (
            int(row["binary_label"]),
            str(row["chromosome"]),
            int(row["position"]),
            str(row["reference_allele"]),
            str(row["alternate_allele"]),
        ),
    )
    write_csv(output_path, rows, OUTPUT_COLUMNS)
    audit = {
        **audit_basic_rows(rows),
        "source_vcf_records_scanned": scanned,
        "clinvar_file_date": release,
        "minimum_review_stars": minimum_stars,
        "exclusion_counts": dict(exclusions),
        "selection_rules": [
            "GRCh38 biallelic single-nucleotide substitutions",
            "ClinVar molecular consequence contains missense_variant",
            "unambiguous benign/likely-benign or pathogenic/likely-pathogenic aggregate label",
            f"ClinVar review status has at least {minimum_stars} star(s)",
            "one deterministic record per genomic variant key",
        ],
    }
    audit_path = PROJECT_ROOT / "audits/02_clinvar_candidates_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "02_extract_candidates",
        [input_path],
        [output_path, audit_path],
        {"minimum_review_stars": minimum_stars, "accepted_significance": labels},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

