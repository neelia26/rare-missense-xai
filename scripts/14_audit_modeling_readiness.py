#!/usr/bin/env python3
"""Audit the derived, leakage-controlled modelling cohort and gene splits."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    audit_basic_rows,
    file_record,
    is_missing,
    read_csv,
    sha256_file,
    write_json,
)


EXPECTED_SOURCE_SHA256 = "50d8ac02ededae696d6867c7e713351d33058d01bad4fab3633055157525494f"
PROHIBITED_UNIPROT_COLUMNS = {
    "uniprot_feature_natural_variant",
    "uniprot_feature_mutagenesis",
    "uniprot_feature_sequence_conflict",
    "uniprot_feature_alternative_sequence",
}
RAW_CONTEXT_COLUMNS = {
    "interpro_overlap_types",
    "uniprot_overlapping_feature_types",
}


def split_report(split_dir: Path) -> tuple[dict[str, Any], dict[str, set[str]]]:
    report: dict[str, Any] = {}
    genes: dict[str, set[str]] = {}
    for split in ("train", "validation", "test"):
        rows = read_csv(split_dir / f"{split}_full.csv")
        labels = Counter(row["binary_label"] for row in rows)
        genes[split] = {row["reconciled_gene_symbol"] for row in rows}
        report[split] = {
            "variants": len(rows),
            "variant_fraction": 0.0,
            "genes": len(genes[split]),
            "class_counts": dict(labels),
            "pathogenic_proportion": labels["1"] / max(len(rows), 1),
        }
    total = sum(item["variants"] for item in report.values())
    for item in report.values():
        item["variant_fraction"] = item["variants"] / max(total, 1)
    return report, genes


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source", default=str(PROJECT_ROOT / "data/processed/final_dataset_full.csv")
    )
    parser.add_argument(
        "--full", default=str(PROJECT_ROOT / "data/processed/modeling_ready_full.csv")
    )
    parser.add_argument(
        "--core", default=str(PROJECT_ROOT / "data/processed/modeling_ready_core.csv")
    )
    parser.add_argument(
        "--no-population",
        default=str(PROJECT_ROOT / "data/processed/modeling_ready_core_no_population.csv"),
    )
    parser.add_argument(
        "--registry", default=str(PROJECT_ROOT / "data/processed/modeling_feature_registry.csv")
    )
    parser.add_argument(
        "--excluded",
        default=str(PROJECT_ROOT / "data/processed/excluded_alternate_gene_variants.csv"),
    )
    parser.add_argument(
        "--split-dir", default=str(PROJECT_ROOT / "data/processed/modeling_splits")
    )
    parser.add_argument("--minimum-evaluation-genes", type=int, default=1100)
    parser.add_argument("--variant-fraction-tolerance", type=float, default=0.01)
    parser.add_argument("--class-proportion-tolerance", type=float, default=0.02)
    args = parser.parse_args()

    source_path = Path(args.source)
    full_path = Path(args.full)
    core_path = Path(args.core)
    no_population_path = Path(args.no_population)
    registry_path = Path(args.registry)
    excluded_path = Path(args.excluded)
    rows = read_csv(full_path)
    core = read_csv(core_path)
    no_population = read_csv(no_population_path)
    registry = read_csv(registry_path)
    excluded = read_csv(excluded_path)
    failures: list[str] = []
    warnings: list[str] = []

    if sha256_file(source_path) != EXPECTED_SOURCE_SHA256:
        failures.append("Frozen 40,000-row source checksum changed")
    if len(rows) != 39552 or len(excluded) != 448:
        failures.append("Alternate-gene exclusion did not produce 39,552 retained and 448 excluded rows")
    if any(row["gene_symbol"].strip() != row["reconciled_gene_symbol"].strip() for row in rows):
        failures.append("Alternate-gene rows remain in the primary cohort")
    if {row["variant_key"] for row in rows} != {row["variant_key"] for row in core}:
        failures.append("Derived full and core tables contain different variant keys")
    if {row["variant_key"] for row in rows} != {row["variant_key"] for row in no_population}:
        failures.append("No-population ablation contains different variant keys")

    core_columns = set(core[0])
    if core_columns & RAW_CONTEXT_COLUMNS:
        failures.append("Raw semicolon-separated protein-context fields remain in the core table")
    if core_columns & PROHIBITED_UNIPROT_COLUMNS:
        failures.append("Leakage-prone UniProt tokens remain in the primary core table")
    if any("missing" in column.lower() or "available" in column.lower() for column in core_columns):
        failures.append("A missingness/availability indicator is present in the primary core table")
    population_columns = {"gnomad_global_af", "gnomad_max_population_af"}
    if population_columns & set(no_population[0]):
        failures.append("Population features remain in the no-population ablation")
    if not population_columns <= core_columns:
        failures.append("Primary core table unexpectedly lacks population evidence")
    if not any(is_missing(row["gnomad_global_af"]) for row in core):
        failures.append("Unavailable population evidence was not preserved as missing")

    interpro_zero = sum(
        row["interpro_lookup_status"] == "FOUND"
        and all(row[column] == "0" for column in core_columns if column.startswith("interpro_type_"))
        for row in rows
    )
    interpro_failures = sum(row["interpro_lookup_status"] != "FOUND" for row in rows)
    failure_rows = [row for row in rows if row["interpro_lookup_status"] != "FOUND"]
    if any(
        not is_missing(row[column])
        for row in failure_rows
        for column in core_columns
        if column.startswith("interpro_type_")
    ):
        failures.append("InterPro lookup failures are not encoded as missing token indicators")

    splits, gene_sets = split_report(Path(args.split_dir))
    overlaps = {
        "train_validation": sorted(gene_sets["train"] & gene_sets["validation"]),
        "train_test": sorted(gene_sets["train"] & gene_sets["test"]),
        "validation_test": sorted(gene_sets["validation"] & gene_sets["test"]),
    }
    if any(overlaps.values()):
        failures.append("Gene overlap remains between partitions")
    for split, target in (("train", 0.70), ("validation", 0.15), ("test", 0.15)):
        if abs(splits[split]["variant_fraction"] - target) > args.variant_fraction_tolerance:
            failures.append(f"{split} variant fraction is outside tolerance")
        if abs(splits[split]["pathogenic_proportion"] - 0.50) > args.class_proportion_tolerance:
            failures.append(f"{split} pathogenic proportion is outside tolerance")
    for split in ("validation", "test"):
        if splits[split]["genes"] < args.minimum_evaluation_genes:
            failures.append(f"{split} has fewer than {args.minimum_evaluation_genes} genes")

    registry_primary = {
        row["field"] for row in registry if row.get("primary_model_feature") == "1"
    }
    if registry_primary != core_columns - {"variant_key", "reconciled_gene_symbol", "binary_label"}:
        failures.append("Feature registry does not exactly match the primary core features")

    population_missing = {
        label: sum(
            is_missing(row["gnomad_global_af"])
            for row in core
            if row["binary_label"] == label
        )
        / max(sum(row["binary_label"] == label for row in core), 1)
        for label in ("0", "1")
    }
    if abs(population_missing["0"] - population_missing["1"]) >= 0.10:
        warnings.append(
            "Population missingness remains class-associated; enforce train-only neutral imputation "
            "without an availability indicator and report the no-population ablation."
        )

    status = "PASS" if not failures else "FAIL"
    report = {
        "status": status,
        "critical_failures": failures,
        "warnings": warnings,
        "frozen_source": file_record(source_path),
        "derived_dataset": audit_basic_rows(rows),
        "alternate_gene_rows_excluded": len(excluded),
        "protein_context": {
            "raw_context_fields_excluded_from_core": sorted(RAW_CONTEXT_COLUMNS),
            "prohibited_uniprot_tokens_excluded": sorted(PROHIBITED_UNIPROT_COLUMNS),
            "interpro_successful_zero_overlap": interpro_zero,
            "interpro_lookup_failures": interpro_failures,
        },
        "population_missing_fraction_by_class": population_missing,
        "population_handling": (
            "No availability indicator is supplied. Fit any imputation/scaling on training genes only; "
            "unavailable evidence must contribute neutrally and be reported as unavailable. A separate "
            "no-population core table is supplied for mandatory ablation."
        ),
        "splits": splits,
        "gene_overlap": overlaps,
        "acceptance_thresholds": {
            "minimum_evaluation_genes": args.minimum_evaluation_genes,
            "variant_fraction_tolerance": args.variant_fraction_tolerance,
            "class_proportion_tolerance": args.class_proportion_tolerance,
        },
        "files": {
            "full": file_record(full_path),
            "core": file_record(core_path),
            "no_population": file_record(no_population_path),
            "registry": file_record(registry_path),
            "excluded": file_record(excluded_path),
        },
        "modelling_readiness": (
            "READY for baseline modelling only if status is PASS and the population-handling warning "
            "is enforced through the documented preprocessing and ablation protocol."
        ),
    }
    output = PROJECT_ROOT / "audits/14_modeling_readiness_audit.json"
    write_json(output, report)
    print(json.dumps(report, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
