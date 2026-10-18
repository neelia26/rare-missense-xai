#!/usr/bin/env python3
"""Run final structural, missingness, leakage and split-governance checks."""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    audit_basic_rows,
    file_record,
    is_missing,
    parse_float,
    read_csv,
    write_json,
)


def missingness(rows: list[dict[str, str]], columns: list[str]) -> list[dict[str, Any]]:
    output = []
    by_label = {
        "0": [row for row in rows if row["binary_label"] == "0"],
        "1": [row for row in rows if row["binary_label"] == "1"],
    }
    for column in columns:
        overall = sum(is_missing(row.get(column)) for row in rows) / max(len(rows), 1)
        rates = {
            label: sum(is_missing(row.get(column)) for row in subset) / max(len(subset), 1)
            for label, subset in by_label.items()
        }
        output.append(
            {
                "field": column,
                "overall_missing_fraction": round(overall, 6),
                "benign_missing_fraction": round(rates["0"], 6),
                "pathogenic_missing_fraction": round(rates["1"], 6),
                "absolute_class_gap": round(abs(rates["0"] - rates["1"]), 6),
                "shortcut_warning": abs(rates["0"] - rates["1"]) >= 0.10,
            }
        )
    return output


def split_sets(split_dir: Path) -> tuple[dict[str, set[str]], dict[str, int]]:
    genes: dict[str, set[str]] = {}
    sizes: dict[str, int] = {}
    for split in ("train", "validation", "test"):
        rows = read_csv(split_dir / f"{split}_full.csv")
        genes[split] = {row["reconciled_gene_symbol"] for row in rows}
        sizes[split] = len(rows)
    return genes, sizes


def range_violations(rows: list[dict[str, str]]) -> dict[str, int]:
    checks = {
        "gnomad_global_af": lambda value: 0 <= value <= 1,
        "gnomad_max_population_af": lambda value: 0 <= value <= 1,
        "aa_relative_protein_position": lambda value: 0 < value <= 1,
        "interpro_overlap_count": lambda value: value >= 0,
    }
    violations: dict[str, int] = {}
    for column, predicate in checks.items():
        if column not in rows[0]:
            continue
        count = 0
        for row in rows:
            value = parse_float(row.get(column))
            if value is not None and not predicate(value):
                count += 1
        violations[column] = count
    return violations


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", default=str(PROJECT_ROOT / "data/processed/final_dataset_full.csv"))
    parser.add_argument("--core", default=str(PROJECT_ROOT / "data/processed/model_table_core.csv"))
    parser.add_argument("--registry", default=str(PROJECT_ROOT / "data/processed/evidence_feature_registry.csv"))
    parser.add_argument("--split-dir", default=str(PROJECT_ROOT / "data/processed/splits"))
    args = parser.parse_args()
    full_path = Path(args.full)
    core_path = Path(args.core)
    registry_path = Path(args.registry)
    rows = read_csv(full_path)
    core = read_csv(core_path)
    registry = read_csv(registry_path)
    model_features = [
        row["field"] for row in registry if row["default_core_model_feature"] == "1"
    ]
    critical_failures: list[str] = []
    warnings: list[str] = []

    basic = audit_basic_rows(rows)
    if basic["duplicate_key_rows"]:
        critical_failures.append("Duplicate variant keys remain")
    if set(basic["class_counts"]) != {"0", "1"}:
        critical_failures.append("Labels are not strictly binary")
    if basic["class_counts"].get("0") != basic["class_counts"].get("1"):
        critical_failures.append("Final classes are not balanced")
    if {row["variant_key"] for row in rows} != {row["variant_key"] for row in core}:
        critical_failures.append("Full and core model tables contain different variant keys")
    if any(column.startswith("predictor_") for column in core[0]):
        critical_failures.append("Composite computational predictor columns leaked into the core table")

    gene_sets, split_sizes = split_sets(Path(args.split_dir))
    overlaps = {
        "train_validation": sorted(gene_sets["train"] & gene_sets["validation"]),
        "train_test": sorted(gene_sets["train"] & gene_sets["test"]),
        "validation_test": sorted(gene_sets["validation"] & gene_sets["test"]),
    }
    if any(overlaps.values()):
        critical_failures.append("At least one gene occurs in more than one data partition")

    missing_report = missingness(core, model_features)
    shortcut_fields = [row["field"] for row in missing_report if row["shortcut_warning"]]
    if shortcut_fields:
        warnings.append(
            "Class-associated missingness (absolute gap >= 0.10) requires controlled missingness experiments: "
            + ", ".join(shortcut_fields)
        )
    all_missing = [column for column in model_features if all(is_missing(row.get(column)) for row in core)]
    constant = [
        column
        for column in model_features
        if len({row.get(column, "") for row in core if not is_missing(row.get(column))}) <= 1
    ]
    if all_missing:
        warnings.append("All-missing default features: " + ", ".join(all_missing))
    if constant:
        warnings.append("Constant or single-valued default features: " + ", ".join(constant))
    violations = range_violations(rows)
    if any(violations.values()):
        critical_failures.append("At least one numeric field violates its valid range")

    status = "PASS" if not critical_failures else "FAIL"
    report = {
        "status": status,
        "critical_failures": critical_failures,
        "warnings": warnings,
        "dataset": basic,
        "files": {
            "full": file_record(full_path),
            "core": file_record(core_path),
            "registry": file_record(registry_path),
        },
        "default_model_features": model_features,
        "missingness_by_class": missing_report,
        "class_associated_missingness_fields": shortcut_fields,
        "all_missing_features": all_missing,
        "constant_features": constant,
        "numeric_range_violations": violations,
        "split_sizes": split_sizes,
        "split_gene_counts": {split: len(values) for split, values in gene_sets.items()},
        "split_gene_overlap": overlaps,
        "governance_decision": (
            "A PASS confirms structural readiness, not scientific validity. Any warning must be "
            "resolved or tested through controlled feature/missingness experiments before modelling."
        ),
    }
    json_path = PROJECT_ROOT / "audits/14_final_dataset_audit.json"
    write_json(json_path, report)
    markdown_path = PROJECT_ROOT / "audits/14_final_dataset_audit.md"
    lines = [
        "# Final Dataset Audit",
        "",
        f"**Status: {status}**",
        "",
        f"- Rows: {basic['rows']:,}",
        f"- Unique variant keys: {basic['unique_keys']:,}",
        f"- Class counts: {basic['class_counts']}",
        f"- Split sizes: {split_sizes}",
        f"- Gene overlap: {'none' if not any(overlaps.values()) else 'present'}",
        "",
        "## Critical failures",
        "",
        *(f"- {item}" for item in critical_failures),
        *( ["- None"] if not critical_failures else [] ),
        "",
        "## Warnings requiring a modelling decision",
        "",
        *(f"- {item}" for item in warnings),
        *( ["- None"] if not warnings else [] ),
        "",
        "Structural PASS does not remove the need for feature-ablation, missingness and circularity experiments.",
    ]
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    if critical_failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
