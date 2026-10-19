#!/usr/bin/env python3
"""Static circularity, leakage, missingness and evidence-overlap audit.

Label-dependent analyses are deliberately restricted to the locked training split.
The validation and test splits are read only for structural integrity checks.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import PROJECT_ROOT, read_csv, sha256_file, write_json


FROZEN_SHA256 = "50d8ac02ededae696d6867c7e713351d33058d01bad4fab3633055157525494f"
METADATA = {"variant_key", "reconciled_gene_symbol", "binary_label"}
MISSING_MARKERS = {"", ".", "na", "nan", "none", "null"}
PROHIBITED_UNIPROT = {
    "uniprot_feature_natural_variant",
    "uniprot_feature_mutagenesis",
    "uniprot_feature_sequence_conflict",
    "uniprot_feature_alternative_sequence",
}
RAW_CONTEXT = {"interpro_overlap_types", "uniprot_overlapping_feature_types"}
POPULATION = {"gnomad_global_af", "gnomad_max_population_af"}
EXPECTED_SPLITS = {
    "train": {"variants": 27727, "genes": 5206, "0": 13741, "1": 13986},
    "validation": {"variants": 5790, "genes": 1115, "0": 2895, "1": 2895},
    "test": {"variants": 6035, "genes": 1115, "0": 3006, "1": 3029},
}
MANDATORY_ABLATIONS = [
    "no-population versus population-inclusive",
    "without gene-constraint evidence",
    "evidence-group removal",
    "missingness-only diagnostic model",
    "optional computational predictors only after their circularity review",
]


def is_missing(value: Any) -> bool:
    return value is None or str(value).strip().lower() in MISSING_MARKERS


def parse_number(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        parsed = float(str(value).strip())
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def write_csv_atomic(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def structural_rows(rows_iterable: Iterable[dict[str, str]], schema: tuple[str, ...]) -> dict[str, Any]:
    """Summarise structure without reading labels from any row."""
    keys: set[str] = set()
    genes: set[str] = set()
    blank_keys = 0
    rows = 0
    for row in rows_iterable:
        rows += 1
        key = row.get("variant_key", "")
        gene = row.get("reconciled_gene_symbol", "")
        blank_keys += int(not key)
        keys.add(key)
        genes.add(gene)
    return {
        "rows": rows,
        "keys": keys,
        "genes": genes,
        "blank_keys": blank_keys,
        "schema": schema,
    }


def structural_split(path: Path) -> dict[str, Any]:
    """Read no labels: only structural fields permitted for validation/test checks."""
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        schema = tuple(reader.fieldnames or [])
        return structural_rows(reader, schema)


def missingness_rows(
    training: list[dict[str, str]], features: list[str]
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    by_label = {
        label: [row for row in training if row["binary_label"] == label]
        for label in ("0", "1")
    }
    for feature in features:
        missing_count = sum(is_missing(row.get(feature)) for row in training)
        rates = {
            label: sum(is_missing(row.get(feature)) for row in rows) / max(len(rows), 1)
            for label, rows in by_label.items()
        }
        gap = abs(rates["0"] - rates["1"])
        output.append(
            {
                "field": feature,
                "training_rows": len(training),
                "missing_count": missing_count,
                "missing_fraction": round(missing_count / max(len(training), 1), 8),
                "benign_missing_fraction": round(rates["0"], 8),
                "pathogenic_missing_fraction": round(rates["1"], 8),
                "absolute_class_gap": round(gap, 8),
                "class_associated_missingness_flag": int(gap >= 0.10),
                "analysis_partition": "train",
            }
        )
    return output


def rank_values(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    index = 0
    while index < len(values):
        end = index + 1
        while end < len(values) and values[order[end]] == values[order[index]]:
            end += 1
        average_rank = (index + end - 1) / 2.0 + 1.0
        ranks[order[index:end]] = average_rank
        index = end
    return ranks


def safe_correlation(left: np.ndarray, right: np.ndarray) -> tuple[float | None, str]:
    if len(left) < 3:
        return None, "INSUFFICIENT_PAIRED_OBSERVATIONS"
    if np.ptp(left) == 0 or np.ptp(right) == 0:
        return None, "CONSTANT_INPUT"
    value = float(np.corrcoef(left, right)[0, 1])
    if not math.isfinite(value):
        return None, "UNDEFINED"
    return value, "VALID"


def numeric_correlations(
    training: list[dict[str, str]], features: list[str]
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for index, left_name in enumerate(features):
        for right_name in features[index + 1 :]:
            pairs = [
                (left, right)
                for row in training
                if (left := parse_number(row.get(left_name))) is not None
                and (right := parse_number(row.get(right_name))) is not None
            ]
            if pairs:
                left = np.asarray([pair[0] for pair in pairs], dtype=float)
                right = np.asarray([pair[1] for pair in pairs], dtype=float)
            else:
                left = right = np.asarray([], dtype=float)
            pearson, pearson_status = safe_correlation(left, right)
            spearman, spearman_status = safe_correlation(
                rank_values(left) if len(left) else left,
                rank_values(right) if len(right) else right,
            )
            valid = pearson is not None or spearman is not None
            redundancy = valid and max(abs(pearson or 0.0), abs(spearman or 0.0)) >= 0.85
            output.append(
                {
                    "field_a": left_name,
                    "field_b": right_name,
                    "paired_observations": len(pairs),
                    "pearson": "" if pearson is None else round(pearson, 8),
                    "spearman": "" if spearman is None else round(spearman, 8),
                    "pearson_status": pearson_status,
                    "spearman_status": spearman_status,
                    "redundancy_candidate": int(redundancy),
                    "analysis_partition": "train",
                }
            )
    return output


def corrected_cramers_v(
    left: Iterable[Any], right: Iterable[Any]
) -> tuple[float | None, int, str]:
    pairs = [
        (str(a).strip(), str(b).strip())
        for a, b in zip(left, right)
        if not is_missing(a) and not is_missing(b)
    ]
    n = len(pairs)
    if n < 3:
        return None, n, "INSUFFICIENT_PAIRED_OBSERVATIONS"
    left_levels = sorted({a for a, _ in pairs})
    right_levels = sorted({b for _, b in pairs})
    if len(left_levels) < 2 or len(right_levels) < 2:
        return None, n, "CONSTANT_INPUT"
    counts = Counter(pairs)
    left_totals = Counter(a for a, _ in pairs)
    right_totals = Counter(b for _, b in pairs)
    chi2 = 0.0
    for a in left_levels:
        for b in right_levels:
            expected = left_totals[a] * right_totals[b] / n
            if expected:
                chi2 += (counts[(a, b)] - expected) ** 2 / expected
    phi2 = chi2 / n
    rows, columns = len(left_levels), len(right_levels)
    phi2_corrected = max(0.0, phi2 - ((columns - 1) * (rows - 1)) / (n - 1))
    rows_corrected = rows - ((rows - 1) ** 2) / (n - 1)
    columns_corrected = columns - ((columns - 1) ** 2) / (n - 1)
    denominator = min(columns_corrected - 1, rows_corrected - 1)
    if denominator <= 0:
        return None, n, "UNDEFINED_AFTER_BIAS_CORRECTION"
    value = math.sqrt(phi2_corrected / denominator)
    return value, n, "VALID"


def categorical_associations(
    training: list[dict[str, str]], features: list[str]
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for index, left in enumerate(features):
        for right in features[index + 1 :]:
            value, paired, status = corrected_cramers_v(
                (row.get(left) for row in training),
                (row.get(right) for row in training),
            )
            output.append(
                {
                    "field_a": left,
                    "field_b": right,
                    "paired_observations": paired,
                    "corrected_cramers_v": "" if value is None else round(value, 8),
                    "status": status,
                    "extremely_sparse": int(paired < max(100, int(0.01 * len(training)))),
                    "redundancy_candidate": int(value is not None and value >= 0.80),
                    "analysis_partition": "train",
                }
            )
    return output


def assignment_summary(rows: list[dict[str, str]]) -> dict[str, Counter[str]]:
    summary: dict[str, Counter[str]] = {}
    for row in rows:
        split = row["split"]
        values = summary.setdefault(split, Counter())
        values["genes"] += 1
        values["variants"] += int(row["variant_count"])
        values["0"] += int(row["benign_count"])
        values["1"] += int(row["pathogenic_count"])
    return summary


def audit(args: argparse.Namespace) -> dict[str, Any]:
    frozen_path = Path(args.frozen)
    full_path = Path(args.full)
    core_path = Path(args.core)
    no_population_path = Path(args.no_population)
    registry_path = Path(args.registry)
    lineage_path = Path(args.lineage)
    assignment_path = Path(args.assignments)
    split_dir = Path(args.split_dir)
    output_dir = Path(args.output_dir)

    failures: list[str] = []
    warnings: list[str] = []
    unresolved: list[str] = []
    if sha256_file(frozen_path) != FROZEN_SHA256:
        failures.append("Frozen-source checksum changed")

    frozen = read_csv(frozen_path)
    full = read_csv(full_path)
    core = read_csv(core_path)
    no_population = read_csv(no_population_path)
    registry = read_csv(registry_path)
    lineage = read_csv(lineage_path)
    assignments = read_csv(assignment_path)
    if len(frozen) != 40000 or Counter(row["binary_label"] for row in frozen) != Counter({"0": 20000, "1": 20000}):
        failures.append("Frozen-source row or class counts differ from the governed values")
    if len(full) != 39552:
        failures.append("Modelling cohort does not contain 39,552 variants")
    for name, rows in (("full", full), ("core", core), ("no_population", no_population)):
        keys = [row.get("variant_key", "") for row in rows]
        if not all(keys) or len(keys) != len(set(keys)):
            failures.append(f"{name} contains blank or duplicate variant keys")
    full_keys = {row["variant_key"] for row in full}
    if full_keys != {row["variant_key"] for row in core} or full_keys != {row["variant_key"] for row in no_population}:
        failures.append("Derived modelling tables contain different variant keys")

    primary_rows = [row for row in registry if row["primary_model_feature"] == "1"]
    primary = sorted(row["field"] for row in primary_rows)
    if set(primary) != set(core[0]) - METADATA:
        failures.append("Primary feature registry does not exactly correspond to core predictors")
    if METADATA & set(primary):
        failures.append("Metadata or target is registered as a predictor")
    if PROHIBITED_UNIPROT & set(primary) or PROHIBITED_UNIPROT & set(core[0]):
        failures.append("Prohibited UniProt curated-variant feature is present in primary inputs")
    if RAW_CONTEXT & set(core[0]):
        failures.append("Raw semicolon-separated context field is present in primary inputs")
    no_population_predictors = set(no_population[0]) - METADATA
    if POPULATION & no_population_predictors:
        failures.append("No-population baseline contains population features")
    if any("missing" in field.lower() or "available" in field.lower() for field in primary):
        failures.append("Primary inputs contain a missingness or availability indicator")

    lineage_primary = sorted(
        row["field"] for row in lineage if row["approved_primary"] == "1"
    )
    if lineage_primary != primary:
        failures.append("Evidence-lineage primary features differ from the modelling registry")
    optional_fields = sorted(
        field for field in structural_split(split_dir / "train_optional_predictors.csv")["schema"]
        if field.startswith("predictor_")
    )
    lineage_fields = {row["field"]: row for row in lineage}
    for field in optional_fields:
        if field not in lineage_fields:
            failures.append(f"Optional predictor missing from lineage registry: {field}")
        elif lineage_fields[field]["approved_primary"] == "1":
            failures.append(f"Optional predictor incorrectly approved as primary: {field}")

    split_structures = {
        split: structural_split(split_dir / f"{split}_full.csv")
        for split in ("train", "validation", "test")
    }
    assignment_counts = assignment_summary(assignments)
    for split, expected in EXPECTED_SPLITS.items():
        structure = split_structures[split]
        observed = assignment_counts.get(split, Counter())
        if structure["rows"] != expected["variants"] or len(structure["genes"]) != expected["genes"]:
            failures.append(f"{split} structural counts differ from locked values")
        if any(observed[key] != expected[key] for key in ("variants", "genes", "0", "1")):
            failures.append(f"{split} assignment summary differs from locked values")
    overlaps = {
        "train_validation_genes": sorted(split_structures["train"]["genes"] & split_structures["validation"]["genes"]),
        "train_test_genes": sorted(split_structures["train"]["genes"] & split_structures["test"]["genes"]),
        "validation_test_genes": sorted(split_structures["validation"]["genes"] & split_structures["test"]["genes"]),
        "train_validation_variants": sorted(split_structures["train"]["keys"] & split_structures["validation"]["keys"]),
        "train_test_variants": sorted(split_structures["train"]["keys"] & split_structures["test"]["keys"]),
        "validation_test_variants": sorted(split_structures["validation"]["keys"] & split_structures["test"]["keys"]),
    }
    if any(overlaps.values()):
        failures.append("Gene or variant overlap exists between locked partitions")

    # This is the only table used for label-dependent analysis.
    training = read_csv(split_dir / "train_core.csv")
    numeric = sorted(row["field"] for row in primary_rows if row["data_type"] == "numeric")
    categorical = sorted(row["field"] for row in primary_rows if row["data_type"] in {"binary", "categorical"})
    missing = missingness_rows(training, primary)
    numeric_rows = numeric_correlations(training, numeric)
    categorical_rows = categorical_associations(training, categorical)

    population_missing = {
        row["field"]: row for row in missing if row["field"] in POPULATION
    }
    if not POPULATION <= set(population_missing) or not all(
        row["class_associated_missingness_flag"] == 1 for row in population_missing.values()
    ):
        failures.append("Expected class-associated population missingness was not confirmed on training data")
    missing_flags = [row for row in missing if row["class_associated_missingness_flag"]]
    if missing_flags:
        warnings.append("Training-only class-associated missingness was detected; flags require controls, not automatic feature removal")

    redundancy: list[dict[str, Any]] = []
    for row in numeric_rows:
        if row["redundancy_candidate"]:
            redundancy.append({"association_type": "numeric", **row})
    for row in categorical_rows:
        if row["redundancy_candidate"]:
            redundancy.append({"association_type": "categorical", **row})
    redundancy.sort(key=lambda row: (row["association_type"], row["field_a"], row["field_b"]))
    if redundancy:
        warnings.append("Strong within-training feature overlap was detected; candidates require grouped ablations, not automatic deletion")

    unresolved = sorted(
        row["field"]
        for row in lineage
        if "UNRESOLVED_REQUIRES_LITERATURE_REVIEW" in (
            row.get("predictor_training_overlap_risk", ""),
            row.get("label_provenance_overlap_risk", ""),
        )
    )
    if unresolved:
        warnings.append("Optional predictor overlap remains unresolved pending literature review")

    governance = []
    for row in lineage:
        if row["prohibited"] == "1":
            decision = "PROHIBITED"
        elif row["approved_sensitivity"] == "1":
            decision = "SENSITIVITY_ONLY"
        elif row["approved_primary"] == "1":
            decision = "APPROVED_PRIMARY"
        else:
            decision = "UNRESOLVED"
        governance.append(
            {
                "field": row["field"],
                "decision": decision,
                "evidence_group": row["evidence_group"],
                "required_control": row["required_control"],
                "decision_reason": row["decision_reason"],
            }
        )
    governance.sort(key=lambda row: row["field"])

    output_dir.mkdir(parents=True, exist_ok=True)
    write_csv_atomic(output_dir / "15_missingness_by_class.csv", missing, list(missing[0]))
    write_csv_atomic(output_dir / "15_numeric_correlations.csv", numeric_rows, list(numeric_rows[0]))
    write_csv_atomic(output_dir / "15_categorical_associations.csv", categorical_rows, list(categorical_rows[0]))
    redundancy_fields = [
        "association_type", "field_a", "field_b", "paired_observations", "pearson",
        "spearman", "corrected_cramers_v", "status", "redundancy_candidate", "analysis_partition",
    ]
    normalized_redundancy = []
    for row in redundancy:
        normalized_redundancy.append({field: row.get(field, "") for field in redundancy_fields})
    write_csv_atomic(output_dir / "15_redundancy_candidates.csv", normalized_redundancy, redundancy_fields)
    write_csv_atomic(output_dir / "15_feature_governance_decisions.csv", governance, list(governance[0]))

    status = "PASS" if not failures else "FAIL"
    report = {
        "structural_status": status,
        "scientific_clearance": "CONDITIONAL",
        "critical_failures": failures,
        "warnings": warnings,
        "unresolved_items": unresolved,
        "required_modelling_controls": MANDATORY_ABLATIONS,
        "analysis_partition_for_label_dependent_results": "train",
        "test_partition_policy": "Structural integrity only; test labels were not read for association analyses",
        "frozen_source_sha256": sha256_file(frozen_path),
        "cohort": {"frozen_rows": len(frozen), "modelling_rows": len(full)},
        "locked_splits": {
            split: {"variants": values["rows"], "genes": len(values["genes"])}
            for split, values in split_structures.items()
        },
        "partition_overlap": {key: len(value) for key, value in overlaps.items()},
        "approved_primary_features": primary,
        "sensitivity_only_features": sorted(row["field"] for row in lineage if row["approved_sensitivity"] == "1"),
        "prohibited_features": sorted(row["field"] for row in lineage if row["prohibited"] == "1"),
        "class_associated_missingness_fields": [row["field"] for row in missing_flags],
        "redundancy_candidate_count": len(redundancy),
        "known_scientific_risks": [
            "ClinVar classifications may incorporate population and computational evidence.",
            "Optional computational predictors may have training-data or methodological overlap.",
            "Gene-constraint values may act as gene-level proxies.",
            "UniProt curated-variant tokens are prohibited from primary inputs.",
            "Population missingness may act as a label shortcut.",
            "Conservation measurements may repeat related evolutionary signals.",
            "InterPro and UniProt context indicators may overlap.",
        ],
    }
    write_json(output_dir / "15_static_circularity_audit.json", report)
    return report


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser()
    value.add_argument("--frozen", default=str(PROJECT_ROOT / "data/processed/final_dataset_full.csv"))
    value.add_argument("--full", default=str(PROJECT_ROOT / "data/processed/modeling_ready_full.csv"))
    value.add_argument("--core", default=str(PROJECT_ROOT / "data/processed/modeling_ready_core.csv"))
    value.add_argument("--no-population", default=str(PROJECT_ROOT / "data/processed/modeling_ready_core_no_population.csv"))
    value.add_argument("--registry", default=str(PROJECT_ROOT / "data/processed/modeling_feature_registry.csv"))
    value.add_argument("--lineage", default=str(PROJECT_ROOT / "data/processed/evidence_lineage_registry.csv"))
    value.add_argument("--assignments", default=str(PROJECT_ROOT / "data/processed/modeling_split_gene_assignments.csv"))
    value.add_argument("--split-dir", default=str(PROJECT_ROOT / "data/processed/modeling_splits"))
    value.add_argument("--output-dir", default=str(PROJECT_ROOT / "audits/circularity"))
    return value


def main() -> None:
    args = parser().parse_args()
    report = audit(args)
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["critical_failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
