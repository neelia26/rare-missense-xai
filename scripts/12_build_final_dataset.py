#!/usr/bin/env python3
"""Apply documented eligibility rules and create balanced full/core/optional tables."""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    load_config,
    read_csv,
    write_csv,
    write_json,
)


FEATURE_REGISTRY = [
    ("gnomad_global_af", "population", "numeric", True, "global allele frequency"),
    ("gnomad_max_population_af", "population", "numeric", True, "maximum population allele frequency"),
    ("gnomad_direct_exome_ac", "population", "numeric", False, "optional direct count"),
    ("gnomad_direct_exome_an", "population", "numeric", False, "optional direct count denominator"),
    ("gnomad_direct_exome_nhomalt", "population", "numeric", False, "optional homozygote count"),
    ("gnomad_direct_genome_ac", "population", "numeric", False, "optional direct count"),
    ("gnomad_direct_genome_an", "population", "numeric", False, "optional direct count denominator"),
    ("gnomad_direct_genome_nhomalt", "population", "numeric", False, "optional homozygote count"),
    ("conservation_gerp_rs", "conservation", "numeric", True, "GERP conservation score"),
    ("conservation_phylop100way", "conservation", "numeric", True, "PhyloP conservation score"),
    ("conservation_phastcons100way", "conservation", "numeric", True, "PhastCons conservation score"),
    ("gnomad_missense_z", "gene_constraint", "numeric", True, "missense constraint"),
    ("gnomad_loeuf", "gene_constraint", "numeric", True, "loss-of-function constraint context"),
    ("aa_grantham_score", "amino_acid_change", "numeric", True, "Grantham distance"),
    ("aa_composition_delta", "amino_acid_change", "numeric", True, "composition change"),
    ("aa_polarity_delta", "amino_acid_change", "numeric", True, "polarity change"),
    ("aa_volume_delta", "amino_acid_change", "numeric", True, "molecular-volume change"),
    ("aa_hydropathy_delta", "amino_acid_change", "numeric", True, "hydropathy change"),
    ("aa_charge_delta", "amino_acid_change", "numeric", True, "approximate charge change"),
    ("aa_charge_changed", "amino_acid_change", "binary", True, "charge-class change"),
    ("aa_polarity_class_changed", "amino_acid_change", "binary", True, "polarity-class change"),
    ("aa_relative_protein_position", "protein_context", "numeric", True, "relative residue position"),
    ("interpro_overlap_count", "protein_context", "numeric", True, "number of overlapping InterPro entries"),
    ("interpro_overlap_types", "protein_context", "categorical", True, "overlapping InterPro entry types"),
    ("uniprot_overlapping_feature_types", "protein_context", "categorical", True, "overlapping UniProt feature types"),
]


def effective_rarity(row: dict[str, str]) -> str:
    direct = row.get("gnomad_direct_rarity_status", "")
    if direct:
        return direct
    return row.get("gnomad_rarity_status", "UNKNOWN")


def exclusion_reason(row: dict[str, str], allow_not_reported: bool) -> str:
    if row.get("binary_label") not in {"0", "1"}:
        return "INVALID_LABEL"
    if row.get("strict_missense_eligible") not in {"1", 1}:
        return "TRANSCRIPT_NOT_STRICT_MISSENSE"
    if row.get("uniprot_validation_status") != "MATCH":
        return "UNIPROT_REFERENCE_NOT_VALIDATED"
    if not (row.get("reconciled_gene_symbol") or row.get("gene_symbol")):
        return "MISSING_GROUP_GENE"
    rarity = effective_rarity(row)
    if rarity == "NOT_RARE":
        return "NOT_RARE"
    if rarity == "NOT_REPORTED" and not allow_not_reported:
        return "FREQUENCY_NOT_REPORTED"
    if rarity not in {"RARE", "NOT_REPORTED"}:
        return "POPULATION_LOOKUP_UNRESOLVED"
    return "ELIGIBLE"


def sample_balanced(
    rows: list[dict[str, str]], target: int, seed: int, allow_smaller: bool
) -> tuple[list[dict[str, str]], int]:
    by_label = {
        "0": [row for row in rows if row["binary_label"] == "0"],
        "1": [row for row in rows if row["binary_label"] == "1"],
    }
    achievable = min(len(by_label["0"]), len(by_label["1"]))
    actual = min(target, achievable) if allow_smaller else target
    if achievable < target and not allow_smaller:
        raise RuntimeError(
            f"Cannot select {target:,} per class after governance checks: "
            f"benign={len(by_label['0']):,}, pathogenic={len(by_label['1']):,}. "
            "Use --allow-smaller for a pilot; do not silently lower the full target."
        )
    rng = random.Random(seed)
    selected = rng.sample(by_label["0"], actual) + rng.sample(by_label["1"], actual)
    rng.shuffle(selected)
    return selected, actual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/11_amino_acid_features.csv"))
    parser.add_argument("--optional-predictors", default=str(PROJECT_ROOT / "data/interim/10_optional_predictor_scores.csv"))
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    parser.add_argument("--target-per-class", type=int)
    parser.add_argument("--allow-smaller", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    input_path = Path(args.input)
    rows = read_csv(input_path)
    allow_not_reported = bool(config["rarity"]["include_not_reported_after_successful_lookup"])
    reasons = Counter(exclusion_reason(row, allow_not_reported) for row in rows)
    eligible = [row for row in rows if exclusion_reason(row, allow_not_reported) == "ELIGIBLE"]
    target = args.target_per_class or int(config["sampling"]["final_per_class"])
    seed = int(config["project"]["random_seed"])
    selected, actual = sample_balanced(eligible, target, seed, args.allow_smaller)
    selected.sort(key=lambda row: row["variant_key"])
    full_path = PROJECT_ROOT / "data/processed/final_dataset_full.csv"
    write_csv(full_path, selected, list(rows[0]))

    registry_rows: list[dict[str, Any]] = []
    available_core_features: list[str] = []
    for feature, group, data_type, default_use, interpretation in FEATURE_REGISTRY:
        available = feature in rows[0]
        if available and default_use:
            available_core_features.append(feature)
        registry_rows.append(
            {
                "field": feature,
                "evidence_group": group,
                "data_type": data_type,
                "available": int(available),
                "default_core_model_feature": int(available and default_use),
                "interpretation": interpretation,
            }
        )
    registry_path = PROJECT_ROOT / "data/processed/evidence_feature_registry.csv"
    write_csv(registry_path, registry_rows, list(registry_rows[0]))

    metadata = ["variant_key", "reconciled_gene_symbol", "binary_label"]
    core_columns = metadata + available_core_features
    core_rows = [{column: row.get(column, "") for column in core_columns} for row in selected]
    core_path = PROJECT_ROOT / "data/processed/model_table_core.csv"
    write_csv(core_path, core_rows, core_columns)

    optional_path = Path(args.optional_predictors)
    experimental_path = PROJECT_ROOT / "data/processed/model_table_with_optional_predictors.csv"
    experimental_created = False
    if optional_path.exists():
        optional_rows = read_csv(optional_path)
        by_key = {row["variant_key"]: row for row in optional_rows}
        predictor_columns = [column for column in optional_rows[0] if column.startswith("predictor_")]
        experimental_rows = []
        for row in core_rows:
            combined = dict(row)
            supplement = by_key.get(row["variant_key"], {})
            combined.update({column: supplement.get(column, "") for column in predictor_columns})
            experimental_rows.append(combined)
        write_csv(experimental_path, experimental_rows, core_columns + predictor_columns)
        experimental_created = True

    audit = {
        **audit_basic_rows(selected),
        "input_rows": len(rows),
        "eligibility_counts": dict(reasons),
        "eligible_before_balancing": dict(Counter(row["binary_label"] for row in eligible)),
        "target_per_class": target,
        "actual_per_class": actual,
        "random_seed": seed,
        "rarity_policy": {
            "reported_af_threshold": config["rarity"]["maximum_global_af"],
            "include_not_reported_after_successful_lookup": allow_not_reported,
            "not_reported_is_never_replaced_with_zero": True,
        },
        "core_model_features": available_core_features,
        "optional_predictor_condition_created": experimental_created,
        "prohibited_as_model_inputs": [
            "ClinVar label text/review/condition fields",
            "variant and database identifiers",
            "gene symbol and transcript identifiers (gene is retained only for grouped splitting)",
            "source, join and lookup-status provenance fields unless a later audit explicitly justifies them",
        ],
    }
    audit_path = PROJECT_ROOT / "audits/12_final_build_audit.json"
    write_json(audit_path, audit)
    outputs = [full_path, core_path, registry_path, audit_path]
    if experimental_created:
        outputs.append(experimental_path)
    append_stage_manifest(
        "12_build_final_dataset",
        [input_path] + ([optional_path] if optional_path.exists() else []),
        outputs,
        {"target_per_class": target, "seed": seed, "allow_smaller": args.allow_smaller},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
