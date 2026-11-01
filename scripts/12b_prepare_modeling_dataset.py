#!/usr/bin/env python3
"""Create leakage-controlled modelling tables from the frozen full dataset."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    read_csv,
    sha256_file,
    write_csv,
    write_json,
)


INTERPRO_TOKENS = (
    "active_site",
    "binding_site",
    "conserved_site",
    "domain",
    "family",
    "homologous_superfamily",
    "ptm",
    "repeat",
)

UNIPROT_EXCLUDED_TOKENS = {
    "Natural variant",
    "Mutagenesis",
    "Sequence conflict",
    "Alternative sequence",
}

BASE_MODEL_FEATURES = (
    "gnomad_global_af",
    "gnomad_max_population_af",
    "conservation_gerp_rs",
    "conservation_phylop100way",
    "conservation_phastcons100way",
    "gnomad_missense_z",
    "gnomad_loeuf",
    "aa_grantham_score",
    "aa_composition_delta",
    "aa_polarity_delta",
    "aa_volume_delta",
    "aa_hydropathy_delta",
    "aa_charge_delta",
    "aa_charge_changed",
    "aa_polarity_class_changed",
    "aa_relative_protein_position",
)

POPULATION_FEATURES = {"gnomad_global_af", "gnomad_max_population_af"}


def tokens(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def gene_matches(row: dict[str, str]) -> bool:
    return row.get("gene_symbol", "").strip() == row.get(
        "reconciled_gene_symbol", ""
    ).strip()


def discover_uniprot_tokens(rows: list[dict[str, str]]) -> list[str]:
    observed: set[str] = set()
    for row in rows:
        observed.update(tokens(row.get("uniprot_overlapping_feature_types", "")))
    return sorted(observed - UNIPROT_EXCLUDED_TOKENS)


def add_context_indicators(
    row: dict[str, str], uniprot_tokens: list[str]
) -> dict[str, Any]:
    derived: dict[str, Any] = dict(row)
    interpro_status = row.get("interpro_lookup_status", "")
    interpro_values = tokens(row.get("interpro_overlap_types", ""))
    for token in INTERPRO_TOKENS:
        column = f"interpro_type_{slug(token)}"
        derived[column] = "" if interpro_status != "FOUND" else int(token in interpro_values)
    uniprot_values = tokens(row.get("uniprot_overlapping_feature_types", ""))
    for token in uniprot_tokens:
        derived[f"uniprot_feature_{slug(token)}"] = int(token in uniprot_values)
    return derived


def registry_rows(
    model_features: list[str], no_population_features: set[str]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for feature in model_features:
        if feature == "aa_relative_protein_position":
            group, data_type, interpretation = (
                "protein_context",
                "numeric",
                "relative residue position within the protein",
            )
        elif feature in {"aa_charge_changed", "aa_polarity_class_changed"}:
            group, data_type, interpretation = (
                "amino_acid_change",
                "binary",
                "amino-acid physicochemical class-change indicator",
            )
        elif feature.startswith("interpro_type_"):
            group, data_type, interpretation = (
                "protein_context",
                "binary",
                "normalized InterPro overlap-type indicator; blank means lookup failure",
            )
        elif feature.startswith("uniprot_feature_"):
            group, data_type, interpretation = (
                "protein_context",
                "binary",
                "normalized UniProt structural/context feature indicator",
            )
        elif feature in POPULATION_FEATURES:
            group, data_type, interpretation = (
                "population",
                "numeric",
                "population frequency; unavailable values stay blank and receive no missingness indicator",
            )
        elif feature.startswith("conservation_"):
            group, data_type, interpretation = "conservation", "numeric", "conservation evidence"
        elif feature.startswith("gnomad_"):
            group, data_type, interpretation = "gene_constraint", "numeric", "gene constraint evidence"
        elif feature.startswith("aa_"):
            group, data_type, interpretation = "amino_acid_change", "numeric", "amino-acid change evidence"
        else:
            group, data_type, interpretation = "other", "numeric", "derived model feature"
        rows.append(
            {
                "field": feature,
                "evidence_group": group,
                "data_type": data_type,
                "primary_model_feature": int(feature not in POPULATION_FEATURES),
                "no_population_ablation_feature": int(feature in no_population_features),
                "interpretation": interpretation,
            }
        )
    for token in sorted(UNIPROT_EXCLUDED_TOKENS):
        rows.append(
            {
                "field": f"uniprot_feature_{slug(token)}",
                "evidence_group": "protein_context",
                "data_type": "binary",
                "primary_model_feature": 0,
                "no_population_ablation_feature": 0,
                "interpretation": "excluded from primary models because it may encode curated variant knowledge",
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        default=str(PROJECT_ROOT / "data/processed/final_dataset_full.csv"),
    )
    parser.add_argument(
        "--expected-source-sha256",
        default="50d8ac02ededae696d6867c7e713351d33058d01bad4fab3633055157525494f",
    )
    args = parser.parse_args()
    source = Path(args.input)
    source_hash = sha256_file(source)
    if args.expected_source_sha256 and source_hash != args.expected_source_sha256:
        raise RuntimeError(
            "Frozen source checksum mismatch: "
            f"expected {args.expected_source_sha256}, observed {source_hash}"
        )
    rows = read_csv(source)
    excluded = [row for row in rows if not gene_matches(row)]
    retained = [row for row in rows if gene_matches(row)]
    uniprot_tokens = discover_uniprot_tokens(retained)
    derived = [add_context_indicators(row, uniprot_tokens) for row in retained]

    interpro_columns = [f"interpro_type_{slug(token)}" for token in INTERPRO_TOKENS]
    uniprot_columns = [f"uniprot_feature_{slug(token)}" for token in uniprot_tokens]
    context_columns = interpro_columns + uniprot_columns
    full_columns = list(rows[0]) + context_columns
    full_path = PROJECT_ROOT / "data/processed/modeling_ready_full.csv"
    write_csv(full_path, derived, full_columns)

    metadata = ["variant_key", "reconciled_gene_symbol", "binary_label"]
    model_features = list(BASE_MODEL_FEATURES) + context_columns
    core_columns = metadata + model_features
    core_rows = [{column: row.get(column, "") for column in core_columns} for row in derived]
    core_path = PROJECT_ROOT / "data/processed/modeling_ready_core.csv"
    write_csv(core_path, core_rows, core_columns)

    no_population_features = set(model_features) - POPULATION_FEATURES
    no_population_columns = metadata + [
        feature for feature in model_features if feature in no_population_features
    ]
    no_population_path = PROJECT_ROOT / "data/processed/modeling_ready_core_no_population.csv"
    write_csv(
        no_population_path,
        ({column: row.get(column, "") for column in no_population_columns} for row in derived),
        no_population_columns,
    )

    registry = registry_rows(model_features, no_population_features)
    registry_path = PROJECT_ROOT / "data/processed/modeling_feature_registry.csv"
    write_csv(registry_path, registry, list(registry[0]))

    excluded_path = PROJECT_ROOT / "data/processed/excluded_alternate_gene_variants.csv"
    excluded_columns = [
        "variant_key",
        "gene_symbol",
        "reconciled_gene_symbol",
        "binary_label",
        "transcript_selection_basis",
    ]
    write_csv(excluded_path, excluded, excluded_columns)

    zero_overlap = sum(
        row.get("interpro_lookup_status") == "FOUND"
        and not tokens(row.get("interpro_overlap_types", ""))
        for row in retained
    )
    lookup_failures = sum(row.get("interpro_lookup_status") != "FOUND" for row in retained)
    report = {
        "status": "PASS",
        "frozen_source_sha256": source_hash,
        "source": audit_basic_rows(rows),
        "derived": audit_basic_rows(derived),
        "alternate_gene_rows_excluded": len(excluded),
        "excluded_class_counts": dict(Counter(row["binary_label"] for row in excluded)),
        "interpro": {
            "normalized_tokens": list(INTERPRO_TOKENS),
            "successful_lookup_zero_overlap": zero_overlap,
            "lookup_failures": lookup_failures,
            "failure_encoding": "blank token indicators; successful zero-overlap is encoded as zeros",
        },
        "uniprot": {
            "included_normalized_tokens": uniprot_tokens,
            "excluded_leakage_risk_tokens": sorted(UNIPROT_EXCLUDED_TOKENS),
        },
        "population_evidence_policy": {
            "missingness_indicators_in_primary_table": False,
            "unavailable_values": "preserved as blank",
            "model_requirement": (
                "Fit imputation and scaling on training genes only. Treat imputed population values "
                "as neutral and report the evidence as unavailable. Evaluate the supplied no-population ablation."
            ),
        },
        "outputs": [str(path) for path in (full_path, core_path, no_population_path, registry_path, excluded_path)],
    }
    audit_path = PROJECT_ROOT / "audits/12b_modeling_preprocessing_audit.json"
    write_json(audit_path, report)
    append_stage_manifest(
        "12b_prepare_modeling_dataset",
        [source],
        [full_path, core_path, no_population_path, registry_path, excluded_path, audit_path],
        {"expected_source_sha256": args.expected_source_sha256},
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
