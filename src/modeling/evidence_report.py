"""Deterministic biological evidence reports for grouped-additive predictions."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .common import FrozenTrainingData, load_frozen_training_data


GROUP_LABELS = {
    "conservation": "Conservation",
    "gene_constraint": "Gene constraint",
    "amino_acid_properties": "Amino-acid properties",
    "protein_context": "Protein context",
}


def _format_value(value: Any) -> str:
    return "unavailable" if pd.isna(value) else f"{float(value):.3g}"


def _feature_label(name: str) -> str:
    for prefix in ("interpro_type_", "uniprot_feature_"):
        if name.startswith(prefix):
            source = "InterPro" if prefix.startswith("interpro") else "UniProt"
            return f"{source} {name.removeprefix(prefix).replace('_', ' ')}"
    return name.replace("_", " ")


def observed_evidence(row: pd.Series, group: str, features: tuple[str, ...]) -> str:
    if group == "conservation":
        return (
            f"phyloP={_format_value(row['conservation_phylop100way'])}; "
            f"PhastCons={_format_value(row['conservation_phastcons100way'])}"
        )
    if group == "gene_constraint":
        return (
            f"missense Z={_format_value(row['gnomad_missense_z'])}; "
            f"LOEUF={_format_value(row['gnomad_loeuf'])}"
        )
    if group == "amino_acid_properties":
        charge_changed = "unavailable" if pd.isna(row["aa_charge_changed"]) else ("yes" if row["aa_charge_changed"] else "no")
        polarity_changed = "unavailable" if pd.isna(row["aa_polarity_class_changed"]) else ("yes" if row["aa_polarity_class_changed"] else "no")
        return (
            f"Grantham={_format_value(row['aa_grantham_score'])}; "
            f"hydropathy change={_format_value(row['aa_hydropathy_delta'])}; "
            f"charge change={_format_value(row['aa_charge_delta'])}; "
            f"charge class changed={charge_changed}; polarity class changed={polarity_changed}"
        )
    present = [_feature_label(feature) for feature in features if not pd.isna(row[feature]) and float(row[feature]) == 1.0]
    missing = sum(pd.isna(row[feature]) for feature in features)
    if present:
        text = "overlapping annotations: " + ", ".join(present)
    else:
        text = "no registered InterPro or UniProt context annotation overlapped the residue"
    if missing:
        text += f"; {missing} context fields unavailable"
    return text


def contribution_thresholds(predictions: pd.DataFrame, groups: list[str]) -> dict[str, dict[str, float]]:
    thresholds: dict[str, dict[str, float]] = {}
    for group in groups:
        values = predictions[f"contribution_{group}"].abs()
        thresholds[group] = {
            "weak_boundary": float(values.quantile(0.25)),
            "strong_boundary": float(values.quantile(0.75)),
        }
    return thresholds


def interpret_contribution(value: float, weak: float, strong: float) -> tuple[str, str]:
    magnitude = abs(value)
    strength = "weakly" if magnitude <= weak else ("strongly" if magnitude >= strong else "moderately")
    direction = "pathogenicity" if value > 0 else "benignity"
    if np.isclose(value, 0.0, atol=1e-12):
        return "negligible", "neither class"
    return strength, direction


def agreement_statement(contributions: dict[str, float], thresholds: dict[str, dict[str, float]]) -> str:
    material = {
        group: value
        for group, value in contributions.items()
        if abs(value) > thresholds[group]["weak_boundary"]
    }
    directions = {"pathogenic" if value > 0 else "benign" for value in material.values()}
    if not material:
        return "All evidence-group effects were weak for this prediction."
    if len(directions) == 1:
        return f"The material evidence-group contributions agreed in supporting {next(iter(directions))}ity."
    positive = [GROUP_LABELS[group] for group, value in material.items() if value > 0]
    negative = [GROUP_LABELS[group] for group, value in material.items() if value < 0]
    return (
        "The model contained conflicting evidence: "
        + ", ".join(positive)
        + " supported pathogenicity, while "
        + ", ".join(negative)
        + " supported benignity."
    )


def create_reports(
    data: FrozenTrainingData,
    predictions: pd.DataFrame,
    selected_threshold: float,
) -> tuple[pd.DataFrame, dict[str, dict[str, float]]]:
    groups = list(data.evidence_groups)
    thresholds = contribution_thresholds(predictions, groups)
    raw = data.frame.drop(columns=["cv_fold"]).copy()
    merged = predictions.merge(raw, on=[data.variant_key, data.gene, data.label], how="left", validate="one_to_one")
    report_rows: list[dict[str, Any]] = []
    for _, row in merged.iterrows():
        predicted = "pathogenic" if row["probability"] >= selected_threshold else "benign"
        contributions = {group: float(row[f"contribution_{group}"]) for group in groups}
        result: dict[str, Any] = {
            data.variant_key: row[data.variant_key],
            data.gene: row[data.gene],
            "observed_label": int(row[data.label]),
            "oof_probability": float(row["probability"]),
            "predicted_class": predicted,
            "classification_threshold": selected_threshold,
        }
        sentences = [
            f"The training out-of-fold model predicted {predicted} with probability {row['probability']:.3f}."
        ]
        missing_total = 0
        for group, features in data.evidence_groups.items():
            value = contributions[group]
            strength, direction = interpret_contribution(
                value,
                thresholds[group]["weak_boundary"],
                thresholds[group]["strong_boundary"],
            )
            observation = observed_evidence(row, group, features)
            missing = int(row[list(features)].isna().sum())
            missing_total += missing
            result[f"{group}_contribution"] = value
            result[f"{group}_interpretation"] = f"{strength} supports {direction}" if strength != "negligible" else "negligible contribution"
            result[f"{group}_observations"] = observation
            result[f"{group}_missing_fields"] = missing
            sentences.append(
                f"{GROUP_LABELS[group]} {result[f'{group}_interpretation']} ({value:+.3f} logit); observed annotations: {observation}."
            )
        agreement = agreement_statement(contributions, thresholds)
        sentences.append(agreement)
        if missing_total:
            sentences.append(
                f"{missing_total} of the 47 primary evidence fields were unavailable and were not described as biological absence."
            )
        sentences.append(
            "These statements describe the model's internal evidence contributions and do not establish causality or independent clinical evidence."
        )
        result["evidence_relationship"] = agreement
        result["total_missing_primary_fields"] = missing_total
        result["narrative"] = " ".join(sentences)
        report_rows.append(result)
    return pd.DataFrame(report_rows), thresholds


def run_evidence_reporting(
    repo_root: str | Path,
    *,
    data_root: str | Path | None = None,
    prediction_path: str | Path = "outputs/grouped_additive_model/02_grouped_additive/oof_predictions_and_group_contributions.csv",
    summary_path: str | Path = "outputs/grouped_additive_model/02_grouped_additive/run_summary.json",
    output_dir: str | Path = "outputs/grouped_additive_model/03_evidence_reports",
) -> dict[str, Any]:
    repo_root = Path(repo_root).resolve()
    data = load_frozen_training_data(repo_root, data_root=data_root)
    predictions = pd.read_csv(repo_root / prediction_path)
    summary = json.loads((repo_root / summary_path).read_text())
    reports, thresholds = create_reports(data, predictions, float(summary["selected_threshold"]))
    output_dir = repo_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    reports.to_csv(output_dir / "training_oof_evidence_reports.csv", index=False)
    with (output_dir / "contribution_interpretation_thresholds.json").open("w") as stream:
        json.dump(thresholds, stream, indent=2, sort_keys=True)
        stream.write("\n")
    examples = reports.sort_values("oof_probability").iloc[[0, len(reports) // 2, -1]]
    lines = ["# Example training out-of-fold evidence reports", ""]
    for _, row in examples.iterrows():
        lines.extend([f"## {row[data.variant_key]} ({row[data.gene]})", "", row["narrative"], ""])
    (output_dir / "example_evidence_reports.md").write_text("\n".join(lines), encoding="utf-8")
    metadata = {
        "status": "training_only_complete",
        "reports": len(reports),
        "deterministic_template": True,
        "llm_used_for_generation": False,
        "validation_loaded": False,
        "test_loaded": False,
        "caution": "model contribution is not causality or independent clinical evidence",
    }
    with (output_dir / "run_metadata.json").open("w") as stream:
        json.dump(metadata, stream, indent=2, sort_keys=True)
        stream.write("\n")
    return metadata
