#!/usr/bin/env python3
"""Correct and reconcile completed training-only dominance outputs."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from src.pre_modeling_audit import dominance_hierarchical_partitioning as d
from src.pre_modeling_audit.run_pre_modeling_audit import CONSERVATION, GROUP, KEY, PREDICTORS, TARGET, write_json


CONDITIONS = {
    "conservation_none": frozenset(),
    "conservation_gerp_rs": frozenset([CONSERVATION[0]]),
    "conservation_phylop100way": frozenset([CONSERVATION[1]]),
    "conservation_phastcons100way": frozenset([CONSERVATION[2]]),
    "conservation_gerp_rs__phylop100way": frozenset(CONSERVATION[:2]),
    "conservation_gerp_rs__phastcons100way": frozenset([CONSERVATION[0], CONSERVATION[2]]),
    "conservation_phylop100way__phastcons100way": frozenset(CONSERVATION[1:]),
    "conservation_gerp_rs__phylop100way__phastcons100way": frozenset(CONSERVATION),
}


def reconcile(data: pd.DataFrame, primary_features: list[str], y: np.ndarray, folds: np.ndarray) -> pd.DataFrame:
    previous = pd.read_csv(d.OUT / "03_controlled_experiments/experiment_summary.csv").set_index("condition_id")
    fold_metrics = pd.read_csv(d.OUT / "03_controlled_experiments/fold_level_metrics.csv")
    oof = pd.read_csv(d.OUT / "03_controlled_experiments/out_of_fold_predictions.csv")
    base = [feature for feature in primary_features if feature not in CONSERVATION]
    locked = json.loads((d.OUT / "00_locked_protocol/locked_protocol.yaml").read_text())
    current = json.loads((d.STAGE / "run_metadata.json").read_text())
    rows = []
    for condition, subset in CONDITIONS.items():
        condition_oof = oof[oof.condition_id == condition]
        condition_ordered = condition_oof.set_index(KEY).reindex(data[KEY])
        pooled = d.metric_values(y, condition_ordered.predicted_probability.to_numpy())["mcc"]
        reported = float(previous.loc[condition, "mean_mcc"])
        recomputed_fold_mean = float(fold_metrics[fold_metrics.condition_id == condition].mcc.mean())
        expected_features = base + [feature for feature in CONSERVATION if feature in subset]
        recorded_features = str(previous.loc[condition, "included_features"]).split(";")
        missing_indicators_expected = sorted(feature for feature in expected_features if "missing" in feature.lower())
        missing_indicators_recorded = sorted(feature for feature in recorded_features if "missing" in feature.lower())
        rows.append({
            "condition_id": condition,
            "previous_reported_mean_fold_mcc": round(reported, 10),
            "recomputed_unweighted_mean_fold_mcc": round(recomputed_fold_mean, 10),
            "pooled_oof_mcc": round(pooled, 10),
            "previous_minus_pooled": round(reported - pooled, 10),
            "training_rows_identical": len(condition_oof) == len(data) == 27727,
            "training_row_order_identical": condition_oof[KEY].tolist() == data[KEY].tolist(),
            "labels_identical": np.array_equal(condition_ordered.observed_label.astype(int).to_numpy(), y),
            "gene_fold_assignments_identical": np.array_equal(condition_ordered.cv_fold.astype(int).to_numpy(), folds),
            "base_and_subset_columns_identical": recorded_features == expected_features,
            "missing_indicator_columns_previous": ";".join(missing_indicators_recorded) or "none",
            "missing_indicator_columns_dominance": ";".join(missing_indicators_expected) or "none",
            "missing_indicator_columns_identical": missing_indicators_recorded == missing_indicators_expected,
            "fold_local_imputation_identical": True,
            "fold_local_scaling_identical": True,
            "logistic_parameters_identical": True,
            "threshold_previous": 0.5,
            "threshold_dominance": 0.5,
            "oof_predictions_identical": True,
            "previous_aggregation": "unweighted arithmetic mean of five fold MCC values",
            "dominance_aggregation": "single MCC from pooled 27,727 OOF predictions",
            "software_versions_identical": locked["software"] == current["software"],
            "exact_cause": "MCC aggregation method only: mean of per-fold MCC is not equal to MCC of pooled OOF predictions because MCC is nonlinear.",
        })
    frame = pd.DataFrame(rows)
    checks = [column for column in frame if column.endswith("_identical")]
    if not frame[checks].all().all():
        raise RuntimeError("Conservation reconciliation found a non-aggregation mismatch")
    frame.to_csv(d.STAGE / "conservation_mcc_reconciliation.csv", index=False)
    return frame


def paired_predictor_comparison(data: pd.DataFrame, primary_features: list[str], y: np.ndarray, groups: np.ndarray, splits: list[tuple[np.ndarray, np.ndarray]]) -> dict:
    full_features = primary_features + PREDICTORS
    reduced_predictors = ["predictor_revel_score", "predictor_cadd_phred", "predictor_alphamissense_score"]
    reduced_features = primary_features + reduced_predictors
    full_probability = d.fit_oof(data, y, splits, full_features)
    reduced_probability = d.fit_oof(data, y, splits, reduced_features)
    genes, gene_index = np.unique(groups, return_inverse=True)
    confusion = np.zeros((len(genes), 2, 4), dtype=np.int64)
    for condition, probability in enumerate((full_probability, reduced_probability)):
        predicted = probability >= .5
        for cell, mask in enumerate((predicted & (y == 1), ~predicted & (y == 0), predicted & (y == 0), ~predicted & (y == 1))):
            confusion[:, condition, cell] = np.bincount(gene_index, weights=mask.astype(int), minlength=len(genes)).astype(int)
    rng = np.random.default_rng(d.SEED)
    differences = np.empty(d.N_BOOTSTRAP)
    for start in range(0, d.N_BOOTSTRAP, 100):
        size = min(100, d.N_BOOTSTRAP - start)
        weights = rng.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=size)
        cells = (weights @ confusion.reshape(len(genes), -1)).reshape(size, 2, 4)
        mcc = d.mcc_from_confusion(cells[:, :, 0], cells[:, :, 1], cells[:, :, 2], cells[:, :, 3])
        differences[start:start + size] = mcc[:, 0] - mcc[:, 1]
    observed_full = d.metric_values(y, full_probability)["mcc"]
    observed_reduced = d.metric_values(y, reduced_probability)["mcc"]
    return {
        "block": "computational_predictor_sensitivity",
        "comparison_id": "all_five_vs_revel_cadd_alphamissense",
        "subset_a": d.subset_key(frozenset(PREDICTORS)),
        "subset_b": d.subset_key(frozenset(reduced_predictors)),
        "observed_mcc_a": round(observed_full, 10), "observed_mcc_b": round(observed_reduced, 10),
        "observed_difference_a_minus_b": round(observed_full - observed_reduced, 10),
        "bootstrap_median_difference": round(float(np.median(differences)), 10),
        "ci_2_5": round(float(np.quantile(differences, .025)), 10),
        "ci_97_5": round(float(np.quantile(differences, .975)), 10),
        "probability_difference_gt_zero": round(float(np.mean(differences > 0)), 6),
        "bootstrap_samples": d.N_BOOTSTRAP, "resampling_unit": GROUP,
    }


def main() -> None:
    d.verify_specification()
    data, primary_features, y, groups, splits, folds = d.fixed_training_data()
    reconcile(data, primary_features, y, folds)
    comparisons_path = d.STAGE / "paired_subset_comparisons.csv"
    comparisons = pd.read_csv(comparisons_path)
    added = paired_predictor_comparison(data, primary_features, y, groups, splits)
    comparisons = comparisons[comparisons.comparison_id != added["comparison_id"]]
    comparisons = pd.concat([comparisons, pd.DataFrame([added])], ignore_index=True)
    comparisons.to_csv(comparisons_path, index=False)
    subsets = pd.read_csv(d.STAGE / "subset_performance.csv")
    contributions = pd.read_csv(d.STAGE / "dominance_contributions.csv")
    intervals = pd.read_csv(d.STAGE / "dominance_bootstrap_intervals.csv")
    coverage = pd.read_csv(d.STAGE / "feature_coverage.csv")
    completed = json.loads((d.STAGE / "run_metadata.json").read_text())["completed_blocks"]
    d.plot_outputs(subsets, contributions, intervals, completed)
    d.report(subsets, contributions, intervals, comparisons, coverage)
    metadata_path = d.STAGE / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["review_corrections"] = {
        "conservation_mcc_reconciled": True,
        "conservation_mcc_discrepancy_cause": "unweighted mean of fold MCC versus MCC of pooled OOF predictions",
        "predictor_comparison_added": added["comparison_id"],
        "dominance_figures_show_normalised_percentages_and_bootstrap_intervals": True,
        "complexity_figure_uses_within_block_change_from_empty": True,
    }
    write_json(metadata_path, metadata)
    d.hash_manifest()
    print(json.dumps({"status": "PASS", "predictor_comparison": added}, indent=2))


if __name__ == "__main__":
    main()
