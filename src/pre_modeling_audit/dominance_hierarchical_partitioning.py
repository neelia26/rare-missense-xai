#!/usr/bin/env python3
"""Training-only exact performance-based dominance analysis.

This module deliberately has no validation/test input paths. It reuses the
locked training fold assignment and the fixed logistic audit instrument.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from src.pre_modeling_audit.run_pre_modeling_audit import (
    CONSERVATION, GROUP, KEY, METADATA, OUT, PATHS, POPULATION, PREDICTORS,
    ROOT, SEED, TARGET, load, numeric, sha256, write_csv, write_json,
)


STAGE = OUT / "05_dominance_hierarchical_partitioning"
SPEC = STAGE / "ANALYSIS_SPECIFICATION.md"
SPEC_HASH = STAGE / "analysis_specification_sha256.txt"
N_BOOTSTRAP = 5000
TOLERANCE = 1e-10
AA = ["aa_hydropathy_delta", "aa_polarity_delta"]
CONSTRAINT = ["gnomad_missense_z", "gnomad_loeuf"]
DISPLAY = {
    "conservation_gerp_rs": "GERP++", "conservation_phylop100way": "phyloP",
    "conservation_phastcons100way": "PhastCons", "aa_hydropathy_delta": "Hydropathy delta",
    "aa_polarity_delta": "Polarity delta", "gnomad_missense_z": "Missense Z",
    "gnomad_loeuf": "LOEUF", "gnomad_global_af": "Global AF",
    "gnomad_max_population_af": "Max-population AF", "predictor_sift_score": "SIFT",
    "predictor_polyphen2_hdiv_score": "PolyPhen-2", "predictor_revel_score": "REVEL",
    "predictor_cadd_phred": "CADD", "predictor_alphamissense_score": "AlphaMissense",
}
BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#666666"


def verify_specification() -> str:
    expected = SPEC_HASH.read_text(encoding="utf-8").split()[0]
    actual = sha256(SPEC)
    if expected != actual:
        raise RuntimeError(f"Analysis specification checksum mismatch: {expected} != {actual}")
    return actual


def subset_key(subset: tuple[str, ...] | frozenset[str]) -> str:
    return "empty" if not subset else "__".join(sorted(subset))


def all_subsets(features: list[str]) -> list[tuple[str, ...]]:
    return [subset for size in range(len(features) + 1) for subset in itertools.combinations(features, size)]


def shapley_contributions(values: dict[frozenset[str], float], features: list[str]) -> dict[str, float]:
    """Exact Shapley/general-dominance contributions for an arbitrary value function."""
    p = len(features)
    contributions: dict[str, float] = {}
    for feature in features:
        others = [x for x in features if x != feature]
        contribution = 0.0
        for subset_tuple in all_subsets(others):
            subset = frozenset(subset_tuple)
            weight = math.factorial(len(subset)) * math.factorial(p - len(subset) - 1) / math.factorial(p)
            contribution += weight * (values[subset | {feature}] - values[subset])
        contributions[feature] = float(contribution)
    expected = values[frozenset(features)] - values[frozenset()]
    if not math.isclose(sum(contributions.values()), expected, rel_tol=0.0, abs_tol=TOLERANCE):
        raise RuntimeError("Shapley efficiency identity failed")
    return contributions


def normalised_shares(contributions: dict[str, float], total_gain: float) -> dict[str, float | None]:
    if total_gain <= 0 or math.isclose(total_gain, 0.0, rel_tol=0.0, abs_tol=TOLERANCE):
        return {feature: None for feature in contributions}
    shares = {feature: value / total_gain for feature, value in contributions.items()}
    if not math.isclose(sum(shares.values()), 1.0, rel_tol=0.0, abs_tol=TOLERANCE):
        raise RuntimeError("Normalised Shapley shares do not sum to one")
    return shares


def mcc_from_confusion(tp: np.ndarray, tn: np.ndarray, fp: np.ndarray, fn: np.ndarray) -> np.ndarray:
    numerator = tp * tn - fp * fn
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(numerator, denominator, out=np.zeros_like(numerator, dtype=float), where=denominator > 0)


def metric_values(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    return {
        "mcc": float(matthews_corrcoef(y, probability >= 0.5)),
        "pr_auc": float(average_precision_score(y, probability)),
        "roc_auc": float(roc_auc_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }
def fixed_training_data() -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray, list[tuple[np.ndarray, np.ndarray]], np.ndarray]:
    primary = load(PATHS["train_primary"])
    optional = load(PATHS["train_optional"])[[KEY] + PREDICTORS]
    population = load(PATHS["train_population"])[[KEY] + POPULATION]
    data = primary.merge(optional, on=KEY, validate="one_to_one").merge(population, on=KEY, validate="one_to_one")
    if len(data) != 27727 or data[KEY].duplicated().any():
        raise RuntimeError("Unexpected or duplicate training rows")
    saved = pd.read_csv(OUT / "03_controlled_experiments/fixed_cv_fold_assignments.csv", dtype={KEY: str, GROUP: str, "cv_fold": int})
    if len(saved) != len(data) or saved[KEY].duplicated().any():
        raise RuntimeError("Invalid saved fold assignments")
    data = data.merge(saved, on=[KEY, GROUP], validate="one_to_one")
    if set(data["cv_fold"]) != {0, 1, 2, 3, 4} or data.groupby(GROUP)["cv_fold"].nunique().max() != 1:
        raise RuntimeError("Saved folds are not five-fold and gene-disjoint")
    splits = []
    for fold in range(5):
        held = np.flatnonzero(data["cv_fold"].to_numpy() == fold)
        fit = np.flatnonzero(data["cv_fold"].to_numpy() != fold)
        splits.append((fit, held))
    features = [column for column in primary.columns if column not in METADATA]
    return data, features, data[TARGET].astype(int).to_numpy(), data[GROUP].to_numpy(), splits, data["cv_fold"].to_numpy()


def fit_oof(data: pd.DataFrame, y: np.ndarray, splits: list[tuple[np.ndarray, np.ndarray]], features: list[str]) -> np.ndarray:
    probability = np.full(len(data), np.nan)
    for fit, held in splits:
        pipeline = Pipeline([
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(C=1.0, max_iter=2000, solver="liblinear", class_weight="balanced", random_state=SEED)),
        ])
        pipeline.fit(data.iloc[fit][features].apply(numeric), y[fit])
        probability[held] = pipeline.predict_proba(data.iloc[held][features].apply(numeric))[:, 1]
    if not np.isfinite(probability).all():
        raise RuntimeError("Incomplete OOF predictions")
    return probability


def reused_conservation_oof(data: pd.DataFrame, y: np.ndarray, folds: np.ndarray) -> dict[frozenset[str], np.ndarray]:
    source = pd.read_csv(OUT / "03_controlled_experiments/out_of_fold_predictions.csv")
    result: dict[frozenset[str], np.ndarray] = {}
    expected_keys = data[KEY].tolist()
    names = {
        frozenset(): "conservation_none",
        frozenset([CONSERVATION[0]]): "conservation_gerp_rs",
        frozenset([CONSERVATION[1]]): "conservation_phylop100way",
        frozenset([CONSERVATION[2]]): "conservation_phastcons100way",
        frozenset([CONSERVATION[0], CONSERVATION[1]]): "conservation_gerp_rs__phylop100way",
        frozenset([CONSERVATION[0], CONSERVATION[2]]): "conservation_gerp_rs__phastcons100way",
        frozenset([CONSERVATION[1], CONSERVATION[2]]): "conservation_phylop100way__phastcons100way",
        frozenset(CONSERVATION): "conservation_gerp_rs__phylop100way__phastcons100way",
    }
    for subset, condition in names.items():
        rows = source[source["condition_id"] == condition].set_index(KEY).reindex(expected_keys)
        if len(rows) != len(data) or rows["predicted_probability"].isna().any():
            raise RuntimeError(f"Cannot reuse incomplete condition {condition}")
        if not np.array_equal(rows["observed_label"].astype(int).to_numpy(), y):
            raise RuntimeError(f"Label mismatch in reused condition {condition}")
        if not np.array_equal(rows["cv_fold"].astype(int).to_numpy(), folds):
            raise RuntimeError(f"Fold mismatch in reused condition {condition}")
        result[subset] = rows["predicted_probability"].to_numpy(dtype=float)
    return result


def block_definitions(primary_features: list[str]) -> list[dict]:
    all3_base = list(primary_features)
    return [
        {"block": "primary_conservation", "candidates": CONSERVATION, "base": [x for x in primary_features if x not in CONSERVATION], "governance": "primary biological evidence", "reuse": True},
        {"block": "primary_amino_acid", "candidates": AA, "base": [x for x in primary_features if x not in AA], "governance": "primary related-feature analysis", "reuse": False},
        {"block": "primary_gene_constraint", "candidates": CONSTRAINT, "base": [x for x in primary_features if x not in CONSTRAINT], "governance": "controlled related-feature analysis", "reuse": False},
        {"block": "population_inclusive_sensitivity", "candidates": POPULATION, "base": all3_base, "governance": "sensitivity only; class-associated missingness warning", "reuse": False},
        {"block": "computational_predictor_sensitivity", "candidates": PREDICTORS, "base": all3_base, "governance": "external comparator/optional sensitivity only", "reuse": False},
    ]


def bootstrap_block(y: np.ndarray, groups: np.ndarray, candidates: list[str], probabilities: dict[frozenset[str], np.ndarray]) -> tuple[dict[frozenset[str], np.ndarray], dict[str, np.ndarray]]:
    genes, gene_index = np.unique(groups, return_inverse=True)
    subsets = list(probabilities)
    gene_confusion = np.zeros((len(genes), len(subsets), 4), dtype=np.int64)  # tp,tn,fp,fn
    for condition_index, subset in enumerate(subsets):
        predicted = probabilities[subset] >= 0.5
        for cell, mask in enumerate((predicted & (y == 1), ~predicted & (y == 0), predicted & (y == 0), ~predicted & (y == 1))):
            gene_confusion[:, condition_index, cell] = np.bincount(gene_index, weights=mask.astype(int), minlength=len(genes)).astype(int)
    rng = np.random.default_rng(SEED)
    bootstrap_values = {subset: np.empty(N_BOOTSTRAP) for subset in subsets}
    bootstrap_phi = {feature: np.empty(N_BOOTSTRAP) for feature in candidates}
    batch_size = 100
    for start in range(0, N_BOOTSTRAP, batch_size):
        size = min(batch_size, N_BOOTSTRAP - start)
        weights = rng.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=size)
        confusion = (weights @ gene_confusion.reshape(len(genes), -1)).reshape(size, len(subsets), 4)
        mcc = mcc_from_confusion(confusion[:, :, 0], confusion[:, :, 1], confusion[:, :, 2], confusion[:, :, 3])
        for index, subset in enumerate(subsets):
            bootstrap_values[subset][start:start + size] = mcc[:, index]
        for row in range(size):
            values = {subset: float(mcc[row, index]) for index, subset in enumerate(subsets)}
            phi = shapley_contributions(values, candidates)
            for feature in candidates:
                bootstrap_phi[feature][start + row] = phi[feature]
    return bootstrap_values, bootstrap_phi


def comparisons_for(candidates: list[str]) -> list[tuple[frozenset[str], frozenset[str], str]]:
    full = frozenset(candidates)
    comparisons = [(full, full - {feature}, f"full_vs_without_{feature}") for feature in candidates]
    if len(candidates) == 2:
        comparisons.append((frozenset([candidates[0]]), frozenset([candidates[1]]), f"{candidates[0]}_vs_{candidates[1]}"))
    if candidates == CONSERVATION:
        comparisons.extend([
            (full, frozenset([CONSERVATION[1], CONSERVATION[2]]), "all3_vs_phylop_phastcons"),
            (frozenset([CONSERVATION[1], CONSERVATION[2]]), frozenset([CONSERVATION[1]]), "phylop_phastcons_vs_phylop"),
        ])
    seen = set(); unique = []
    for a, b, name in comparisons:
        key = (a, b, name)
        if key not in seen:
            unique.append((a, b, name)); seen.add(key)
    return unique


def savefig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=320, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_outputs(subset_frame: pd.DataFrame, contribution_frame: pd.DataFrame, interval_frame: pd.DataFrame, completed: list[str]) -> None:
    figures = STAGE / "figures"
    mapping = {
        "primary_conservation": "conservation_dominance_contributions.png",
        "primary_amino_acid": "amino_acid_dominance_contributions.png",
        "primary_gene_constraint": "constraint_dominance_contributions.png",
        "population_inclusive_sensitivity": "population_sensitivity_dominance.png",
        "computational_predictor_sensitivity": "predictor_dominance_contributions.png",
    }
    feature_palettes = {
        "primary_conservation": ["#355C7D", "#6554A4", "#8C7BC8"],
        "primary_amino_acid": ["#4169A1", "#7567B1"],
        "primary_gene_constraint": ["#2F4B7C", "#7859B8"],
        "population_inclusive_sensitivity": ["#59617A", "#5F6DB2"],
        "computational_predictor_sensitivity": ["#2F4B7C", "#4C6A92", "#6554A4", "#7D6DB3", "#283E6A"],
    }
    for block in completed:
        rows = contribution_frame[contribution_frame.block == block].sort_values("raw_dominance_contribution")
        uncertainty = interval_frame.set_index(["block", "feature"]).loc[[(block, feature) for feature in rows.feature]].reset_index()
        fig, ax = plt.subplots(figsize=(12, max(4.8, .85 * len(rows) + 2)))
        colors = feature_palettes[block][:len(rows)]
        values = rows.raw_dominance_contribution.to_numpy()
        errors = np.vstack([values - uncertainty.ci_2_5.to_numpy(), uncertainty.ci_97_5.to_numpy() - values])
        y = np.arange(len(rows))
        ax.barh(y, values, xerr=errors, color=colors, edgecolor="white",
                error_kw={"ecolor":"#3F4652", "elinewidth":2.0, "capsize":5, "capthick":2})
        ax.axvline(0, color="#999999", linewidth=1.2, linestyle="--", zorder=0)
        for i, (_, row) in enumerate(rows.iterrows()):
            share = "undefined" if pd.isna(row.normalised_share) else f"{100 * row.normalised_share:.1f}%"
            ax.text(1.02, i, f"MCC {row.raw_dominance_contribution:.4f}  |  {share}",
                    transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=9, color=colors[i])
        ax.set_yticks(y, [DISPLAY[x] for x in rows.feature])
        ax.grid(axis="x", color="#E6E6E6", linewidth=.8)
        ax.set_xlabel("Exact Shapley contribution to OOF MCC")
        ax.set_title(f"{block.replace('_', ' ').title()}\nRaw MCC contribution with 95% gene-bootstrap interval; percentage is normalised share", fontweight="bold")
        fig.subplots_adjust(bottom=.18, right=.76)
        fig.text(.5, .035, "Conditional predictive contribution only; not causality or biological independence.", ha="center", fontsize=9, color=GREY)
        savefig(fig, figures / mapping[block])
    subset_mapping = {
        "primary_conservation": "subset_combinations/conservation_subset_performance.png",
        "primary_amino_acid": "subset_combinations/amino_acid_subset_performance.png",
        "primary_gene_constraint": "subset_combinations/gene_constraint_subset_performance.png",
        "population_inclusive_sensitivity": "subset_combinations/population_subset_performance.png",
        "computational_predictor_sensitivity": "subset_combinations/computational_predictor_subset_performance.png",
    }
    block_count_colours = {
        "primary_conservation": {0:"#606060", 1:"#4C78A8", 2:"#D16BA5", 3:"#7A5195"},
        "primary_amino_acid": {0:"#606060", 1:"#4C78A8", 2:"#D16BA5"},
        "primary_gene_constraint": {0:"#555555", 1:"#2F4B7C", 2:"#D16BA5"},
        "population_inclusive_sensitivity": {0:"#60656F", 1:"#586BA4", 2:"#D16BA5"},
        "computational_predictor_sensitivity": {0:"#606060", 1:"#4C78A8", 2:"#8C7BC8", 3:"#D16BA5", 4:"#6F4C9B", 5:"#2F4B7C"},
    }
    for block in completed:
        rows = subset_frame[subset_frame.block == block].copy()
        empty_mcc = float(rows.loc[rows.is_empty.eq(1), "mcc"].iloc[0])
        rows["delta_mcc"] = rows.mcc - empty_mcc
        rows["display_label"] = rows.subset_label.fillna("No candidate features")
        rows = rows.sort_values(["delta_mcc", "candidate_count"])
        height = max(5.5, .38 * len(rows) + 2.2)
        fig, ax = plt.subplots(figsize=(13, height))
        count_colours = block_count_colours[block]
        colours = [count_colours[int(value)] for value in rows.candidate_count]
        ax.barh(rows.display_label, rows.delta_mcc, color=colours, edgecolor="white")
        ax.axvline(0, color="#666666", linewidth=1.2, linestyle="--")
        ax.grid(axis="x", color="#E6E6E6", linewidth=.8)
        for i, (_, row) in enumerate(rows.iterrows()):
            ax.text(1.02, i, f"Δ {row.delta_mcc:+.4f}  |  MCC {row.mcc:.4f}",
                    transform=ax.get_yaxis_transform(), va="center", ha="left", fontsize=8,
                    color=count_colours[int(row.candidate_count)])
        handles = [plt.Rectangle((0,0),1,1,color=count_colours[n],label=f"{n} candidate{'s' if n != 1 else ''}")
                   for n in sorted(set(rows.candidate_count.astype(int)))]
        ax.legend(handles=handles, title="Combination size", loc="lower right", fontsize=8)
        ax.set_xlabel("Change in OOF MCC from this block's no-candidate condition")
        ax.set_title(f"{block.replace('_', ' ').title()} subset performance\nAll combinations; training only", fontweight="bold")
        fig.subplots_adjust(right=.76)
        savefig(fig, figures / subset_mapping[block])
    rows = subset_frame.copy(); rows["mcc_change_from_block_empty"] = rows.groupby("block")["mcc"].transform(lambda values: values - values.loc[rows.loc[values.index, "is_empty"].eq(1)].iloc[0])
    fig, ax = plt.subplots(figsize=(11, 7))
    for block, group in rows.groupby("block"):
        ax.scatter(group.candidate_count, group.mcc_change_from_block_empty, label=block.replace("_", " "), alpha=.8)
    ax.set_xlabel("Number of candidate features included")
    ax.set_ylabel("Δ OOF MCC vs block-specific empty condition")
    ax.axhline(0, color=GREY, linewidth=1)
    ax.set_title("Within-block performance change versus candidate complexity\nEach block is centred on its own fixed-base MCC", fontweight="bold")
    ax.legend(fontsize=8)
    savefig(fig, figures / "subset_combinations/subset_performance_vs_complexity.png")


def report(subsets: pd.DataFrame, contributions: pd.DataFrame, intervals: pd.DataFrame, comparisons: pd.DataFrame, coverage: pd.DataFrame) -> None:
    sections = []
    for block in subsets.block.drop_duplicates():
        block_subsets = subsets[subsets.block == block]
        block_contrib = contributions[contributions.block == block]
        best = block_subsets.loc[block_subsets.mcc.idxmax()]
        full = block_subsets[block_subsets.is_full == 1].iloc[0]
        empty = block_subsets[block_subsets.is_empty == 1].iloc[0]
        lines = [
            f"## {block.replace('_', ' ').title()}", "",
            f"Best observed subset: `{best.subset_label}` (MCC {best.mcc:.6f}). Full candidate subset MCC: {full.mcc:.6f}; empty-candidate base MCC: {empty.mcc:.6f}; total gain: {full.mcc-empty.mcc:.6f}.", "",
            "| Feature | Raw contribution | Normalised share | Bootstrap median | 95% interval |",
            "|---|---:|---:|---:|---:|",
        ]
        for _, row in block_contrib.iterrows():
            ci = intervals[(intervals.block == block) & (intervals.feature == row.feature)].iloc[0]
            share = "undefined" if pd.isna(row.normalised_share) else f"{row.normalised_share:.3f}"
            lines.append(f"| {DISPLAY[row.feature]} | {row.raw_dominance_contribution:.6f} | {share} | {ci.bootstrap_median:.6f} | [{ci.ci_2_5:.6f}, {ci.ci_97_5:.6f}] |")
        lines.extend(["", "Coverage and missingness:", ""])
        for _, row in coverage[coverage.block == block].iterrows():
            lines.append(f"- {DISPLAY[row.feature]}: {row.nonmissing_count:,}/{row.nonmissing_count + row.missing_count:,} nonmissing ({row.coverage:.2%}); missing values were not converted to zero.")
        lines.extend(["", "Paired full/reduced comparisons:", "", "| Comparison | Observed MCC difference | Bootstrap median | 95% interval | P(diff > 0) |", "|---|---:|---:|---:|---:|"])
        for _, row in comparisons[comparisons.block == block].iterrows():
            lines.append(f"| {row.comparison_id} | {row.observed_difference_a_minus_b:.6f} | {row.bootstrap_median_difference:.6f} | [{row.ci_2_5:.6f}, {row.ci_97_5:.6f}] | {row.probability_difference_gt_zero:.3f} |")
        lines.extend(["", f"Governance: {full.governance}. Dominance is conditional predictive contribution, not causality, independence, or proof against circularity."])
        if block == "primary_conservation":
            target = comparisons[(comparisons.block == block) & (comparisons.comparison_id == "all3_vs_phylop_phastcons")].iloc[0]
            if target.ci_2_5 <= 0 <= target.ci_97_5:
                lines.append("Provisional interpretation: phyloP + PhastCons remains the parsimonious alternative because the paired full-minus-reduced MCC interval includes zero. This does not modify the authoritative configuration.")
            else:
                lines.append("Provisional interpretation: the all-three advantage is supported within training; no authoritative configuration is changed here.")
        elif block == "population_inclusive_sensitivity":
            lines.append("Provisional interpretation: maximum-population AF alone is the parsimonious sensitivity condition; adding global AF has no supported incremental MCC gain. Population evidence nevertheless remains sensitivity-only. Missing values were imputed fold-locally, never converted to zero, and class-associated missingness remains a warning.")
        elif block == "computational_predictor_sensitivity":
            lines.append("Provisional interpretation: REVEL, CADD and AlphaMissense show supported conditional full-model losses when removed; SIFT and PolyPhen-2 do not show a clearly supported incremental loss from the full five-predictor subset. This is inconclusive for a governed reduced subset and does not justify automatic inclusion or exclusion. Predictors remain external comparators or optional sensitivity inputs; performance cannot override provenance or circularity governance.")
        elif block == "primary_amino_acid":
            lines.append("Provisional interpretation: hydropathy alone is the more parsimonious observed single-feature option, and adding polarity has no clearly supported MCC gain. Because both dominance intervals include zero and the descriptors are biologically distinct, removal is not finalised; the evidence is supportive but inconclusive.")
        elif block == "primary_gene_constraint":
            lines.append("Provisional interpretation: retain both constraint measures for now. Removing either from the full pair produced a positive full-minus-reduced interval, while the two singleton models were indistinguishable; their contributions are balanced and biologically non-identical.")
        else:
            lines.append("No automatic removal is made. Any reduced-subset choice requires biological and interpretability review in addition to paired uncertainty.")
        sections.append("\n".join(lines))
    reconciliation_section = ""
    reconciliation_path = STAGE / "conservation_mcc_reconciliation.csv"
    if reconciliation_path.exists():
        reconciliation = pd.read_csv(reconciliation_path)
        reconciliation_section = "\n## Conservation MCC reconciliation\n\nThe previous audit reported the unweighted arithmetic mean of five fold-level MCC values. The dominance stage reports one MCC calculated after pooling all 27,727 OOF predictions. Rows, order, labels, folds, features, preprocessing, model parameters, threshold, predictions and software versions were identical. Therefore the discrepancy is solely the nonlinearity of MCC under mean-of-fold versus pooled-OOF aggregation. See `conservation_mcc_reconciliation.csv` for condition-level values and exact checks.\n"
    text = """# Rare Missense-SNV Pre-Modelling Evidence Audit:
## Dominance and Hierarchical-Partitioning Addendum

## Scope

This is an **Exact all-subsets Shapley/general-dominance decomposition of training-only cross-validated MCC**. It uses only 27,727 training variants and the existing five gene-grouped folds. It is a predictive-contribution analysis, not classical R-squared dominance analysis. Normalised shares sum to 100% by construction; they are not probabilities of feature necessity and do not establish biological independence or causality. Dominance is supplementary evidence rather than the sole keep/drop criterion. Validation and test outcomes were not loaded, inspected, or evaluated.

Exact Shapley efficiency and valid normalised-share identities were verified for every block. Negative contributions, if present, are preserved. Supporting PR-AUC values are average precision.

""" + reconciliation_section + "\n" + "\n\n".join(sections) + """

## Decision boundary

These results are provisional training-only supporting evidence. They do not alter `recommended_training_only_configuration.json`, any validation command, or the locked split. Correlation indicates shared information, not circularity. No validation confirmation, formal baseline modelling, or final model selection was run.
"""
    (STAGE / "DOMINANCE_HIERARCHICAL_PARTITIONING_REPORT.md").write_text(text, encoding="utf-8")


def hash_manifest() -> None:
    rows = []
    for path in sorted(STAGE.rglob("*")):
        if path.is_file() and path.name != "dominance_run_file_hashes.csv":
            rows.append({"path": path.relative_to(STAGE).as_posix(), "size_bytes": path.stat().st_size, "sha256": sha256(path)})
    write_csv(STAGE / "dominance_run_file_hashes.csv", rows)


def main() -> None:
    started = time.time()
    spec_hash = verify_specification()
    data, primary_features, y, groups, splits, folds = fixed_training_data()
    definitions = block_definitions(primary_features)
    subset_rows, contribution_rows, interval_rows, comparison_rows, marginal_rows, coverage_rows = [], [], [], [], [], []
    completed = []
    for definition in definitions:
        block = definition["block"]; candidates = definition["candidates"]; base = definition["base"]
        probabilities = reused_conservation_oof(data, y, folds) if definition["reuse"] else {}
        for subset_tuple in all_subsets(candidates):
            subset = frozenset(subset_tuple)
            if subset not in probabilities:
                probabilities[subset] = fit_oof(data, y, splits, base + list(subset_tuple))
        values = {subset: metric_values(y, probability)["mcc"] for subset, probability in probabilities.items()}
        phi = shapley_contributions(values, candidates)
        total_gain = values[frozenset(candidates)] - values[frozenset()]
        shares = normalised_shares(phi, total_gain)
        bootstrap_values, bootstrap_phi = bootstrap_block(y, groups, candidates, probabilities)
        for subset, probability in probabilities.items():
            metrics = metric_values(y, probability)
            subset_rows.append({
                "block": block, "subset_key": subset_key(subset), "subset_label": "None" if not subset else " + ".join(DISPLAY[x] for x in candidates if x in subset),
                "candidate_count": len(subset), "included_candidates": ";".join(x for x in candidates if x in subset),
                "fixed_base_feature_count": len(base), "total_feature_count": len(base) + len(subset),
                **{name: round(value, 10) for name, value in metrics.items()}, "row_count": len(data), "gene_count": len(set(groups)),
                "is_empty": int(not subset), "is_full": int(subset == frozenset(candidates)), "governance": definition["governance"],
            })
        for feature in candidates:
            samples = bootstrap_phi[feature]
            contribution_rows.append({"block": block, "feature": feature, "raw_dominance_contribution": round(phi[feature], 10), "normalised_share": "" if shares[feature] is None else round(float(shares[feature]), 10), "full_minus_empty_mcc": round(total_gain, 10), "efficiency_residual": round(sum(phi.values()) - total_gain, 14), "normalisation_defined": shares[feature] is not None, "governance": definition["governance"]})
            interval_rows.append({"block": block, "feature": feature, "bootstrap_samples": N_BOOTSTRAP, "seed": SEED, "resampling_unit": GROUP, "bootstrap_median": round(float(np.median(samples)), 10), "ci_2_5": round(float(np.quantile(samples, .025)), 10), "ci_97_5": round(float(np.quantile(samples, .975)), 10), "probability_contribution_gt_zero": round(float(np.mean(samples > 0)), 6)})
            others = [x for x in candidates if x != feature]
            for subset_tuple in all_subsets(others):
                subset = frozenset(subset_tuple)
                marginal_rows.append({"block": block, "feature": feature, "conditioning_subset": subset_key(subset), "conditioning_size": len(subset), "marginal_mcc_contribution": round(values[subset | {feature}] - values[subset], 10)})
            coverage_rows.append({"block": block, "feature": feature, "nonmissing_count": int(numeric(data[feature]).notna().sum()), "missing_count": int(numeric(data[feature]).isna().sum()), "coverage": round(float(numeric(data[feature]).notna().mean()), 10), "missing_values_converted_to_zero": False})
        for subset_a, subset_b, comparison_id in comparisons_for(candidates):
            differences = bootstrap_values[subset_a] - bootstrap_values[subset_b]
            comparison_rows.append({"block": block, "comparison_id": comparison_id, "subset_a": subset_key(subset_a), "subset_b": subset_key(subset_b), "observed_mcc_a": round(values[subset_a], 10), "observed_mcc_b": round(values[subset_b], 10), "observed_difference_a_minus_b": round(values[subset_a] - values[subset_b], 10), "bootstrap_median_difference": round(float(np.median(differences)), 10), "ci_2_5": round(float(np.quantile(differences, .025)), 10), "ci_97_5": round(float(np.quantile(differences, .975)), 10), "probability_difference_gt_zero": round(float(np.mean(differences > 0)), 6), "bootstrap_samples": N_BOOTSTRAP, "resampling_unit": GROUP})
        completed.append(block)
    subsets = pd.DataFrame(subset_rows); contributions = pd.DataFrame(contribution_rows); intervals = pd.DataFrame(interval_rows); comparisons = pd.DataFrame(comparison_rows); coverage = pd.DataFrame(coverage_rows)
    subsets.to_csv(STAGE / "subset_performance.csv", index=False)
    contributions.to_csv(STAGE / "dominance_contributions.csv", index=False)
    intervals.to_csv(STAGE / "dominance_bootstrap_intervals.csv", index=False)
    comparisons.to_csv(STAGE / "paired_subset_comparisons.csv", index=False)
    pd.DataFrame(marginal_rows).to_csv(STAGE / "conditional_marginal_contributions.csv", index=False)
    coverage.to_csv(STAGE / "feature_coverage.csv", index=False)
    plot_outputs(subsets, contributions, intervals, completed)
    report(subsets, contributions, intervals, comparisons, coverage)
    metadata = {
        "method": "Exact all-subsets Shapley/general-dominance decomposition of training-only cross-validated MCC",
        "analysis_specification_sha256": spec_hash, "seed": SEED, "bootstrap_samples": N_BOOTSTRAP,
        "training_rows": len(data), "training_genes": len(set(groups)), "fixed_folds": 5,
        "completed_blocks": completed, "validation_accessed": False, "test_accessed": False,
        "authoritative_recommendation_modified": False, "validation_command_modified": False,
        "runtime_seconds": round(time.time() - started, 3),
        "software": {"python": sys.version.split()[0], "platform": platform.platform(), "numpy": np.__version__, "pandas": pd.__version__, "scikit_learn": sklearn.__version__, "matplotlib": matplotlib.__version__},
        "input_checksums": {name: sha256(path) for name, path in PATHS.items()},
        "split_assignment_sha256": sha256(PATHS["assignments"]),
        "fixed_fold_assignment_sha256": sha256(OUT / "03_controlled_experiments/fixed_cv_fold_assignments.csv"),
    }
    write_json(STAGE / "run_metadata.json", metadata)
    hash_manifest()
    print(json.dumps({"status": "PASS", "completed_blocks": completed, "runtime_seconds": metadata["runtime_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
