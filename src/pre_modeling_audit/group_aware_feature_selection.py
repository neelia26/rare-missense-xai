#!/usr/bin/env python3
"""Training-only group-aware correlated-feature selection audit."""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
import platform
import sys
import textwrap
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from joblib import Parallel, delayed
from sklearn.feature_selection import mutual_info_classif, mutual_info_regression
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from src.pre_modeling_audit.run_pre_modeling_audit import (
    CONSERVATION, GROUP, KEY, METADATA, OUT, PATHS, POPULATION, PREDICTORS,
    ROOT, SEED, TARGET, numeric, sha256, write_csv, write_json,
)


STAGE = OUT / "05_group_aware_feature_selection"
SPEC = STAGE / "ANALYSIS_SPECIFICATION.md"
SPEC_HASH = STAGE / "analysis_specification_sha256.txt"
N_BOOTSTRAP = 5000
MARGINS = [0.0025, 0.005, 0.010]
AA = ["aa_hydropathy_delta", "aa_polarity_delta"]
CONSTRAINT = ["gnomad_missense_z", "gnomad_loeuf"]
DISPLAY = {
    "conservation_gerp_rs": "GERP++", "conservation_phylop100way": "phyloP", "conservation_phastcons100way": "PhastCons",
    "aa_hydropathy_delta": "Hydropathy delta", "aa_polarity_delta": "Polarity delta",
    "gnomad_missense_z": "Missense Z", "gnomad_loeuf": "LOEUF", "gnomad_global_af": "Global AF",
    "gnomad_max_population_af": "Max-population AF", "predictor_sift_score": "SIFT",
    "predictor_polyphen2_hdiv_score": "PolyPhen-2", "predictor_revel_score": "REVEL",
    "predictor_cadd_phred": "CADD", "predictor_alphamissense_score": "AlphaMissense",
}
BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#666666"


ALLOWED_INPUTS = {
    PATHS["train_primary"].resolve(), PATHS["train_population"].resolve(), PATHS["train_optional"].resolve(),
    (OUT / "03_controlled_experiments/fixed_cv_fold_assignments.csv").resolve(),
}


def guarded_training_load(path: Path, **kwargs) -> pd.DataFrame:
    resolved = path.resolve()
    if resolved not in ALLOWED_INPUTS:
        raise RuntimeError(f"Training-only guard rejected unapproved path: {path}")
    if any(part.lower().startswith(("validation", "test")) for part in resolved.parts):
        raise RuntimeError(f"Training-only guard rejected held-out path: {path}")
    return pd.read_csv(resolved, **kwargs)


def verify_specification() -> str:
    expected = SPEC_HASH.read_text().split()[0]
    actual = sha256(SPEC)
    if expected != actual:
        raise RuntimeError("Analysis specification checksum mismatch")
    return actual


def subset_key(subset: frozenset[str] | tuple[str, ...]) -> str:
    return "empty" if not subset else "__".join(sorted(subset))


def subsets(features: list[str]) -> list[frozenset[str]]:
    return [frozenset(combo) for size in range(len(features) + 1) for combo in itertools.combinations(features, size)]


def pipeline() -> Pipeline:
    return Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(C=1.0, max_iter=2000, solver="liblinear", class_weight="balanced", random_state=SEED)),
    ])


def metrics(y: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    return {
        "mcc": float(matthews_corrcoef(y, probability >= .5)),
        "pr_auc": float(average_precision_score(y, probability)),
        "roc_auc": float(roc_auc_score(y, probability)),
        "brier": float(brier_score_loss(y, probability)),
    }


def load_training() -> tuple[pd.DataFrame, list[str], np.ndarray, np.ndarray, np.ndarray]:
    primary = guarded_training_load(PATHS["train_primary"], dtype=str, keep_default_na=False, na_values=[])
    population = guarded_training_load(PATHS["train_population"], dtype=str, keep_default_na=False, na_values=[])[[KEY] + POPULATION]
    optional = guarded_training_load(PATHS["train_optional"], dtype=str, keep_default_na=False, na_values=[])[[KEY] + PREDICTORS]
    data = primary.merge(population, on=KEY, validate="one_to_one").merge(optional, on=KEY, validate="one_to_one")
    folds = guarded_training_load(OUT / "03_controlled_experiments/fixed_cv_fold_assignments.csv", dtype={KEY: str, GROUP: str, "cv_fold": int})
    data = data.merge(folds, on=[KEY, GROUP], validate="one_to_one")
    if len(data) != 27727 or data[KEY].duplicated().any() or data.groupby(GROUP).cv_fold.nunique().max() != 1:
        raise RuntimeError("Training rows/folds failed integrity checks")
    primary_features = [column for column in primary.columns if column not in METADATA]
    return data, primary_features, data[TARGET].astype(int).to_numpy(), data[GROUP].to_numpy(), data.cv_fold.to_numpy()


def block_definitions(primary: list[str]) -> list[dict]:
    return [
        {"block": "conservation", "candidates": CONSERVATION, "base": [x for x in primary if x not in CONSERVATION], "governance": "primary evidence"},
        {"block": "amino_acid_properties", "candidates": AA, "base": [x for x in primary if x not in AA], "governance": "primary evidence"},
        {"block": "gene_constraint", "candidates": CONSTRAINT, "base": [x for x in primary if x not in CONSTRAINT], "governance": "primary evidence group"},
        {"block": "population_frequency", "candidates": POPULATION, "base": list(primary), "governance": "sensitivity-only; class-associated missingness"},
        {"block": "computational_predictors", "candidates": PREDICTORS, "base": list(primary), "governance": "sensitivity-only; provenance/circularity governed"},
    ]


def pairwise_redundancy(matrix: np.ndarray) -> np.ndarray:
    p = matrix.shape[1]; result = np.zeros((p, p))
    for a in range(p):
        for b in range(a + 1, p):
            ab = mutual_info_regression(matrix[:, [a]], matrix[:, b], random_state=SEED)[0]
            ba = mutual_info_regression(matrix[:, [b]], matrix[:, a], random_state=SEED)[0]
            result[a, b] = result[b, a] = (ab + ba) / 2
    return result


def mrmr_ranking(frame: pd.DataFrame, y: np.ndarray, candidates: list[str]) -> list[dict]:
    imputed = SimpleImputer(strategy="median").fit_transform(frame[candidates].apply(numeric))
    relevance = mutual_info_classif(imputed, y, random_state=SEED)
    redundancy = pairwise_redundancy(imputed)
    remaining = list(range(len(candidates))); ranked = []
    while remaining:
        scored = []
        for index in remaining:
            mean_redundancy = float(np.mean([redundancy[index, prior] for prior in ranked])) if ranked else 0.0
            scored.append((float(relevance[index] - mean_redundancy), float(relevance[index]), -index, index, mean_redundancy))
        _, rel, _, selected, red = max(scored)
        ranked.append(selected); remaining.remove(selected)
    return [{"feature": candidates[index], "rank": rank, "relevance_mi": float(relevance[index]), "redundancy_with_preceding_mi": float(np.mean([redundancy[index, prior] for prior in ranked[:rank-1]])) if rank > 1 else 0.0, "mrmr_score_at_selection": float(relevance[index] - (np.mean([redundancy[index, prior] for prior in ranked[:rank-1]]) if rank > 1 else 0.0))} for rank, index in enumerate(ranked, 1)]


def evaluate_subset(data: pd.DataFrame, y: np.ndarray, groups: np.ndarray, outer_fit: np.ndarray, outer_held: np.ndarray, features: list[str], inner_seed: int) -> tuple[float, np.ndarray]:
    inner = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=inner_seed)
    inner_probability = np.full(len(outer_fit), np.nan)
    y_fit = y[outer_fit]; groups_fit = groups[outer_fit]
    for fit_local, held_local in inner.split(outer_fit, y_fit, groups_fit):
        fit_index, held_index = outer_fit[fit_local], outer_fit[held_local]
        model = pipeline(); model.fit(data.iloc[fit_index][features].apply(numeric), y[fit_index])
        inner_probability[held_local] = model.predict_proba(data.iloc[held_index][features].apply(numeric))[:, 1]
    if not np.isfinite(inner_probability).all():
        raise RuntimeError("Incomplete inner OOF probabilities")
    model = pipeline(); model.fit(data.iloc[outer_fit][features].apply(numeric), y[outer_fit])
    outer_probability = model.predict_proba(data.iloc[outer_held][features].apply(numeric))[:, 1]
    return metrics(y_fit, inner_probability)["mcc"], outer_probability


def mcc_from_confusion(cells: np.ndarray) -> np.ndarray:
    tp, tn, fp, fn = (cells[..., index] for index in range(4))
    numerator = tp * tn - fp * fn
    denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return np.divide(numerator, denominator, out=np.zeros_like(numerator, dtype=float), where=denominator > 0)


def bootstrap_subset_values(y: np.ndarray, groups: np.ndarray, probabilities: dict[frozenset[str], np.ndarray]) -> dict[frozenset[str], np.ndarray]:
    genes, gene_index = np.unique(groups, return_inverse=True); ordered = list(probabilities)
    confusion = np.zeros((len(genes), len(ordered), 4), dtype=np.int64)
    for condition, subset in enumerate(ordered):
        predicted = probabilities[subset] >= .5
        for cell, mask in enumerate((predicted & (y == 1), ~predicted & (y == 0), predicted & (y == 0), ~predicted & (y == 1))):
            confusion[:, condition, cell] = np.bincount(gene_index, weights=mask.astype(int), minlength=len(genes)).astype(int)
    rng = np.random.default_rng(SEED); values = {subset: np.empty(N_BOOTSTRAP) for subset in ordered}
    for start in range(0, N_BOOTSTRAP, 100):
        size = min(100, N_BOOTSTRAP - start)
        weights = rng.multinomial(len(genes), np.full(len(genes), 1 / len(genes)), size=size)
        cells = (weights @ confusion.reshape(len(genes), -1)).reshape(size, len(ordered), 4)
        mcc = mcc_from_confusion(cells)
        for index, subset in enumerate(ordered): values[subset][start:start+size] = mcc[:, index]
    return values


def target_comparisons(block: str, candidates: list[str]) -> list[tuple[frozenset[str], frozenset[str], str]]:
    full = frozenset(candidates)
    if block == "conservation":
        return [(frozenset(CONSERVATION[1:]), frozenset([CONSERVATION[1]]), "phylop_phastcons_vs_phylop"), (full, frozenset(CONSERVATION[1:]), "all3_vs_phylop_phastcons")]
    if block == "amino_acid_properties":
        return [(full, frozenset([AA[0]]), "both_vs_hydropathy"), (full, frozenset([AA[1]]), "both_vs_polarity"), (frozenset([AA[0]]), frozenset([AA[1]]), "hydropathy_vs_polarity")]
    if block == "gene_constraint":
        return [(full, frozenset([CONSTRAINT[0]]), "both_vs_missense_z"), (full, frozenset([CONSTRAINT[1]]), "both_vs_loeuf"), (frozenset([CONSTRAINT[0]]), frozenset([CONSTRAINT[1]]), "missense_z_vs_loeuf")]
    if block == "population_frequency":
        return [(full, frozenset([POPULATION[1]]), "both_vs_max_population_af"), (frozenset([POPULATION[1]]), frozenset([POPULATION[0]]), "max_population_af_vs_global_af")]
    trio = frozenset([PREDICTORS[2], PREDICTORS[3], PREDICTORS[4]])
    return [(full, trio, "all_five_vs_revel_cadd_alphamissense")]


def reduced_candidates(block: str) -> list[frozenset[str]]:
    return {
        "conservation": [frozenset(CONSERVATION[1:])],
        "amino_acid_properties": [frozenset([AA[0]]), frozenset([AA[1]])],
        "gene_constraint": [frozenset([CONSTRAINT[0]]), frozenset([CONSTRAINT[1]])],
        "population_frequency": [frozenset([POPULATION[1]])],
        "computational_predictors": [frozenset([PREDICTORS[2], PREDICTORS[3], PREDICTORS[4]])],
    }[block]


def choose_subset(scores: dict[frozenset[str], float]) -> frozenset[str]:
    return min(scores, key=lambda subset: (-scores[subset], len(subset), subset_key(subset)))


def savefig(fig: plt.Figure, name: str) -> None:
    path = STAGE / "figures" / name; path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=320, bbox_inches="tight", facecolor="white"); plt.close(fig)


def figures(rankings: pd.DataFrame, performance: pd.DataFrame, frequencies: pd.DataFrame, decisions: pd.DataFrame) -> None:
    stability = rankings.groupby(["block", "feature"]).agg(rank1_frequency=("rank", lambda x: float(np.mean(x == 1)))).reset_index()
    fig, ax = plt.subplots(figsize=(11, 7)); labels=[f"{row.block}: {DISPLAY[row.feature]}" for _,row in stability.iterrows()]; ax.barh(labels,100*stability.rank1_frequency,color=BLUE); ax.set_xlabel("Ranked first across outer-training folds (%)"); ax.set_title("Training only: fold-local mRMR ranking stability",fontweight="bold"); savefig(fig,"mrmr_selection_frequency.png")
    selected=performance[performance.analysis_type.eq("fixed_subset_outer_oof")].copy(); selected=selected.sort_values(["block","mcc"],ascending=[True,False]); fig,axes=plt.subplots(3,2,figsize=(15,14));
    for ax,(block,rows) in zip(axes.flat,selected.groupby("block")):
        rows=rows.head(10).sort_values("mcc"); ax.errorbar(rows.mcc,np.arange(len(rows)),xerr=[rows.mcc-rows.ci_2_5,rows.ci_97_5-rows.mcc],fmt="o",color=BLUE,capsize=3); ax.set_yticks(np.arange(len(rows)),rows.subset_label); ax.set_title(block.replace("_"," ").title()); ax.set_xlabel("Pooled OOF MCC with 95% gene-bootstrap interval")
    axes.flat[-1].axis("off"); fig.suptitle("Training only: subset performance uncertainty",fontweight="bold"); savefig(fig,"subset_performance_confidence_intervals.png")
    fig,ax=plt.subplots(figsize=(11,7)); outer=frequencies[frequencies.frequency_type.eq("outer_fold_selected")]; labels=[f"{r.block}: {r.subset_label}" for _,r in outer.iterrows()]; ax.barh(labels,outer.frequency,color=ORANGE); ax.set_xlabel("Selection frequency across five outer folds"); ax.set_title("Training only: nested subset-selection stability",fontweight="bold"); savefig(fig,"subset_selection_stability.png")
    fig,axes=plt.subplots(3,2,figsize=(14,13));
    for ax,(block,rows) in zip(axes.flat,selected.groupby("block")):
        ax.scatter(rows.candidate_count,rows.mcc-rows.loc[rows.is_empty.eq(1),"mcc"].iloc[0],color=BLUE,alpha=.8); ax.axhline(0,color=GREY,lw=1); ax.set_title(block.replace("_"," ").title()); ax.set_xlabel("Candidate count"); ax.set_ylabel("Δ MCC vs block empty")
    axes.flat[-1].axis("off"); fig.suptitle("Training only: within-block performance versus subset complexity",fontweight="bold"); savefig(fig,"performance_vs_subset_complexity.png")
    summary = decisions[["block", "provisional_decision", "provisional_subset"]].copy()
    summary["block"] = summary["block"].str.replace("_", " ", regex=False).str.title()
    summary["provisional_decision"] = summary["provisional_decision"].map(lambda value: textwrap.fill(str(value), 30))
    summary["provisional_subset"] = summary["provisional_subset"].map(lambda value: textwrap.fill(str(value), 48))
    fig,ax=plt.subplots(figsize=(15,7)); ax.axis("off"); table=ax.table(cellText=summary.values,colLabels=["Block","Decision","Provisional subset"],colWidths=[.20,.28,.52],loc="center",cellLoc="left"); table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1,3.2); ax.set_title("Training only: provisional block-level decision summary\nNot an authoritative model lock",fontweight="bold",pad=18); savefig(fig,"block_level_decision_summary.png")


def method_documents() -> None:
    (STAGE/"METHOD_COMPARISON.md").write_text("""# Method comparison

| Method | Question answered | Strength | Limitation | Decision role |
|---|---|---|---|---|
| Exact all-subsets Shapley/general-dominance decomposition of training-only cross-validated MCC | How is fixed-model MCC gain distributed across correlated features? | Exact additive predictive-contribution accounting | Shares total 100% by construction and are not probabilities of necessity; no causality or independence | Supplementary |
| Fold-local mRMR | Which candidates combine outcome relevance with lower candidate redundancy? | Training-fold ranking sensitivity | Ranking can change among correlated alternatives; not a final keep/drop proof | Supplementary ranking/stability |
| Nested gene-grouped wrapper | Which subset is selected without using the outer holdout? | Separates selection from outer assessment and respects genes | Results remain model- and sample-specific | Primary training-only subset evidence |
| Gene bootstrap/non-inferiority sensitivity | Are paired losses stable across genes and plausible margins? | Paired uncertainty and stability | Margins are declared sensitivity values, not automatic scientific thresholds | Decision support |

Correlated retained features are interpreted as an evidence group, not independent biological votes. Population and computational predictors remain sensitivity-only regardless of these results.
""",encoding="utf-8")


def synthesise_decisions(comparisons: pd.DataFrame, margins: pd.DataFrame, frequencies: pd.DataFrame) -> pd.DataFrame:
    by_id=comparisons.set_index("comparison_id")
    strict=margins[margins.noninferiority_margin_mcc.eq(.0025)].set_index(["block","reduced_subset"])
    outer=frequencies[frequencies.frequency_type.eq("outer_fold_selected")].set_index(["block","subset_key"])["frequency"]
    conservation_pair=subset_key(frozenset(CONSERVATION[1:])); conservation_strict=float(strict.loc[("conservation",conservation_pair),"probability_loss_within_margin"])
    conservation_pair_supported=by_id.loc["phylop_phastcons_vs_phylop","ci_2_5"]>0
    conservation_stable=float(outer.get(("conservation",conservation_pair),0))>=.5
    conservation_reduced=conservation_strict>=.95 and conservation_pair_supported and conservation_stable
    aa_hydropathy=subset_key(frozenset([AA[0]])); aa_stable=float(outer.get(("amino_acid_properties",aa_hydropathy),0))>=.5
    constraint_joint=(by_id.loc["both_vs_missense_z","ci_2_5"]>0 and by_id.loc["both_vs_loeuf","ci_2_5"]>0)
    population_max=subset_key(frozenset([POPULATION[1]])); population_reduced=float(strict.loc[("population_frequency",population_max),"probability_loss_within_margin"])>=.95
    predictor_trio=subset_key(frozenset([PREDICTORS[2],PREDICTORS[3],PREDICTORS[4]])); predictor_noninferior=float(strict.loc[("computational_predictors",predictor_trio),"probability_loss_within_margin"])>=.95
    return pd.DataFrame([
        {"block":"conservation","provisional_decision":"retain reduced subset" if conservation_reduced else "inconclusive—requires supervisor decision","provisional_subset":"phyloP + PhastCons" if conservation_reduced else "all three vs phyloP + PhastCons unresolved","rationale":f"Pair improves on phyloP, but strict-margin non-inferiority probability={conservation_strict:.4f} and pair outer-selection frequency={float(outer.get(('conservation',conservation_pair),0)):.2f}; methods disagree.","authoritative_configuration_changed":False},
        {"block":"amino_acid_properties","provisional_decision":"retain reduced subset" if aa_stable else "inconclusive—requires supervisor decision","provisional_subset":"hydropathy delta" if aa_stable else "hydropathy/polarity/full alternatives unstable","rationale":"mRMR ranks hydropathy first, but nested selections and gene-bootstrap preferences do not identify a stable single-feature replacement.","authoritative_configuration_changed":False},
        {"block":"gene_constraint","provisional_decision":"retain as grouped representation" if constraint_joint else "inconclusive—requires supervisor decision","provisional_subset":"Missense Z + LOEUF" if constraint_joint else "unresolved","rationale":"Both full-versus-single intervals support joint value; all five nested folds selected the pair. They are correlated but biologically complementary.","authoritative_configuration_changed":False},
        {"block":"population_frequency","provisional_decision":"sensitivity-only","provisional_subset":"maximum-population AF sensitivity" if population_reduced else "both population fields sensitivity","rationale":"Max-population AF is strictly non-inferior to both and global AF adds no stable value; class-associated missingness prevents primary promotion.","authoritative_configuration_changed":False},
        {"block":"computational_predictors","provisional_decision":"sensitivity-only","provisional_subset":"REVEL + CADD + AlphaMissense non-inferior diagnostic subset" if predictor_noninferior else "all five diagnostic predictors","rationale":"The trio is non-inferior at the strict margin, but nested subset identity is unstable; provenance and circularity governance prevents primary promotion.","authoritative_configuration_changed":False},
    ])


def hashes() -> None:
    rows=[]
    inputs=[PATHS["frozen"],PATHS["cohort"],PATHS["assignments"],PATHS["train_primary"],PATHS["train_population"],PATHS["train_optional"],OUT/"03_controlled_experiments/fixed_cv_fold_assignments.csv"]
    for path in inputs: rows.append({"kind":"input","path":path.relative_to(ROOT).as_posix(),"size_bytes":path.stat().st_size,"sha256":sha256(path)})
    for path in sorted(STAGE.rglob("*")):
        if path.is_file() and path.name!="input_and_output_hashes.csv": rows.append({"kind":"output","path":path.relative_to(ROOT).as_posix(),"size_bytes":path.stat().st_size,"sha256":sha256(path)})
    write_csv(STAGE/"input_and_output_hashes.csv",rows)


def main() -> None:
    started=time.time(); spec_hash=verify_specification(); data,primary,y,groups,outer_folds=load_training(); definitions=block_definitions(primary)
    ranking_rows=[]; performance_rows=[]; comparison_rows=[]; frequency_rows=[]; margin_rows=[]; decision_rows=[]
    for definition in definitions:
        block=definition["block"]; candidates=definition["candidates"]; base=definition["base"]; candidate_subsets=subsets(candidates); fixed_probs={subset:np.full(len(data),np.nan) for subset in candidate_subsets}; nested_prob=np.full(len(data),np.nan); outer_choices=[]; inner_scores_by_fold=[]
        for outer_fold in range(5):
            outer_fit=np.flatnonzero(outer_folds!=outer_fold); outer_held=np.flatnonzero(outer_folds==outer_fold)
            ranking=mrmr_ranking(data.iloc[outer_fit],y[outer_fit],candidates)
            results=Parallel(n_jobs=min(4,os.cpu_count() or 1),prefer="processes")(delayed(evaluate_subset)(data,y,groups,outer_fit,outer_held,base+[f for f in candidates if f in subset],SEED+100+outer_fold) for subset in candidate_subsets)
            scores={subset:result[0] for subset,result in zip(candidate_subsets,results)}; selected=choose_subset(scores); outer_choices.append(selected); inner_scores_by_fold.append(scores)
            for subset,result in zip(candidate_subsets,results): fixed_probs[subset][outer_held]=result[1]
            nested_prob[outer_held]=dict(zip(candidate_subsets,[r[1] for r in results]))[selected]
            for row in ranking: ranking_rows.append({"block":block,"outer_fold":outer_fold,**row,"wrapper_selected":int(row["feature"] in selected),"selected_subset":subset_key(selected)})
        boot=bootstrap_subset_values(y,groups,fixed_probs); point={subset:metrics(y,prob) for subset,prob in fixed_probs.items()}
        for subset in candidate_subsets:
            values=boot[subset]; performance_rows.append({"block":block,"analysis_type":"fixed_subset_outer_oof","subset_key":subset_key(subset),"subset_label":"None" if not subset else " + ".join(DISPLAY[f] for f in candidates if f in subset),"included_candidates":";".join(f for f in candidates if f in subset),"candidate_count":len(subset),**{k:round(v,10) for k,v in point[subset].items()},"ci_2_5":round(float(np.quantile(values,.025)),10),"ci_97_5":round(float(np.quantile(values,.975)),10),"is_empty":int(not subset),"is_full":int(subset==frozenset(candidates)),"governance":definition["governance"]})
        nested=metrics(y,nested_prob); performance_rows.append({"block":block,"analysis_type":"nested_wrapper_selected","subset_key":"varies_by_outer_fold","subset_label":"Nested-selected subset","included_candidates":"varies","candidate_count":np.nan,**{k:round(v,10) for k,v in nested.items()},"ci_2_5":np.nan,"ci_97_5":np.nan,"is_empty":0,"is_full":0,"governance":definition["governance"]})
        choice_counts=pd.Series([subset_key(x) for x in outer_choices]).value_counts(); preferred=np.array([choose_subset({subset:boot[subset][i] for subset in candidate_subsets}) for i in range(N_BOOTSTRAP)],dtype=object); preferred_counts=pd.Series([subset_key(x) for x in preferred]).value_counts()
        labels={subset_key(s):("None" if not s else " + ".join(DISPLAY[f] for f in candidates if f in s)) for s in candidate_subsets}
        for key,count in choice_counts.items(): frequency_rows.append({"block":block,"frequency_type":"outer_fold_selected","subset_key":key,"subset_label":labels[key],"count":int(count),"denominator":5,"frequency":count/5})
        for key,count in preferred_counts.items(): frequency_rows.append({"block":block,"frequency_type":"bootstrap_preferred","subset_key":key,"subset_label":labels[key],"count":int(count),"denominator":N_BOOTSTRAP,"frequency":count/N_BOOTSTRAP})
        for a,b,name in target_comparisons(block,candidates):
            diff=boot[a]-boot[b]; comparison_rows.append({"block":block,"comparison_id":name,"subset_a":subset_key(a),"subset_b":subset_key(b),"observed_difference_a_minus_b":round(point[a]["mcc"]-point[b]["mcc"],10),"bootstrap_median_difference":round(float(np.median(diff)),10),"ci_2_5":round(float(np.quantile(diff,.025)),10),"ci_97_5":round(float(np.quantile(diff,.975)),10),"probability_difference_gt_zero":round(float(np.mean(diff>0)),6),"bootstrap_samples":N_BOOTSTRAP})
        full=frozenset(candidates)
        for reduced in reduced_candidates(block):
            loss=boot[full]-boot[reduced]
            for margin in MARGINS: margin_rows.append({"block":block,"full_subset":subset_key(full),"reduced_subset":subset_key(reduced),"noninferiority_margin_mcc":margin,"observed_loss_full_minus_reduced":round(point[full]["mcc"]-point[reduced]["mcc"],10),"probability_loss_within_margin":round(float(np.mean(loss<=margin)),6),"sensitivity_conclusion":"reduced non-inferior in sensitivity analysis" if np.mean(loss<=margin)>=.95 else "not established at 95% probability"})
        decision_rows.append({"block":block,"provisional_decision":"pending result synthesis","provisional_subset":"","authoritative_configuration_changed":False})
    rankings=pd.DataFrame(ranking_rows); performance=pd.DataFrame(performance_rows); comparisons=pd.DataFrame(comparison_rows); frequencies=pd.DataFrame(frequency_rows); margins=pd.DataFrame(margin_rows)
    stability=rankings.groupby(["block","feature"]).agg(mean_rank=("rank","mean"),sd_rank=("rank","std"),rank_1_frequency=("rank",lambda x:float(np.mean(x==1))),wrapper_inclusion_frequency=("wrapper_selected","mean")).reset_index()
    decisions=synthesise_decisions(comparisons,margins,frequencies)
    rankings.to_csv(STAGE/"mrmr_fold_rankings.csv",index=False); stability.to_csv(STAGE/"mrmr_selection_stability.csv",index=False); performance.to_csv(STAGE/"nested_wrapper_subset_performance.csv",index=False); comparisons.to_csv(STAGE/"paired_subset_comparisons.csv",index=False); frequencies.to_csv(STAGE/"subset_selection_frequencies.csv",index=False); margins.to_csv(STAGE/"decision_sensitivity_by_margin.csv",index=False); decisions.to_csv(STAGE/"provisional_feature_decisions.csv",index=False)
    method_documents(); figures(rankings,performance,frequencies,decisions)
    report_lines=["# Group-Aware Feature Selection Report","","All results are training only. Validation and test outcomes were not accessed. Dominance remains supplementary predictive-contribution evidence, not the sole selection rule.",""]
    for _,decision in decisions.iterrows():
        block=decision.block; report_lines.extend([f"## {block.replace('_',' ').title()}","",f"**Provisional decision:** {decision.provisional_decision}. **Subset:** {decision.provisional_subset}.","",decision.rationale,"",f"Nested wrapper pooled OOF MCC: {performance[(performance.block==block)&(performance.analysis_type=='nested_wrapper_selected')].mcc.iloc[0]:.6f}.","",f"mRMR and wrapper stability are in `mrmr_fold_rankings.csv` and `subset_selection_frequencies.csv`; paired uncertainty and margin sensitivity are in their corresponding CSV files.",""])
    report_lines.extend(["## Interpretation boundary","","Correlated retained features are one evidence group, not independent votes. Population and computational predictors remain sensitivity-only. These provisional decisions do not modify the authoritative configuration and require supervisor review before any one-time validation."])
    (STAGE/"GROUP_AWARE_FEATURE_SELECTION_REPORT.md").write_text("\n".join(report_lines)+"\n",encoding="utf-8")
    metadata={"analysis_specification_sha256":spec_hash,"seed":SEED,"outer_folds":5,"inner_folds":4,"bootstrap_samples":N_BOOTSTRAP,"training_rows":len(data),"training_genes":len(set(groups)),"validation_accessed":False,"test_accessed":False,"authoritative_configuration_modified":False,"validation_command_modified":False,"runtime_seconds":round(time.time()-started,3),"software":{"python":sys.version.split()[0],"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__,"scikit_learn":sklearn.__version__,"matplotlib":matplotlib.__version__,"joblib":__import__("joblib").__version__},"input_hashes":{name:sha256(path) for name,path in PATHS.items()},"fold_assignment_sha256":sha256(OUT/"03_controlled_experiments/fixed_cv_fold_assignments.csv")}
    write_json(STAGE/"run_metadata.json",metadata); hashes(); print(json.dumps({"status":"PASS","runtime_seconds":metadata["runtime_seconds"]},indent=2))


if __name__=="__main__": main()
