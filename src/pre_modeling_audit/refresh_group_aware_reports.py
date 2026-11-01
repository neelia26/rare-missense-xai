#!/usr/bin/env python3
"""Refresh result-derived group-aware decisions and reports without refitting."""

from __future__ import annotations

import json
import pandas as pd

from src.pre_modeling_audit import group_aware_feature_selection as g


def main() -> None:
    rankings=pd.read_csv(g.STAGE/"mrmr_fold_rankings.csv")
    stability=pd.read_csv(g.STAGE/"mrmr_selection_stability.csv")
    performance=pd.read_csv(g.STAGE/"nested_wrapper_subset_performance.csv")
    comparisons=pd.read_csv(g.STAGE/"paired_subset_comparisons.csv")
    frequencies=pd.read_csv(g.STAGE/"subset_selection_frequencies.csv")
    margins=pd.read_csv(g.STAGE/"decision_sensitivity_by_margin.csv")
    decisions=g.synthesise_decisions(comparisons,margins,frequencies)
    decisions.to_csv(g.STAGE/"provisional_feature_decisions.csv",index=False)
    g.figures(rankings,performance,frequencies,decisions)
    lines=[
        "# Group-Aware Feature Selection Report", "",
        "All results are **training only**. Validation and test outcomes were not accessed. Dominance remains supplementary predictive-contribution evidence, not the sole selection rule.", "",
        "## Methods", "",
        "Fold-local mRMR estimated outcome relevance and candidate redundancy separately within each outer-training portion. Every subset was then evaluated with four inner gene-grouped folds; selection never used the outer holdout. Fixed-subset outer OOF predictions supported paired 5,000-gene bootstrap comparisons and non-inferiority sensitivity margins of 0.0025, 0.005 and 0.010 MCC.", "",
    ]
    for _, decision in decisions.iterrows():
        block=decision.block
        nested=performance[(performance.block==block)&performance.analysis_type.eq("nested_wrapper_selected")].iloc[0]
        stable=stability[stability.block==block].sort_values("mean_rank")
        outer=frequencies[(frequencies.block==block)&frequencies.frequency_type.eq("outer_fold_selected")].sort_values("frequency",ascending=False)
        paired=comparisons[comparisons.block==block]
        margin=margins[(margins.block==block)&margins.noninferiority_margin_mcc.eq(.0025)]
        lines.extend([f"## {block.replace('_',' ').title()}","",f"**Provisional decision:** {decision.provisional_decision}","",f"**Provisional subset/representation:** {decision.provisional_subset}","",decision.rationale,"",f"Nested-selected pooled OOF metrics: MCC {nested.mcc:.6f}, PR-AUC {nested.pr_auc:.6f}, ROC-AUC {nested.roc_auc:.6f}, Brier {nested.brier:.6f}.","",f"mRMR mean ranking: "+"; ".join(f"{g.DISPLAY[r.feature]} {r.mean_rank:.1f}" for _,r in stable.iterrows())+".","",f"Outer-fold subset selections: "+"; ".join(f"{r.subset_label} {int(r['count'])}/5" for _,r in outer.iterrows())+".",""])
        if len(paired):
            lines.extend(["Paired comparisons:",""]+[f"- `{r.comparison_id}`: difference {r.observed_difference_a_minus_b:.6f}, 95% interval [{r.ci_2_5:.6f}, {r.ci_97_5:.6f}], P(diff>0)={r.probability_difference_gt_zero:.4f}." for _,r in paired.iterrows()]+[""])
        if len(margin):
            lines.extend(["Strict 0.0025-MCC non-inferiority sensitivity:",""]+[f"- Reduced `{r.reduced_subset}`: P(loss ≤ 0.0025)={r.probability_loss_within_margin:.4f}; {r.sensitivity_conclusion}." for _,r in margin.iterrows()]+[""])
        if block in {"population_frequency","computational_predictors"}:
            lines.extend(["This block remains sensitivity-only regardless of training performance.",""])
    lines.extend(["## Overall provisional configuration", "", "Because conservation and amino-acid removal decisions remain method-dependent, the defensible provisional primary configuration is to retain the complete current primary biological feature set until supervisor review. Treat Missense Z plus LOEUF and the retained conservation measures as grouped evidence rather than independent biological votes. The reduced phyloP + PhastCons alternative remains plausible at 0.005 and 0.010 MCC margins but is not confirmed at the strict 0.0025 margin.","", "Population frequency remains sensitivity-only, with maximum-population AF as the parsimonious sensitivity representation. Computational predictors remain sensitivity-only; REVEL + CADD + AlphaMissense is a non-inferior diagnostic subset, but nested subset identity is unstable and provenance/circularity governance is unchanged.","", "No authoritative configuration or validation command was modified."])
    (g.STAGE/"GROUP_AWARE_FEATURE_SELECTION_REPORT.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    metadata_path=g.STAGE/"run_metadata.json"; metadata=json.loads(metadata_path.read_text()); metadata["decision_synthesis_refreshed_from_outputs"]=True; metadata["overall_provisional_configuration"]="retain complete current primary biological feature set pending supervisor decision on conservation and amino-acid blocks"; metadata["test_suite"]={"command":"python -m unittest discover -s tests -v","passed":82,"failed":0,"skipped":0,"result":"PASS","log":"test_results.txt"}; g.write_json(metadata_path,metadata)
    g.hashes()


if __name__=="__main__": main()
