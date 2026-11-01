# Group-aware correlated-feature selection specification

## Scope and safeguards

This specification is fixed before analysis. The work uses only the 27,727-row training partition and the existing five fixed gene-grouped outer folds. Validation and test files, rows, labels, predictions and outcomes are prohibited. The frozen dataset, locked split, authoritative configuration and validation command are immutable. No model-family or hyperparameter search is permitted.

Hard path guards reject filenames containing `validation` or `test`. Only the three governed training tables and pre-existing training fold assignments may be loaded.

## Relationship to dominance analysis

The existing method is **Exact all-subsets Shapley/general-dominance decomposition of training-only cross-validated MCC.** It describes distribution of predictive contribution across correlated features. Its normalised shares sum to 100% by construction; they are not probabilities of feature necessity and do not establish biological independence or causality. Dominance is supplementary evidence, not the sole keep/drop criterion.

## Candidate blocks

- Primary conservation: GERP++, phyloP, PhastCons.
- Primary amino-acid properties: hydropathy delta, polarity delta.
- Primary gene constraint: missense Z, LOEUF.
- Population-frequency sensitivity: global AF, maximum-population AF.
- Computational-predictor sensitivity: SIFT, PolyPhen-2 HDIV, REVEL, CADD, AlphaMissense.

Population and computational-predictor blocks remain sensitivity-only regardless of performance. Computational predictors retain existing provenance and circularity governance.

## Fold-local mRMR sensitivity analysis

Within each outer-training portion, candidate relevance is estimated with deterministic mutual information against the training outcome. Pairwise candidate redundancy is estimated with deterministic mutual information between fold-locally median-imputed candidate features. A greedy mRMR ranking uses relevance minus mean redundancy with already ranked features. Rankings are descriptive sensitivity evidence; every non-empty subset is evaluated by the wrapper, so the analysis does not rely on the greedy ranking alone.

## Nested gene-grouped wrapper

The existing fixed folds are outer folds. For each outer fold, four deterministic `StratifiedGroupKFold` inner folds are created exclusively from outer-training genes. Every candidate subset, including the empty subset, is evaluated by pooled inner OOF MCC using the fixed logistic audit model (`C=1`, `liblinear`, balanced class weights, seed 20261018). Median imputation and standard scaling are fitted independently inside every inner fitting fold. The highest inner MCC selects the subset, with fewer candidates and lexical subset key as deterministic tie-breakers. The selected subset is then fitted on the complete outer-training portion and evaluated once on the untouched outer holdout. Combined outer OOF predictions provide MCC, average precision, ROC-AUC and Brier score.

Fixed-subset outer OOF results are also reported for all subsets. These comparisons use the same pre-specified subset in each outer fold and therefore do not use outer outcomes for selection.

## Stability and uncertainty

Paired uncertainty uses 5,000 deterministic gene-cluster bootstrap samples of fixed-subset outer OOF predictions, preserving all variants of sampled genes. Selection frequencies are reported across outer folds and bootstrap samples. Non-inferiority sensitivity is reported at MCC margins 0.0025, 0.005 and 0.010 without locking a configuration.

## Decision framework

Each block receives one provisional status: retain complete block; retain reduced subset; retain as grouped representation; sensitivity-only; or inconclusive—requires supervisor decision. Conclusions jointly consider biology, provenance/circularity, mRMR relevance and redundancy, nested wrapper performance, paired uncertainty and stability. No arbitrary post-result removal threshold is introduced, no dominance percentage is used as a final rule, and correlated retained features are interpreted as an evidence group rather than independent biological votes.
