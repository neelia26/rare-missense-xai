# Training-only dominance and hierarchical-partitioning analysis specification

## Status and timing

This specification is written and SHA-256 locked before the new subset experiments are executed. It extends the completed pre-modelling evidence audit without modifying the frozen dataset, locked split, original experiment results, authoritative recommendation, or validation command.

## Method

The method is **Exact performance-based general dominance analysis using Shapley decomposition of gene-grouped cross-validated MCC**. It is a classification-performance analogue of hierarchical partitioning/general dominance analysis; it is not classical R-squared dominance analysis.

For candidate set `F`, fixed base features `B`, and candidate subset `S`, the value function is:

`v(S) = MCC from one out-of-fold prediction per training variant, using B union S.`

For each candidate feature `j`, its exact contribution is:

`phi_j = sum_{S subset F without j} [|S|! (p-|S|-1)! / p!] [v(S union {j}) - v(S)]`.

The implementation must verify `sum(phi_j) = v(F) - v(empty)` within `1e-10`. When the total gain is positive and non-zero, `normalised_share_j = phi_j / [v(F)-v(empty)]` and shares must sum to one within `1e-10`. Negative contributions are retained. If total gain is zero or negative, normalisation is undefined; raw contributions remain reported.

## Common controls

- Training partition only: 27,727 variants.
- Existing fixed five-fold assignments grouped by `reconciled_gene_symbol`.
- Seed: 20261018.
- Fixed logistic-regression audit model: `C=1`, `solver=liblinear`, `class_weight=balanced`, `max_iter=2000`.
- Fold-local median imputation and standard scaling, fitted only on each fitting fold.
- Existing governed feature tables and missing-value representation; missing biological evidence is never changed to zero.
- One out-of-fold probability per variant per subset.
- Primary value: aggregate out-of-fold MCC at threshold 0.5.
- Supporting metrics: average precision (reported as PR-AUC), ROC-AUC and Brier score.
- No model-family search, hyperparameter search, validation use, or test use.
- Existing conservation OOF predictions may be reused only after exact key, fold, label, feature, model and preprocessing equivalence checks.

## Governed blocks

1. **Primary conservation:** candidates GERP++, phyloP and PhastCons; fixed base is the primary no-population biological set without these candidates; all eight subsets.
2. **Primary amino-acid related pair:** hydropathy delta and polarity delta; fixed base is the primary no-population set without them; all four subsets.
3. **Primary gene-constraint related pair:** gnomAD missense Z and LOEUF; fixed base is the primary no-population set without them; all four subsets.
4. **Population sensitivity:** global AF and maximum-population AF; fixed base is the governed primary no-population biological model; all four subsets, labelled `population_inclusive_sensitivity`. Missing frequency is not absence and is never set to zero. Class-associated missingness remains a governance warning.
5. **Computational-predictor sensitivity:** SIFT, PolyPhen-2 HDIV, REVEL, CADD and AlphaMissense; fixed base is the primary no-population biological configuration with all three conservation scores; all 32 subsets. Contributions are conditional on this base and cannot override provenance or circularity governance.

## Paired gene-cluster bootstrap

Each completed block uses 5,000 deterministic paired bootstrap samples. Genes are sampled with replacement and every variant belonging to a sampled gene is retained with the gene's multiplicity. The same resample is applied to all subsets within a block. Subset MCC values and exact Shapley contributions are recalculated for every sample. Outputs include percentile 95% intervals, medians, and probabilities that defined paired performance differences are greater than zero.

## Decision restrictions

Dominance quantifies predictive contribution under this fixed audit model. It does not establish biological independence, causal importance, or absence of circularity. Correlation is shared information, not circularity. A feature is not removed solely because it has the smallest share, and no post-result percentage threshold will be introduced. Provisional interpretations must jointly consider paired uncertainty, magnitude of loss, parsimony, biological distinction, missingness, provenance, and interpretability. Inconclusive evidence remains explicitly inconclusive.

SIFT, PolyPhen-2 and REVEL remain external comparators. CADD and AlphaMissense remain optional sensitivity inputs. Population evidence remains sensitivity-only. No result from this stage modifies the authoritative configuration or triggers validation.
