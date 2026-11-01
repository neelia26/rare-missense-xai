# Rare Missense-SNV Pre-Modelling Evidence Audit:
## Dominance and Hierarchical-Partitioning Addendum

## Scope

This is an **Exact all-subsets Shapley/general-dominance decomposition of training-only cross-validated MCC**. It uses only 27,727 training variants and the existing five gene-grouped folds. It is a predictive-contribution analysis, not classical R-squared dominance analysis. Normalised shares sum to 100% by construction; they are not probabilities of feature necessity and do not establish biological independence or causality. Dominance is supplementary evidence rather than the sole keep/drop criterion. Validation and test outcomes were not loaded, inspected, or evaluated.

Exact Shapley efficiency and valid normalised-share identities were verified for every block. Negative contributions, if present, are preserved. Supporting PR-AUC values are average precision.


## Conservation MCC reconciliation

The previous audit reported the unweighted arithmetic mean of five fold-level MCC values. The dominance stage reports one MCC calculated after pooling all 27,727 OOF predictions. Rows, order, labels, folds, features, preprocessing, model parameters, threshold, predictions and software versions were identical. Therefore the discrepancy is solely the nonlinearity of MCC under mean-of-fold versus pooled-OOF aggregation. See `conservation_mcc_reconciliation.csv` for condition-level values and exact checks.

## Primary Conservation

Best observed subset: `GERP++ + phyloP + PhastCons` (MCC 0.721727). Full candidate subset MCC: 0.721727; empty-candidate base MCC: 0.551119; total gain: 0.170608.

| Feature | Raw contribution | Normalised share | Bootstrap median | 95% interval |
|---|---:|---:|---:|---:|
| GERP++ | 0.041675 | 0.244 | 0.041634 | [0.034225, 0.049609] |
| phyloP | 0.086702 | 0.508 | 0.086301 | [0.070331, 0.105433] |
| PhastCons | 0.042230 | 0.248 | 0.042227 | [0.035552, 0.049273] |

Coverage and missingness:

- GERP++: 27,687/27,727 nonmissing (99.86%); missing values were not converted to zero.
- phyloP: 27,726/27,727 nonmissing (100.00%); missing values were not converted to zero.
- PhastCons: 27,726/27,727 nonmissing (100.00%); missing values were not converted to zero.

Paired full/reduced comparisons:

| Comparison | Observed MCC difference | Bootstrap median | 95% interval | P(diff > 0) |
|---|---:|---:|---:|---:|
| full_vs_without_conservation_gerp_rs | 0.001315 | 0.001296 | [-0.000772, 0.003411] | 0.888 |
| full_vs_without_conservation_phylop100way | 0.036976 | 0.036792 | [0.025340, 0.049517] | 1.000 |
| full_vs_without_conservation_phastcons100way | 0.003598 | 0.003585 | [0.000314, 0.007055] | 0.983 |
| all3_vs_phylop_phastcons | 0.001315 | 0.001296 | [-0.000772, 0.003411] | 0.888 |
| phylop_phastcons_vs_phylop | 0.004998 | 0.005032 | [0.001286, 0.008965] | 0.994 |

Governance: primary biological evidence. Dominance is conditional predictive contribution, not causality, independence, or proof against circularity.
Provisional interpretation: phyloP + PhastCons remains the parsimonious alternative because the paired full-minus-reduced MCC interval includes zero. This does not modify the authoritative configuration.

## Primary Amino Acid

Best observed subset: `Hydropathy delta + Polarity delta` (MCC 0.721727). Full candidate subset MCC: 0.721727; empty-candidate base MCC: 0.718722; total gain: 0.003005.

| Feature | Raw contribution | Normalised share | Bootstrap median | 95% interval |
|---|---:|---:|---:|---:|
| Hydropathy delta | 0.001894 | 0.630 | 0.001891 | [-0.000412, 0.004149] |
| Polarity delta | 0.001111 | 0.370 | 0.001115 | [-0.000628, 0.002856] |

Coverage and missingness:

- Hydropathy delta: 27,727/27,727 nonmissing (100.00%); missing values were not converted to zero.
- Polarity delta: 27,727/27,727 nonmissing (100.00%); missing values were not converted to zero.

Paired full/reduced comparisons:

| Comparison | Observed MCC difference | Bootstrap median | 95% interval | P(diff > 0) |
|---|---:|---:|---:|---:|
| full_vs_without_aa_hydropathy_delta | 0.001017 | 0.001020 | [-0.001213, 0.003174] | 0.817 |
| full_vs_without_aa_polarity_delta | 0.000234 | 0.000224 | [-0.000562, 0.001076] | 0.716 |
| aa_hydropathy_delta_vs_aa_polarity_delta | 0.000784 | 0.000788 | [-0.001496, 0.003076] | 0.746 |

Governance: primary related-feature analysis. Dominance is conditional predictive contribution, not causality, independence, or proof against circularity.
Provisional interpretation: hydropathy alone is the more parsimonious observed single-feature option, and adding polarity has no clearly supported MCC gain. Because both dominance intervals include zero and the descriptors are biologically distinct, removal is not finalised; the evidence is supportive but inconclusive.

## Primary Gene Constraint

Best observed subset: `Missense Z + LOEUF` (MCC 0.721727). Full candidate subset MCC: 0.721727; empty-candidate base MCC: 0.712218; total gain: 0.009509.

| Feature | Raw contribution | Normalised share | Bootstrap median | 95% interval |
|---|---:|---:|---:|---:|
| Missense Z | 0.004743 | 0.499 | 0.004444 | [-0.002086, 0.012619] |
| LOEUF | 0.004766 | 0.501 | 0.004805 | [-0.000947, 0.010308] |

Coverage and missingness:

- Missense Z: 25,492/27,727 nonmissing (91.94%); missing values were not converted to zero.
- LOEUF: 25,434/27,727 nonmissing (91.73%); missing values were not converted to zero.

Paired full/reduced comparisons:

| Comparison | Observed MCC difference | Bootstrap median | 95% interval | P(diff > 0) |
|---|---:|---:|---:|---:|
| full_vs_without_gnomad_missense_z | 0.009883 | 0.009626 | [0.000397, 0.020483] | 0.979 |
| full_vs_without_gnomad_loeuf | 0.009905 | 0.009933 | [0.001847, 0.017640] | 0.994 |
| gnomad_missense_z_vs_gnomad_loeuf | -0.000023 | -0.000230 | [-0.008465, 0.008711] | 0.479 |

Governance: controlled related-feature analysis. Dominance is conditional predictive contribution, not causality, independence, or proof against circularity.
Provisional interpretation: retain both constraint measures for now. Removing either from the full pair produced a positive full-minus-reduced interval, while the two singleton models were indistinguishable; their contributions are balanced and biologically non-identical.

## Population Inclusive Sensitivity

Best observed subset: `Max-population AF` (MCC 0.774576). Full candidate subset MCC: 0.774505; empty-candidate base MCC: 0.721727; total gain: 0.052778.

| Feature | Raw contribution | Normalised share | Bootstrap median | 95% interval |
|---|---:|---:|---:|---:|
| Global AF | 0.020193 | 0.383 | 0.020222 | [0.017022, 0.023423] |
| Max-population AF | 0.032585 | 0.617 | 0.032649 | [0.027033, 0.038312] |

Coverage and missingness:

- Global AF: 16,567/27,727 nonmissing (59.75%); missing values were not converted to zero.
- Max-population AF: 16,567/27,727 nonmissing (59.75%); missing values were not converted to zero.

Paired full/reduced comparisons:

| Comparison | Observed MCC difference | Bootstrap median | 95% interval | P(diff > 0) |
|---|---:|---:|---:|---:|
| full_vs_without_gnomad_global_af | -0.000071 | -0.000071 | [-0.000785, 0.000639] | 0.423 |
| full_vs_without_gnomad_max_population_af | 0.012322 | 0.012357 | [0.008516, 0.016212] | 1.000 |
| gnomad_global_af_vs_gnomad_max_population_af | -0.012392 | -0.012459 | [-0.016323, -0.008540] | 0.000 |

Governance: sensitivity only; class-associated missingness warning. Dominance is conditional predictive contribution, not causality, independence, or proof against circularity.
Provisional interpretation: maximum-population AF alone is the parsimonious sensitivity condition; adding global AF has no supported incremental MCC gain. Population evidence nevertheless remains sensitivity-only. Missing values were imputed fold-locally, never converted to zero, and class-associated missingness remains a warning.

## Computational Predictor Sensitivity

Best observed subset: `SIFT + PolyPhen-2 + REVEL + CADD + AlphaMissense` (MCC 0.875863). Full candidate subset MCC: 0.875863; empty-candidate base MCC: 0.721727; total gain: 0.154136.

| Feature | Raw contribution | Normalised share | Bootstrap median | 95% interval |
|---|---:|---:|---:|---:|
| SIFT | 0.012480 | 0.081 | 0.012479 | [0.010432, 0.014487] |
| PolyPhen-2 | 0.011329 | 0.073 | 0.011320 | [0.008926, 0.013636] |
| REVEL | 0.055668 | 0.361 | 0.055882 | [0.047020, 0.064540] |
| CADD | 0.024015 | 0.156 | 0.024012 | [0.018758, 0.028834] |
| AlphaMissense | 0.050644 | 0.329 | 0.050626 | [0.043864, 0.057494] |

Coverage and missingness:

- SIFT: 26,890/27,727 nonmissing (96.98%); missing values were not converted to zero.
- PolyPhen-2: 25,575/27,727 nonmissing (92.24%); missing values were not converted to zero.
- REVEL: 27,249/27,727 nonmissing (98.28%); missing values were not converted to zero.
- CADD: 27,726/27,727 nonmissing (100.00%); missing values were not converted to zero.
- AlphaMissense: 27,553/27,727 nonmissing (99.37%); missing values were not converted to zero.

Paired full/reduced comparisons:

| Comparison | Observed MCC difference | Bootstrap median | 95% interval | P(diff > 0) |
|---|---:|---:|---:|---:|
| full_vs_without_predictor_sift_score | 0.000076 | 0.000077 | [-0.000699, 0.000867] | 0.601 |
| full_vs_without_predictor_polyphen2_hdiv_score | 0.000503 | 0.000495 | [-0.000788, 0.001858] | 0.768 |
| full_vs_without_predictor_revel_score | 0.025600 | 0.025619 | [0.018402, 0.033150] | 1.000 |
| full_vs_without_predictor_cadd_phred | 0.005996 | 0.005950 | [0.002997, 0.008931] | 1.000 |
| full_vs_without_predictor_alphamissense_score | 0.020159 | 0.020033 | [0.014699, 0.025655] | 1.000 |
| all_five_vs_revel_cadd_alphamissense | 0.000215 | 0.000216 | [-0.001060, 0.001530] | 0.618 |

Governance: external comparator/optional sensitivity only. Dominance is conditional predictive contribution, not causality, independence, or proof against circularity.
Provisional interpretation: REVEL, CADD and AlphaMissense show supported conditional full-model losses when removed; SIFT and PolyPhen-2 do not show a clearly supported incremental loss from the full five-predictor subset. This is inconclusive for a governed reduced subset and does not justify automatic inclusion or exclusion. Predictors remain external comparators or optional sensitivity inputs; performance cannot override provenance or circularity governance.

## Decision boundary

These results are provisional training-only supporting evidence. They do not alter `recommended_training_only_configuration.json`, any validation command, or the locked split. Correlation indicates shared information, not circularity. No validation confirmation, formal baseline modelling, or final model selection was run.
