# Grouped additive model: frozen starting point

This directory records the starting point for grouped-additive model development. It does not contain validation or test results.

## Frozen primary condition

- Cohort: 39,552 variants with the existing gene-separated assignments.
- Development partition: 27,727 training variants only.
- Predictors: 48 approved no-population predictors.
- Evidence groups: conservation, gene constraint, amino-acid properties and protein context.
- Conservation: phyloP and PhastCons; GERP++ is excluded.
- Population frequency: sensitivity only.
- Composite computational predictors: external comparators or sensitivity only.
- Missingness indicators: excluded from the primary condition.

The machine-readable definition is `config/grouped_additive_primary_input.json`.

## Frozen experiment protocol

Model development uses the existing fixed five-fold gene-grouped training assignments. Preprocessing is fitted inside each fold using fold-training rows only. A logistic-regression baseline is implemented first to verify the shared modelling pipeline. The grouped additive model must produce one scalar logit contribution for each evidence group, and those contributions plus the intercept must reconstruct every saved model logit.

Validation remains unavailable during development. The test partition remains locked until one final architecture has been selected after the later validation comparison.

The machine-readable protocol is `config/grouped_additive_experiment_protocol.json`.

Any later change requires a numbered protocol amendment written before the affected validation or test outcomes are accessed.
