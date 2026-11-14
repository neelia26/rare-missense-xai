# Grouped additive model: frozen starting point

This directory records the starting point for grouped-additive model development. It does not contain validation or test results.

[Open the executed baseline notebook in Google Colab](https://colab.research.google.com/github/neelia26/rare-missense-xai/blob/feature/grouped-additive-model/notebooks/grouped_additive_baseline.ipynb)

## Frozen primary condition

- Cohort: 39,552 variants with the existing gene-separated assignments.
- Development partition: 27,727 training variants only.
- Predictors: 47 frozen no-population predictors.
- Evidence groups: conservation, gene constraint, amino-acid properties and protein context.
- Conservation: phyloP and PhastCons; GERP++ is excluded.
- Amino-acid properties: hydropathy delta is retained; continuous polarity delta is excluded; the distinct polarity-class-change indicator is retained.
- Population frequency: sensitivity only.
- Composite computational predictors: external comparators or sensitivity only.
- Missingness indicators: excluded from the primary condition.

The machine-readable definition is `config/grouped_additive_primary_input.json`.

## Frozen experiment protocol

Model development uses the existing fixed five-fold gene-grouped training assignments. Preprocessing is fitted inside each fold using fold-training rows only. A logistic-regression baseline is implemented first to verify the shared modelling pipeline. The grouped additive model must produce one scalar logit contribution for each evidence group, and those contributions plus the intercept must reconstruct every saved model logit.

Validation remains unavailable during development. The test partition remains locked until one final architecture has been selected after the later validation comparison.

The machine-readable protocol is `config/grouped_additive_experiment_protocol.json`.

Any later change requires a numbered protocol amendment written before the affected validation or test outcomes are accessed.

The registered population sensitivity uses maximum-population AF, an availability-only diagnostic and a same-covered-variants comparison. This is intended to distinguish genuine frequency information from the strong class-associated database-coverage pattern.

The first training-only baseline and its saved outputs are documented in `notebooks/grouped_additive_baseline.ipynb` and `outputs/grouped_additive_model/01_logistic_baseline/`.

## First grouped-additive training-only result

The compact grouped-additive search selected one hidden layer with eight tanh units per evidence group. Its pooled training OOF MCC was 0.7452, compared with 0.7221 for the fixed logistic baseline. The model saves one signed contribution for each of the four evidence groups for every OOF variant. The largest observed difference between a saved logit and the intercept-plus-contributions reconstruction was below 0.000001.

[Open the executed grouped-additive explanation notebook in Google Colab](https://colab.research.google.com/github/neelia26/rare-missense-xai/blob/feature/grouped-additive-model/notebooks/grouped_additive_model.ipynb)

These are training-only development results. They do not select a final project architecture and do not contain validation or test outcomes.

## Evidence-reporting layer

`scripts/generate_evidence_reports.py` converts each faithful four-group contribution into a deterministic report containing the observed biological annotations, the direction and relative strength of each model contribution, missing evidence and agreement or conflict between groups. The text explicitly separates model behaviour from causality or independent clinical evidence. An LLM is not used to generate these reports.
