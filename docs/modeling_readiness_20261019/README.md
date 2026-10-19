# Leakage-controlled modelling cohort

The frozen `dataset-v1.0-40000` source remains unchanged. Derived post-processing removes 448 rows whose ClinVar and VEP-reconciled genes differ, leaving 39,552 unique variants.

## Primary modelling rules

- `variant_key` is metadata only and must never be supplied to a model.
- `reconciled_gene_symbol` is metadata used only for gene-grouped splitting, gene-level bootstrapping and reporting. It must never be supplied as a predictor.
- `binary_label` is the target only.
- Models must select predictors from `data/processed/modeling_feature_registry.csv` where `primary_model_feature = 1`.
- Use `data/processed/modeling_ready_core_no_population.csv` as the primary leakage-controlled baseline.
- Treat `data/processed/modeling_ready_core.csv`, which includes population evidence, as a sensitivity experiment until neutral population handling is implemented.
- Fit every imputer, scaler, encoder, feature selector and model using training genes only.
- Do not add population-frequency availability or missingness indicators.
- Never pass raw population NaNs directly to Random Forest, XGBoost or another model capable of learning missing-value routing.
- Unavailable population evidence must contribute neutrally and must be reported as unavailable.
- The no-population ablation remains mandatory.
- Do not use the raw semicolon-separated InterPro or UniProt context strings.
- Do not use UniProt `Natural variant`, `Mutagenesis`, `Sequence conflict` or `Alternative sequence` tokens in primary models.
- Use validation genes for model selection and test genes once for the final locked evaluation.
- This readiness workflow prepares data and governance artifacts only; it does not train a model.

## Corrected split

| Split | Variants | Genes | Pathogenic proportion |
|---|---:|---:|---:|
| Train | 27,727 | 5,206 | 50.44% |
| Validation | 5,790 | 1,115 | 50.00% |
| Test | 6,035 | 1,115 | 50.19% |

Gene overlap is zero. The modelling-readiness audit passed with no critical failures. Population-frequency missingness remains label-associated, so the neutral train-only preprocessing rule and the no-population ablation are mandatory.

## Reproduce derived outputs

```bash
python scripts/12b_prepare_modeling_dataset.py

python scripts/13_split_by_gene.py \
  --full data/processed/modeling_ready_full.csv \
  --core data/processed/modeling_ready_core.csv \
  --optional data/processed/model_table_with_optional_predictors.csv \
  --no-population data/processed/modeling_ready_core_no_population.csv \
  --output-dir data/processed/modeling_splits \
  --assignment-output data/processed/modeling_split_gene_assignments.csv \
  --audit-output audits/13_modeling_gene_split_audit.json \
  --search-attempts 2000 \
  --minimum-evaluation-genes 1100

python scripts/14_audit_modeling_readiness.py
```

Generated datasets and split CSVs are intentionally ignored by Git. The accompanying audit and checksum manifest identify the verified local outputs.
