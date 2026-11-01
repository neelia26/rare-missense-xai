# Pre-modelling evidence audit

## Purpose and scope

This work separates (A) reproducible cohort/preprocessing verification, (B) a pre-modelling provenance and redundancy audit, and (C) preliminary modelling for evidence selection. Phase C is not exploratory data analysis and is not the final model comparison.

## Dataset and split confirmation

The frozen dataset contains 40,000 unique variants. The reconciled modelling cohort contains 39,552, while 448 gene-mapping disagreements remain preserved separately. Retained and excluded keys are disjoint and account exactly for the frozen source. The existing train/validation/test assignment remains gene- and variant-disjoint. No source dataset or split was modified. Aggregate validation/test class counts were used only from the locked assignment summary for integrity verification; outcome rows were not analysed.

## Preprocessing and leakage

The script inventory, schema, missingness and leakage reports are under `00_locked_protocol/cohort_validation`. Target, clinical text, identifiers, gene/transcript identifiers, database IDs and lookup metadata are excluded. Natural variant, Mutagenesis, Sequence conflict and Alternative sequence indicators remain prohibited. Fold-local median imputation and scaling were used only in preliminary training CV; missing biology was never converted to zero.

## Circularity versus redundancy

Circularity concerns reuse of clinical labels/variants or downstream information. Redundancy concerns shared information between features. Low correlation cannot prove absence of circularity, and high correlation does not by itself prove it. Heatmaps are therefore presented only as training-data shared-variance summaries.

## Provenance findings

GERP++, phyloP and phastCons are evolutionary measures and are classified as not applicable to supervised-label circularity, although they can double-count evolutionary information. PolyPhen-2 HDIV uses supervised disease/divergence labels and UniProt/Swiss-Prot annotations. REVEL uses an HGMD-derived pathogenic training set and incorporates SIFT, PolyPhen and multiple conservation scores; exact overlap is unverifiable because the full pathogenic set is restricted. CADD uses proxy labels rather than clinical labels but overlaps broad annotations. AlphaMissense was evaluated on ClinVar rather than directly trained on it, while still sharing evolutionary/protein and weak population information. Exact overlap was not forced where safe harmonisation or public records were unavailable.

## Coverage and redundancy

Training-only coverage, class-stratified coverage, pairwise counts, Spearman/Pearson matrices, VIF and flagged pairs are under `02_redundancy`. Score direction matters: lower SIFT indicates greater predicted damage, so negative correlations may still represent shared information. Flags trigger investigation, not automatic deletion.

## Preliminary controlled experiments

All eight conservation subsets were compared with a fixed logistic-regression audit instrument using identical five-fold gene-grouped training-only folds. Each optional predictor was evaluated alone and as a one-at-a-time addition. The population-inclusive condition was sensitivity-only. No model-family or hyperparameter search was performed.

The training-only recommended condition is `conservation_gerp_rs__phylop100way__phastcons100way` with 49 features. This recommendation is provisional and must be confirmed once on validation using the saved command. Optional predictors remain external comparators/sensitivity inputs because predictive improvement cannot resolve circularity. The no-population condition remains primary.

## Limitations

Exact external training overlap remains unavailable or unsafe to harmonise for several scores. The reference model is an audit instrument, not the final thesis model. No causal or non-circularity conclusion follows from correlations, PCA (not used), or performance.

## Commands and next step

Main command: `python -m src.pre_modeling_audit.run_pre_modeling_audit`. The exact future validation command is in `03_controlled_experiments/validation_confirmation_command.txt`; it was not executed. Validation outcomes were not inspected. Test outcomes were not inspected and the test partition remains locked.

## Output index

Outputs follow the requested `00_locked_protocol`, `01_provenance_circularity`, `02_redundancy`, `03_controlled_experiments`, and `04_report` directories. CSV/JSON files contain machine-readable evidence; PNG files are correlation visualisations only.


## Complete execution log

```text
git status --short
git branch --show-current
git worktree list
python -m py_compile src/pre_modeling_audit/run_pre_modeling_audit.py src/pre_modeling_audit/validate_final_cohort.py tests/test_pre_modeling_audit.py
python -m unittest tests.test_pre_modeling_audit -v
python -m src.pre_modeling_audit.validate_final_cohort
python -m src.pre_modeling_audit.run_pre_modeling_audit
python -m unittest discover -s tests -v
```
