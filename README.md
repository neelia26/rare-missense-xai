# Rare Missense-SNV Pre-Modelling Evidence Audit

Branch: `feature/pre-modeling-evidence-audit`

This branch documents and implements the training-only evidence audit performed before formal model development. It starts from the already frozen rare-missense dataset and locked gene-separated split; it does not rebuild, rebalance, regenerate or resplit the dataset.

The purpose of this branch is to determine which evidence groups are scientifically and statistically suitable for controlled modelling, while keeping validation and test outcomes protected.

## Branch status

- Frozen source: 40,000 variants.
- Reconciled modelling cohort: 39,552 variants.
- Preserved gene-mapping exclusions: 448 variants.
- Train/validation/test: 27,727 / 5,790 / 6,035 variants.
- Split design: gene-disjoint and variant-disjoint.
- Feature analyses: training partition only.
- Validation outcomes: not inspected for feature selection or methodological decisions.
- Test outcomes: not inspected; the test partition remains locked.
- Formal model development: not started on this branch.

No final feature-lock decision from later review conversations is encoded here unless it appears in the committed audit outputs.

## What this branch adds

The audit is organised under `outputs/pre_modeling_evidence_audit/`:

| Stage | Purpose | Main outputs |
|---|---|---|
| `00_locked_protocol` | Verify the frozen cohort, schemas, preprocessing policy and locked split | Cohort reconciliation, split-integrity checks, leakage inventory and protocol checksum |
| `01_provenance_circularity` | Review where each evidence field came from and whether supervised-label or predictor-training overlap is possible | Provenance register, overlap checks, governance decisions and limitations |
| `02_redundancy` | Measure shared information without confusing correlation with circularity | Coverage, Pearson/Spearman correlations, VIF and redundancy figures |
| `03_controlled_experiments` | Compare predefined training-only evidence conditions using one fixed audit model | Fixed gene folds, OOF predictions, conservation experiments, predictor comparators and population sensitivity |
| `04_report` | Present the complete pre-modelling evidence review | Revised audit report, supervisor summary, figure index and initial hash manifest |
| `05_dominance_hierarchical_partitioning` | Distribute cross-validated MCC gain across every within-block feature subset | Exact Shapley/general-dominance contributions, paired comparisons and 5,000-gene bootstrap intervals |
| `05_group_aware_feature_selection` | Test whether candidate subsets are stable across unseen genes | Fold-local mRMR, nested gene-grouped wrapper selection and non-inferiority sensitivity analyses |

## Evidence blocks examined

The controlled analyses considered:

- conservation: GERP++, phyloP and PhastCons;
- amino-acid properties: hydropathy delta and polarity delta;
- gene constraint: missense Z and LOEUF;
- population sensitivity: global AF and maximum-population AF; and
- computational-predictor sensitivity: SIFT, PolyPhen-2, REVEL, CADD and AlphaMissense.

Population features and composite computational predictors remain governed sensitivity/comparator conditions in the committed audit. Predictive performance alone does not override missingness, provenance or circularity concerns.

## Common ML controls

All controlled feature experiments used the same audit instrument:

- fixed five-fold cross-validation grouped by reconciled gene;
- fold-local median imputation;
- fold-local standard scaling;
- logistic regression with `C=1`, `liblinear`, balanced class weights and a fixed seed;
- one out-of-fold probability for every training variant;
- pooled OOF MCC at threshold 0.5 as the primary metric; and
- PR-AUC, ROC-AUC and Brier score as supporting metrics.

This logistic regression is a controlled comparison instrument, not the final thesis model.

## Dominance and hierarchical partitioning

Every possible candidate subset was evaluated within each evidence block while the remaining governed biological features were held fixed. Exact Shapley/general-dominance decomposition then allocated the block's full-minus-empty OOF MCC gain across its candidates.

The implementation verifies that raw contributions sum exactly to the observed block gain. Normalised shares sum to 100% by construction and must not be interpreted as causal effects, probabilities of necessity or proof of biological independence.

See:

- `outputs/pre_modeling_evidence_audit/05_dominance_hierarchical_partitioning/ANALYSIS_SPECIFICATION.md`
- `outputs/pre_modeling_evidence_audit/05_dominance_hierarchical_partitioning/DOMINANCE_HIERARCHICAL_PARTITIONING_REPORT.md`

## Group-aware feature selection

The second feature-selection check asked whether conclusions remained stable when applied to genes excluded from model fitting. It combined:

- fold-local minimum-redundancy maximum-relevance (mRMR) ranking;
- nested gene-grouped evaluation of every within-block subset;
- fixed-subset outer-fold OOF comparisons;
- 5,000-sample paired gene bootstrap intervals; and
- non-inferiority sensitivity margins of 0.0025, 0.005 and 0.010 MCC.

Dominance and group-aware selection answer different questions. Dominance allocates fixed-model performance gain; group-aware selection evaluates subset stability and generalisation across genes.

See:

- `outputs/pre_modeling_evidence_audit/05_group_aware_feature_selection/ANALYSIS_SPECIFICATION.md`
- `outputs/pre_modeling_evidence_audit/05_group_aware_feature_selection/GROUP_AWARE_FEATURE_SELECTION_REPORT.md`
- `outputs/pre_modeling_evidence_audit/05_group_aware_feature_selection/METHOD_COMPARISON.md`

## Primary reports

- `outputs/pre_modeling_evidence_audit/04_report/PRE_MODELING_EVIDENCE_AUDIT_REVISED.md`
- `outputs/pre_modeling_evidence_audit/04_report/SUPERVISOR_PROGRESS_SUMMARY.md`
- `outputs/pre_modeling_evidence_audit/04_report/FIGURE_INDEX.md`
- `docs/circularity_audit_summary.md`

Machine-readable results are preserved beside the reports as CSV and JSON files. Audit figures are committed as PNG files.

## Reproduce and verify

From the repository root, using the project environment:

```bash
source .venv/bin/activate
python -m src.pre_modeling_audit.validate_final_cohort
python -m src.pre_modeling_audit.run_pre_modeling_audit
python -m src.pre_modeling_audit.dominance_hierarchical_partitioning
python -m src.pre_modeling_audit.group_aware_feature_selection
python -m unittest discover -s tests -v
```

These analysis commands regenerate audit outputs. They do not authorise validation or test inspection. Do not execute a saved validation-confirmation command until the training-only configuration has been formally reviewed and locked.

## Why this README differs across branches

Each long-lived research branch represents a distinct reproducible milestone and should explain that milestone rather than copy one generic README:

| Branch | Intended landing-page focus |
|---|---|
| `feature/data-pipeline` | Numbered dataset-construction stages and governed source integration |
| `feature/full-api-automation` | Resumable API execution, checkpoints and the completed full dataset rebuild |
| `feature/circularity-audit` | Evidence lineage, leakage controls, feature governance and circularity findings |
| `feature/pre-modeling-evidence-audit` | Training-only evidence review, controlled experiments and feature-selection audits |
| `feature/modeling` | The clean modelling foundation after configurations are formally locked |

Shared installation or source documentation should live in dedicated files and be linked where needed. The branch README should remain a concise scientific landing page: what changed, why it exists, what is complete, what remains locked and where the authoritative artifacts live.
