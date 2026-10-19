# Static circularity-audit protocol

## Scope

This audit evaluates structural integrity, missingness shortcuts, within-group feature overlap and documented circularity risks without training a predictive model. Statistical flags identify candidates for controlled ablation; they do not prove leakage and do not automatically remove features.

## Partition policy

All label-dependent calculations use `data/processed/modeling_splits/train_core.csv` only. Validation is not used to choose features in this audit. Test files are opened only to verify row counts, schemas, unique keys, and the absence of gene or variant overlap. Test labels are not read for association analyses.

## Missingness

Blank, `.`, `NA`, `NaN`, `None` and `null` are treated as missing. Per-feature training-only missingness is reported overall and by class. An absolute class gap of at least 0.10 is flagged. Flags require modelling controls and do not cause automatic exclusion.

## Numeric overlap

Pearson and Spearman correlations use pairwise complete training observations. The paired count is reported. Constant, all-missing and insufficient pairs are marked invalid. An absolute coefficient of at least 0.85 is a redundancy candidate.

## Binary and categorical overlap

Corrected Cramér's V is computed from pairwise complete training observations. Constant, all-missing, insufficient and extremely sparse comparisons are explicitly identified. A corrected value of at least 0.80 is a redundancy candidate.

## Governance interpretation

The evidence-lineage registry distinguishes approved primary evidence, sensitivity-only evidence, prohibited fields and unresolved literature questions. Required later controls include no-population comparison, removal of gene-constraint evidence, evidence-group removal, a missingness-only diagnostic and optional-predictor experiments only after literature review. Scientific clearance remains conditional.
