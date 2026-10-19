# Static circularity-audit protocol

## Scope

This audit evaluates structural integrity, missingness shortcuts, within-group feature overlap and documented circularity risks without training a predictive model. Statistical flags identify candidates for controlled ablation; they do not prove leakage and do not automatically remove features.

## Partition policy

All label-dependent calculations use the locked training partition only. Primary-feature analyses use `train_core_no_population.csv`; the two population fields are examined separately from `train_core.csv` as a sensitivity condition. Validation is not used to choose features. Test labels are not used for feature associations, feature selection or methodological decisions. Aggregate test class counts are accessed only from the pre-existing gene-assignment summary to verify the locked assignment; split CSVs are read structurally for schemas, keys, genes and overlap checks.

## Missingness

Blank, `.`, `NA`, `NaN`, `None` and `null` are treated as missing. Per-feature training-only missingness is reported overall and by class. An absolute class gap of at least 0.10 is flagged. Flags require modelling controls and do not cause automatic exclusion.

## Numeric overlap

Pearson and Spearman correlations use pairwise complete training observations. Primary correlations use the 49-feature no-population table. The two population fields are also correlated separately and labelled as a population-inclusive sensitivity result so their redundancy is documented before that experiment. The paired count is reported. Constant, all-missing and insufficient pairs are marked invalid. An absolute coefficient of at least 0.85 is a redundancy candidate.

## Binary and categorical overlap

Corrected Cramér's V is computed from pairwise complete training observations. Every binary feature has positive and negative counts and prevalence reported. A binary minority level is flagged when its exact count is below 20 or its unrounded proportion is below 1%; contingency tables are flagged as sparse when any observed level combination has fewer than five observations. These flags use level/cell counts rather than the total paired count. Constant, all-missing and insufficient comparisons are explicitly identified. A corrected value of at least 0.80 is a redundancy candidate. Sparse features are not automatically removed, but their associations require cautious interpretation.

## Mixed numeric–categorical overlap

Correlation ratio eta is computed for every primary numeric–categorical pair using pairwise complete training observations. Category counts and rare-level flags accompany eta so associations involving sparse levels can be interpreted cautiously. No sparse feature is automatically removed.

## Governance interpretation

The 49 predictors in `modeling_ready_core_no_population.csv` are the approved primary set. The two population-frequency fields and five optional computational predictors are sensitivity-only, and the primary and sensitivity sets must be disjoint. Required later controls include population-inclusive sensitivity comparison, removal of gene-constraint evidence, evidence-group removal, a missingness-only diagnostic and optional-predictor experiments only after literature review. Scientific clearance remains conditional.
