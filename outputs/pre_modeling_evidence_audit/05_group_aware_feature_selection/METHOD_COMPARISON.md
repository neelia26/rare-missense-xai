# Method comparison

| Method | Question answered | Strength | Limitation | Decision role |
|---|---|---|---|---|
| Exact all-subsets Shapley/general-dominance decomposition of training-only cross-validated MCC | How is fixed-model MCC gain distributed across correlated features? | Exact additive predictive-contribution accounting | Shares total 100% by construction and are not probabilities of necessity; no causality or independence | Supplementary |
| Fold-local mRMR | Which candidates combine outcome relevance with lower candidate redundancy? | Training-fold ranking sensitivity | Ranking can change among correlated alternatives; not a final keep/drop proof | Supplementary ranking/stability |
| Nested gene-grouped wrapper | Which subset is selected without using the outer holdout? | Separates selection from outer assessment and respects genes | Results remain model- and sample-specific | Primary training-only subset evidence |
| Gene bootstrap/non-inferiority sensitivity | Are paired losses stable across genes and plausible margins? | Paired uncertainty and stability | Margins are declared sensitivity values, not automatic scientific thresholds | Decision support |

Correlated retained features are interpreted as an evidence group, not independent biological votes. Population and computational predictors remain sensitivity-only regardless of these results.
