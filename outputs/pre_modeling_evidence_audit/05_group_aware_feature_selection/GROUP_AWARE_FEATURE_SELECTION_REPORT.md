# Group-Aware Feature Selection Report

All results are **training only**. Validation and test outcomes were not accessed. Dominance remains supplementary predictive-contribution evidence, not the sole selection rule.

## Methods

Fold-local mRMR estimated outcome relevance and candidate redundancy separately within each outer-training portion. Every subset was then evaluated with four inner gene-grouped folds; selection never used the outer holdout. Fixed-subset outer OOF predictions supported paired 5,000-gene bootstrap comparisons and non-inferiority sensitivity margins of 0.0025, 0.005 and 0.010 MCC.

## Conservation

**Provisional decision:** inconclusive—requires supervisor decision

**Provisional subset/representation:** all three vs phyloP + PhastCons unresolved

Pair improves on phyloP, but strict-margin non-inferiority probability=0.8656 and pair outer-selection frequency=0.40; methods disagree.

Nested-selected pooled OOF metrics: MCC 0.721791, PR-AUC 0.925600, ROC-AUC 0.932953, Brier 0.101117.

mRMR mean ranking: phyloP 1.0; PhastCons 2.0; GERP++ 3.0.

Outer-fold subset selections: GERP++ + phyloP + PhastCons 3/5; phyloP + PhastCons 2/5.

Paired comparisons:

- `phylop_phastcons_vs_phylop`: difference 0.004998, 95% interval [0.001286, 0.008965], P(diff>0)=0.9936.
- `all3_vs_phylop_phastcons`: difference 0.001315, 95% interval [-0.000772, 0.003411], P(diff>0)=0.8878.

Strict 0.0025-MCC non-inferiority sensitivity:

- Reduced `conservation_phastcons100way__conservation_phylop100way`: P(loss ≤ 0.0025)=0.8656; not established at 95% probability.

## Amino Acid Properties

**Provisional decision:** inconclusive—requires supervisor decision

**Provisional subset/representation:** hydropathy/polarity/full alternatives unstable

mRMR ranks hydropathy first, but nested selections and gene-bootstrap preferences do not identify a stable single-feature replacement.

Nested-selected pooled OOF metrics: MCC 0.718915, PR-AUC 0.924169, ROC-AUC 0.932129, Brier 0.101683.

mRMR mean ranking: Hydropathy delta 1.0; Polarity delta 2.0.

Outer-fold subset selections: Polarity delta 2/5; nan 2/5; Hydropathy delta + Polarity delta 1/5.

Paired comparisons:

- `both_vs_hydropathy`: difference 0.000234, 95% interval [-0.000562, 0.001076], P(diff>0)=0.7160.
- `both_vs_polarity`: difference 0.001017, 95% interval [-0.001213, 0.003174], P(diff>0)=0.8170.
- `hydropathy_vs_polarity`: difference 0.000784, 95% interval [-0.001496, 0.003076], P(diff>0)=0.7464.

Strict 0.0025-MCC non-inferiority sensitivity:

- Reduced `aa_hydropathy_delta`: P(loss ≤ 0.0025)=1.0000; reduced non-inferior in sensitivity analysis.
- Reduced `aa_polarity_delta`: P(loss ≤ 0.0025)=0.9046; not established at 95% probability.

## Gene Constraint

**Provisional decision:** retain as grouped representation

**Provisional subset/representation:** Missense Z + LOEUF

Both full-versus-single intervals support joint value; all five nested folds selected the pair. They are correlated but biologically complementary.

Nested-selected pooled OOF metrics: MCC 0.721727, PR-AUC 0.925531, ROC-AUC 0.932965, Brier 0.101120.

mRMR mean ranking: Missense Z 1.0; LOEUF 2.0.

Outer-fold subset selections: Missense Z + LOEUF 5/5.

Paired comparisons:

- `both_vs_missense_z`: difference 0.009905, 95% interval [0.001847, 0.017640], P(diff>0)=0.9940.
- `both_vs_loeuf`: difference 0.009883, 95% interval [0.000397, 0.020483], P(diff>0)=0.9790.
- `missense_z_vs_loeuf`: difference -0.000023, 95% interval [-0.008465, 0.008711], P(diff>0)=0.4790.

Strict 0.0025-MCC non-inferiority sensitivity:

- Reduced `gnomad_missense_z`: P(loss ≤ 0.0025)=0.0352; not established at 95% probability.
- Reduced `gnomad_loeuf`: P(loss ≤ 0.0025)=0.0728; not established at 95% probability.

## Population Frequency

**Provisional decision:** sensitivity-only

**Provisional subset/representation:** maximum-population AF sensitivity

Max-population AF is strictly non-inferior to both and global AF adds no stable value; class-associated missingness prevents primary promotion.

Nested-selected pooled OOF metrics: MCC 0.774646, PR-AUC 0.945935, ROC-AUC 0.951511, Brier 0.083509.

mRMR mean ranking: Max-population AF 1.0; Global AF 2.0.

Outer-fold subset selections: Max-population AF 3/5; Global AF + Max-population AF 2/5.

Paired comparisons:

- `both_vs_max_population_af`: difference -0.000071, 95% interval [-0.000785, 0.000639], P(diff>0)=0.4234.
- `max_population_af_vs_global_af`: difference 0.012392, 95% interval [0.008540, 0.016323], P(diff>0)=1.0000.

Strict 0.0025-MCC non-inferiority sensitivity:

- Reduced `gnomad_max_population_af`: P(loss ≤ 0.0025)=1.0000; reduced non-inferior in sensitivity analysis.

This block remains sensitivity-only regardless of training performance.

## Computational Predictors

**Provisional decision:** sensitivity-only

**Provisional subset/representation:** REVEL + CADD + AlphaMissense non-inferior diagnostic subset

The trio is non-inferior at the strict margin, but nested subset identity is unstable; provenance and circularity governance prevents primary promotion.

Nested-selected pooled OOF metrics: MCC 0.876369, PR-AUC 0.979231, ROC-AUC 0.981487, Brier 0.047903.

mRMR mean ranking: REVEL 1.0; AlphaMissense 2.0; PolyPhen-2 3.4; CADD 3.8; SIFT 4.8.

Outer-fold subset selections: SIFT + PolyPhen-2 + REVEL + CADD + AlphaMissense 2/5; PolyPhen-2 + REVEL + CADD + AlphaMissense 2/5; SIFT + REVEL + CADD + AlphaMissense 1/5.

Paired comparisons:

- `all_five_vs_revel_cadd_alphamissense`: difference 0.000215, 95% interval [-0.001060, 0.001530], P(diff>0)=0.6184.

Strict 0.0025-MCC non-inferiority sensitivity:

- Reduced `predictor_alphamissense_score__predictor_cadd_phred__predictor_revel_score`: P(loss ≤ 0.0025)=0.9998; reduced non-inferior in sensitivity analysis.

This block remains sensitivity-only regardless of training performance.

## Overall provisional configuration

Because conservation and amino-acid removal decisions remain method-dependent, the defensible provisional primary configuration is to retain the complete current primary biological feature set until supervisor review. Treat Missense Z plus LOEUF and the retained conservation measures as grouped evidence rather than independent biological votes. The reduced phyloP + PhastCons alternative remains plausible at 0.005 and 0.010 MCC margins but is not confirmed at the strict 0.0025 margin.

Population frequency remains sensitivity-only, with maximum-population AF as the parsimonious sensitivity representation. Computational predictors remain sensitivity-only; REVEL + CADD + AlphaMissense is a non-inferior diagnostic subset, but nested subset identity is unstable and provenance/circularity governance is unchanged.

No authoritative configuration or validation command was modified.
