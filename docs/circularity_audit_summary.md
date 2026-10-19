# Static circularity-audit summary

## Decision

- **Structural status:** PASS
- **Scientific clearance:** CONDITIONAL
- **Critical failures:** none
- **Label-dependent analysis partition:** training only
- **Test use:** structural integrity checks only; test labels were not read for feature association or methodological decisions

The frozen 40,000-variant source checksum, the 39,552-row modelling cohort, locked split sizes, schemas, unique keys and zero gene/variant overlap were verified. The audit did not train a model and did not alter any dataset or split assignment.

## Class-associated missingness

Two fields exceeded the pre-specified absolute training-class gap of 0.10:

| Field | Training missing | Benign missing | Pathogenic missing | Absolute gap |
|---|---:|---:|---:|---:|
| `gnomad_global_af` | 40.25% | 11.88% | 68.13% | 56.25% |
| `gnomad_max_population_af` | 40.25% | 11.88% | 68.13% | 56.25% |

This confirms that population availability could act as a label shortcut. The no-population table therefore remains the primary leakage-controlled condition. No missingness/availability indicator is permitted, population-inclusive modelling remains a sensitivity experiment, and later preprocessing must be fitted on training genes only.

## Redundancy candidates

Two training-only numeric pairs met the absolute Pearson or Spearman threshold of 0.85:

| Field A | Field B | Paired observations | Pearson | Spearman | Trigger |
|---|---|---:|---:|---:|---|
| `aa_hydropathy_delta` | `aa_polarity_delta` | 27,727 | -0.8651 | -0.8287 | Pearson |
| `gnomad_global_af` | `gnomad_max_population_af` | 16,567 | 0.6988 | 0.8764 | Spearman |

These are candidates for grouped interpretation and ablation, not automatic deletion. No binary/categorical pair reached corrected Cramér's V of 0.80.

## Governance outcome

The evidence-lineage registry records 60 governed fields. The decision output contains 49 approved-primary features, seven sensitivity-only fields and four prohibited UniProt fields. The sensitivity-only set comprises the two population-frequency fields and five optional computational predictors. The prohibited set comprises Natural variant, Mutagenesis, Sequence conflict and Alternative sequence indicators.

The following optional predictor questions remain unresolved pending literature review:

- SIFT training-data and methodological overlap;
- PolyPhen-2 HDIV training-data and methodological overlap;
- REVEL training-data and methodological overlap;
- CADD training-data and methodological overlap;
- AlphaMissense training-data and methodological overlap.

No overlap claim was inferred where the repository did not provide sufficient evidence.

## Mandatory later controls

The modelling study must include:

1. no-population versus population-inclusive conditions;
2. an ablation without gene-constraint evidence;
3. evidence-group removal experiments;
4. a missingness-only diagnostic model;
5. optional computational predictors only after their circularity review.

Scientific clearance remains conditional until these controls, neutral population handling and the predictor literature review are completed. The locked test partition must be used once for final evaluation only.
