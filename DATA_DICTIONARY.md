# Output and Evidence Dictionary

## Stable identifiers and targets

| Field | Purpose | Model input? |
|---|---|---:|
| `variant_key` | `GRCh38:chromosome:position:ref:alt` join key | No |
| `binary_label` | ClinVar-derived target: 0 benign-class, 1 pathogenic-class | Target only |
| `reconciled_gene_symbol` | Gene-group partitioning and provenance | No |
| `reconciled_transcript_id` | Transcript provenance | No |

## Candidate evidence groups

| Group | Main fields | Interpretation status |
|---|---|---|
| Population | gnomAD AF; optional AC/AN/homozygote counts | Candidate; missing/not-reported stays explicit |
| Conservation | GERP, PhyloP, PhastCons | Candidate; audit overlap before modelling |
| Gene constraint | missense Z and LOEUF | Candidate; audit gene-identity shortcuts |
| Amino-acid change | Grantham, polarity, volume, hydropathy, charge | Deterministic candidate descriptors |
| Protein context | UniProt feature and InterPro overlap | Candidate; absence differs from failed annotation |

The exact model inputs are defined by `data/processed/evidence_feature_registry.csv`. Do not select columns by taking “everything numeric.”

## Missingness states

| State | Meaning |
|---|---|
| `FOUND` | Source returned a usable value or record |
| `NOT_REPORTED` | Lookup succeeded, but the requested value/allele was not reported |
| `NOT_FOUND` | A source entity such as an accession was not found |
| `LOOKUP_FAILED` / `ERROR` | Technical retrieval or parsing failure |
| blank numeric value | No numeric measurement is available; never silently converted to zero |

## Tables

| File | Purpose |
|---|---|
| `final_dataset_full.csv` | Complete governed record for scientific inspection |
| `model_table_core.csv` | Whitelisted core candidate evidence condition |
| `model_table_with_optional_predictors.csv` | Separate composite-score experiment |
| `split_gene_assignments.csv` | Gene-to-partition mapping |
| `splits/train_core.csv` | Training data; fit preprocessing here only |
| `splits/validation_core.csv` | Model/configuration choice and explanation checks |
| `splits/test_core.csv` | One-time locked final evaluation |

