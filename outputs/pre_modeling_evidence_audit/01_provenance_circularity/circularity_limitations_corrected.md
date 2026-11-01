# Corrected circularity limitations

GERP++, phyloP, PhastCons and SIFT do not have conventional clinically labelled training-variant overlap; those checks are `NOT_APPLICABLE`. Their material concern is overlapping evolutionary signal.

PolyPhen-2 HumDiv was protein-matched using accession, amino-acid position, reference and alternate residues. Status: **TESTED**; exact matches: **940**; protein overlap: **819**; mapped-gene overlap: **820**; training records not mapped to a current cohort protein: **1316**. No genomic coordinate conversion was attempted.

REVEL remains `UNKNOWN/UNVERIFIABLE` for exact pathogenic overlap because the HGMD-derived portion is restricted. CADD uses proxy labels, not clinical labels, but its exact release is `UNKNOWN EXACT RELEASE`. AlphaMissense used population-frequency information for fine-tuning and ClinVar for evaluation; correlation cannot determine exact training overlap. Unavailable data never demonstrate absence of overlap.
