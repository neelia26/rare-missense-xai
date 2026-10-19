# Dataset construction

## 1. Data source and provenance

The dataset was constructed as a governed, reproducible collection of rare missense single-nucleotide variants on the GRCh38 assembly. The primary labels and genomic variants were obtained from the ClinVar GRCh38 VCF frozen on 11 October 2026. Each subsequent source was recorded through a stage manifest, an audit report and file checksums. Ensembl VEP release 116, UniProt release 2026_03, InterPro release 110.0, gnomAD v4.1 and dbNSFP 4.8a accessed through MyVariant.info build 20250624 supplied the principal annotation evidence. The final source and processing records are summarised in `dataset_sources.csv`. [Citation to be linked from proposal reference list]

## 2. ClinVar label definition and review-status filtering

ClinVar records were restricted to biallelic GRCh38 single-nucleotide substitutions annotated as missense variants. Benign, likely benign and combined benign/likely benign classifications formed class 0. Pathogenic, likely pathogenic and combined pathogenic/likely pathogenic classifications formed class 1. Uncertain, conflicting and other non-binary classifications were excluded. A minimum ClinVar review status of one star was required, and one deterministic record was retained per genomic variant key. The filtering stage scanned 4,594,764 source VCF records and retained 209,183 eligible labelled variants: 147,254 class-0 and 61,929 class-1 variants. [Citation to be linked from proposal reference list]

## 3. Candidate selection and class balancing

The eligible ClinVar set was deterministically sampled with random seed 20261018 to form an oversized candidate pool of 75,000 variants. This pool contained 50,000 class-0 and 25,000 class-1 variants. Oversampling before biological validation allowed the final balanced cohort to be selected only after rarity, transcript and protein-quality controls had been applied. The governed final target was 20,000 variants per class.

## 4. VEP consequence and transcript reconciliation

The candidate pool was submitted in content-addressed batches to the Ensembl REST Variant Effect Predictor. Release 116 was verified both before and after annotation. Transcript consequences were reconciled deterministically, prioritising an appropriate gene-matched MANE Select transcript, followed by documented canonical and transcript-selection rules. A row remained eligible only when the selected transcript represented a strict missense consequence with a single reconciled amino-acid substitution and protein position. The reconciliation process retained explicit transcript-selection provenance rather than silently substituting an alternative consequence. [Citation to be linked from proposal reference list]

## 5. UniProt protein validation

VEP Swiss-Prot accessions were preferred for protein validation. The reconciled reference amino acid and residue position were checked against the selected UniProt sequence from release 2026_03. Only rows with a matching reference residue remained eligible for the final frozen cohort. Missing accessions, API errors, out-of-range positions and reference-residue mismatches were recorded as distinct outcomes. UniProt validation therefore acted as a biological consistency control rather than as a source of labels. [Citation to be linked from proposal reference list]

## 6. InterPro and UniProt protein-context enrichment

Validated protein positions were enriched with InterPro release 110.0 and UniProt release 2026_03 context. InterPro overlap types were normalised into token-level binary indicators, distinguishing a successful lookup with no overlap from a lookup failure. UniProt context was likewise represented by token-level indicators rather than by complete semicolon-separated strings. Four curated tokens—Natural variant, Mutagenesis, Sequence conflict and Alternative sequence—were prohibited from primary modelling because they could encode prior variant or experimental knowledge. The remaining structural and positional context features were retained subject to evidence-group ablation. [Citation to be linked from proposal reference list]

## 7. Population-frequency and gene-constraint evidence

Population-frequency fields were obtained from the coordinate- and allele-reconciled annotation record and governed by a maximum allele-frequency threshold of 0.01. A successfully completed lookup with no reported frequency remained missing and was never converted to zero. Unresolved population lookup failures were excluded under the governed pipeline rule. Gene-level missense constraint and LOEUF values were merged from gnomAD v4.1 using documented gene/transcript matching. Population evidence was retained only as a sensitivity condition because its missingness differed strongly by class; the no-population table was designated the primary leakage-controlled baseline. Gene constraint was retained with a mandatory no-constraint ablation because it could act as a gene-level proxy. [Citation to be linked from proposal reference list]

## 8. Conservation and optional computational predictors

Coordinate- and allele-verified dbNSFP 4.8a records were accessed through MyVariant.info build 20250624. GERP++, phyloP100way and phastCons100way values formed the conservation evidence group. When multiple valid values were returned, the documented reduction rule selected the maximum valid numeric value. Optional SIFT, PolyPhen-2 HDIV, REVEL, CADD PHRED and AlphaMissense scores were separated from primary evidence. They were not approved for primary modelling because their training-data and methodological overlap required a separate literature and circularity review. [Citation to be linked from proposal reference list]

## 9. Amino-acid-change features

Physicochemical features were calculated deterministically from the reconciled reference and alternate amino acids. These included Grantham distance, composition, polarity, volume, hydropathy and charge differences, charge-class and polarity-class change indicators, and relative residue position within the validated protein. Relative protein position was registered as numeric protein-context evidence; the two class-change fields were registered as binary amino-acid-change evidence.

## 10. Missing-data handling

Missing values were preserved as missing rather than reinterpreted as biological absence. The governed missing representations were blank, `.`, `NA`, `NaN`, `None` and `null`. A successful InterPro lookup with no overlap was encoded as zero for each token, whereas a lookup failure remained missing. Population-frequency availability was not supplied to a model as an indicator. Raw population NaNs were not considered safe for models capable of learning missing-value routing. Any later imputation, scaling or encoding must be fitted on training genes only, and unavailable population evidence must contribute neutrally and be described as unavailable in explanations.

## 11. Quality-control and exclusion rules

The 75,000-row candidate pool was assessed for binary label validity, strict transcript-level missense eligibility, UniProt reference-residue agreement, a valid grouping gene and governed rarity status. The mutually exclusive recorded exclusions were 1,195 transcript-ineligible rows, 13,076 rows without a validated UniProt reference residue and 7,734 variants that were not rare. This left 52,995 eligible variants before balancing: 30,606 class-0 and 22,389 class-1 variants. Identifiers, variant keys and gene symbols were retained as metadata only and were not approved predictors; the binary label was the target only.

## 12. Frozen 40,000-variant dataset

Deterministic class balancing selected 20,000 eligible variants from each class, producing 40,000 unique variants with no blank or duplicated variant keys. The frozen file contained 104 columns and had SHA-256 `50d8ac02ededae696d6867c7e713351d33058d01bad4fab3633055157525494f`. The dataset-build audit passed without critical failures. The frozen file was preserved unchanged as the source for all modelling-readiness derivations.

## 13. Derivation of the 39,552-variant modelling cohort

The modelling cohort excluded 448 variants for which the ClinVar gene and VEP-reconciled gene differed. These comprised 358 class-0 and 90 class-1 variants and were retained in a separate exclusion log for sensitivity analysis. The resulting cohort contained 39,552 unique variants: 19,642 class-0 and 19,910 class-1 variants. Protein-context strings were normalised, curated UniProt leakage-risk tokens were excluded, and a no-population primary baseline was created without modifying the frozen source.

## 14. Gene-separated train, validation and test partitions

The modelling cohort was partitioned by reconciled gene so that no gene or variant appeared in more than one partition. The training set contained 27,727 variants across 5,206 genes, with 13,741 class-0 and 13,986 class-1 variants. Validation contained 5,790 variants across 1,115 genes, evenly divided between the two classes. Test contained 6,035 variants across 1,115 genes, with 3,006 class-0 and 3,029 class-1 variants. The test partition was locked for one final evaluation. Test labels were not used for missingness, feature-association, feature-selection or methodological decisions during the static circularity audit.

## 15. Reproducibility and governance

Every pipeline stage produced a manifest containing input and output checksums, parameters, software context and completion evidence. The frozen dataset was tagged as `dataset-v1.0-40000`; the modelling cohort, feature registry, alternate-gene exclusions and locked gene assignments were governed separately. Predictor selection must follow the feature and evidence-lineage registries. All preprocessing must be learned from training genes only. The validation partition may confirm decisions established from training data, while the test partition remains reserved for one final locked evaluation. Static correlation or association flags identify redundancy candidates but do not justify automatic feature removal. Scientific clearance remains conditional on population handling, optional-predictor literature review and the prescribed diagnostic ablations.
