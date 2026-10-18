# Governed Rare-Missense-SNV Dataset Rebuild

This package reconstructs the research dataset as a sequence of small, inspectable stages. It does **not** depend on the earlier 40,000-row master file. That file should be retained only as a historical checkpoint and schema reference.

The workflow is deliberately split into numbered scripts so that you can explain:

1. which source entered at each stage;
2. which records were removed and why;
3. which key joined two sources;
4. which database release was used;
5. what was missing rather than silently converting missing data to zero; and
6. which fields are permitted in the core model.

## What the pipeline produces

- A versioned, balanced rare-missense cohort.
- A full provenance table for scientific review.
- A restricted core model table containing candidate biological evidence only.
- A separate optional table containing composite pathogenicity-predictor scores.
- Gene-separated train, validation and test partitions.
- JSON audit records and a final readable audit report.

It does **not** construct PS1/PM5-like label-derived features before splitting. If these are investigated later, they must be calculated from training labels only and treated as a separate experiment.

## 1. Set up the project

On macOS or Linux:

```bash
cd rare_missense_dataset_pipeline
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
python run_pipeline.py run 00
```

On Windows PowerShell, activation is:

```powershell
.venv\Scripts\Activate.ps1
```

The environment report is written to `audits/00_environment_report.json`.

## 2. Decide between a pilot and the full rebuild

Use `config/pilot_config.json` first. It selects 100 candidate variants per class and uses one VEP REST worker. This lets you learn the order and inspect each output before spending time on the full cohort.

Use `config/pipeline_config.json` for the governed full run. It selects an oversampled candidate pool of 50,000 benign-class and 25,000 pathogenic-class records before enrichment, then attempts to select 20,000 validated records per class at the end. These values are configuration choices, not hidden assumptions.

Do not treat pilot outputs as full-run outputs. Archive the successful pilot before starting the full rebuild. The `run-full` command uses a configuration-bound state file, so it does not mistake existing pilot files for completed full stages.

## 3. Freeze one ClinVar release

Preferred: download a specific archived monthly GRCh38 ClinVar VCF, then register it:

```bash
python run_pipeline.py run 01 -- \
  --local-vcf /absolute/path/to/clinvar_GRCh38_release.vcf.gz
```

If an archived file is not yet available, the script can download the live monthly file and immediately freeze its header date and SHA-256:

```bash
python run_pipeline.py run 01
```

Never overwrite the registered file during the same study. `--force` exists only for an intentional restart.

ClinVar publishes the GRCh38 short-variant VCF monthly and archives monthly releases. Record the exact file date, URL and SHA-256 in the study log.

## 4. Extract eligible ClinVar candidates

```bash
python run_pipeline.py run 02 -- --config config/pilot_config.json
python run_pipeline.py run 03 -- --config config/pilot_config.json
```

For the full run, replace `pilot_config.json` with `pipeline_config.json`.

Stage 02 retains:

- GRCh38 biallelic A/C/G/T substitutions;
- a ClinVar `missense_variant` molecular consequence;
- unambiguous benign/likely-benign or pathogenic/likely-pathogenic labels;
- at least the configured ClinVar review level; and
- one deterministic record per genomic variant key.

Stage 03 uses a fixed random seed to choose the pre-enrichment pool. Balancing is **not** finalised here because rarity, transcript and residue validation can remove different proportions from each class.

## 5. Annotate and reconcile variants with Ensembl VEP

Prepare VCF input:

```bash
python run_pipeline.py run 04
```

### Pilot REST run

```bash
python run_pipeline.py run 05 -- \
  --mode rest \
  --config config/pilot_config.json
```

REST results are saved in resumable batch checkpoints.

### Full checkpointed REST run

The governed full configuration uses Ensembl REST in batches of no more than 200 variants. Stage 05 first queries the official metadata endpoint and stops if the live Ensembl release differs from the expected release in the configuration. Each response is cached under a request-derived hash containing the variants, endpoint, options, assembly and verified release. Therefore, a pilot checkpoint cannot be silently reused for a different full-run batch.

Inspect the REST settings without making a request:

```bash
python run_pipeline.py run 05 -- \
  --mode rest \
  --config config/pipeline_config.json \
  --dry-run
```

Stage 05 is normally invoked by `run-full`. Local or externally generated VEP JSONL remains available as a contingency and can be registered explicitly:

```bash
python run_pipeline.py run 05 -- \
  --mode register \
  --existing-jsonl /absolute/path/to/vep_output.jsonl
```

Parse and reconcile the results:

```bash
python run_pipeline.py run 06 -- --config config/pilot_config.json
```

The transcript selection is deterministic:

1. require a coding missense consequence;
2. require the VEP gene to match a ClinVar gene;
3. prefer the primary ClinVar gene;
4. prefer MANE Select, then canonical; and
5. apply a deterministic final ordering if more than one candidate remains.

VEP population results have three different states:

- `FOUND`: a matching gnomAD frequency was returned;
- `NOT_REPORTED`: VEP succeeded but no matching frequency was returned; the value stays missing, never zero; and
- `LOOKUP_FAILED`: the annotation result itself was unavailable.

## 6. Validate the protein and add protein context

```bash
python run_pipeline.py run 07 -- --config config/pilot_config.json
python run_pipeline.py run 08 -- --config config/pilot_config.json
```

Stage 07 prefers a VEP-linked Swiss-Prot accession and checks that the amino acid at the selected protein position matches the VEP reference residue. Only `MATCH` records are eligible for the final cohort.

Stage 08 queries InterPro once per validated protein and records whether the variant residue overlaps a documented entry. The cached API results make the stage resumable.

## 7. Optionally add direct gnomAD counts

VEP supplies the population frequencies used by the default full workflow. Direct AC, AN and homozygote counts are an optional experimental enrichment because repeated remote indexed-VCF queries can be slow or unstable at full scale. They are disabled by default and are not part of the core model table.

To test the configured remote sources manually:

```bash
BASE="https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf"
python run_pipeline.py run 08b -- \
  --exomes-vcf "$BASE/exomes/gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz" \
  --genomes-vcf "$BASE/genomes/gnomad.genomes.v4.1.sites.{chrom}.vcf.bgz" \
  --source-release "gnomAD v4.1" \
  --smoke-test 10
```

Set `gnomad.direct_counts_enabled` to `true` only after a reviewed full-scale access plan. Absence from a sites VCF is stored as `NOT_REPORTED`, not AF=0. A lookup error has a separate state.

## 8. Merge gnomAD gene constraint

Download and preserve the official gene-constraint table. If direct counts were not used:

```bash
python run_pipeline.py run 09 -- \
  --constraint-file /absolute/path/to/gnomad_constraint.tsv \
  --source-release "gnomAD v4.x"
```

If Stage 08b was used, add:

```bash
--input data/interim/08b_gnomad_counts.csv
```

The script matches the reconciled Ensembl transcript first. It falls back to a gene only when the source table has exactly one row for that gene; ambiguous gene matches are retained as missing and audited.

## 9. Add dbNSFP conservation values through MyVariant.info

No local dbNSFP archive is required. First run a separate ten-row smoke test:

```bash
python run_pipeline.py run 10 -- \
  --config config/pilot_config.json \
  --smoke-test 10
```

Then run the complete stage through the chosen configuration:

```bash
python run_pipeline.py run 10 -- --config config/pipeline_config.json
```

Stage 10 records the MyVariant build and API-reported dbNSFP release, verifies every returned coordinate and allele, and uses content-addressed raw-response caches. GERP, PhyloP and PhastCons enter the core evidence table. SIFT, PolyPhen-2, REVEL, CADD and AlphaMissense remain in `10_optional_predictor_scores.csv` and are never silently inserted into the core feature set.

## 10. Calculate amino-acid descriptors

```bash
python run_pipeline.py run 11
```

This stage calculates Grantham distance and transparent composition, polarity, volume, hydropathy, charge and relative-position descriptors. They are deterministic descriptors, not proof of pathogenicity.

## 11. Build the final dataset

Pilot:

```bash
python run_pipeline.py run 12 -- \
  --config config/pilot_config.json \
  --allow-smaller
```

Full rebuild:

```bash
python run_pipeline.py run 12 -- --config config/pipeline_config.json
```

Final eligibility requires a strict gene-matched missense transcript, a matching UniProt reference residue and a resolved rarity state. The configuration currently permits a successfully queried but unreported gnomAD frequency while retaining it as missing; lookup failures are excluded. Change this policy only through a documented protocol amendment.

Outputs:

- `final_dataset_full.csv`: labels, evidence, identifiers and provenance;
- `model_table_core.csv`: only the default core candidate features plus key, split gene and label;
- `model_table_with_optional_predictors.csv`: separate experimental condition, if scores exist; and
- `evidence_feature_registry.csv`: the explicit modelling whitelist.

## 12. Create gene-separated partitions

```bash
python run_pipeline.py run 13 -- --config config/pilot_config.json
```

Use the full configuration for the full run. All variants from one reconciled gene remain in one partition. Fit every learned transformation on the training genes only. Use validation genes for choices and tuning. Use test genes once after the complete method is locked.

## 13. Run the final audit

```bash
python run_pipeline.py run 14
python run_pipeline.py status
```

Read both:

- `audits/14_final_dataset_audit.json`; and
- `audits/14_final_dataset_audit.md`.

A structural `PASS` is not permission to ignore warnings. Large label-specific missingness gaps require controlled feature/missingness experiments before model development.

## 14. Run the code tests

```bash
python -m unittest discover -s tests -v
```

## 15. Run the governed full rebuild automatically

The full run requires the frozen ClinVar file and its provenance JSON, plus the configured gnomAD constraint table. First print and validate the complete plan:

```bash
python run_pipeline.py run-full \
  --config config/pipeline_config.json \
  --dry-run
```

Start the full rebuild:

```bash
python run_pipeline.py run-full \
  --config config/pipeline_config.json
```

If the network, terminal or computer interrupts the run, resume it with:

```bash
python run_pipeline.py run-full \
  --config config/pipeline_config.json \
  --resume
```

Completed stages are skipped only when the state file was created from the same configuration and the recorded output SHA-256 still matches. API stages also retain their own content-addressed caches. The command constructs and audits the dataset; it does not train a model.

## Files that should not be committed to a public Git repository

- full ClinVar, gnomAD or dbNSFP source files;
- API caches containing entire returned records;
- generated full/intermediate datasets; and
- any source with redistribution restrictions.

Commit the scripts, configuration, source URLs, release identifiers, checksums, audit reports and a small non-sensitive test fixture instead.

## Official source pages

- ClinVar access and releases: <https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/>
- Ensembl VEP: <https://www.ensembl.org/info/docs/tools/vep/index.html>
- UniProt REST: <https://rest.uniprot.org/>
- InterPro API: <https://www.ebi.ac.uk/interpro/api/>
- gnomAD downloads: <https://gnomad.broadinstitute.org/downloads>
- MyVariant.info API: <https://docs.myvariant.info/>
- dbNSFP: <https://www.dbnsfp.org/>

