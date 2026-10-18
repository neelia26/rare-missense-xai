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

Do not mix pilot and full outputs in one directory. After completing the pilot, archive or remove the pilot `data/`, `audits/` and `logs/` outputs before starting the full run.

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

### Full local VEP run

Install the Ensembl VEP release and matching GRCh38 cache specified in the configuration. Confirm the command before execution:

```bash
python run_pipeline.py run 05 -- \
  --mode local \
  --config config/pipeline_config.json \
  --dry-run
```

Then run the same command without `--dry-run`. The command-line tool is the recommended route for the full cohort because it is intended for large-volume annotation. If VEP is run in Docker or on another computer, copy its JSON-lines output back and register it:

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

VEP supplies frequencies but not a complete set of AC, AN and homozygote counts. If you have local, indexed gnomAD site VCFs, run:

```bash
python run_pipeline.py run 08b -- \
  --exomes-vcf /absolute/path/to/gnomad.exomes.sites.vcf.bgz \
  --genomes-vcf /absolute/path/to/gnomad.genomes.sites.vcf.bgz \
  --source-release "gnomAD v4.x"
```

This stage is optional. Absence from a sites VCF is stored as `NOT_REPORTED`, not AF=0. A lookup error has a separate state.

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

## 9. Add dbNSFP conservation values

Prepare a licensed/local dbNSFP variant file as bgzip plus tabix. Then run:

```bash
python run_pipeline.py run 10 -- \
  --dbnsfp-file /absolute/path/to/dbNSFP_variant.tsv.gz \
  --source-release "dbNSFP x.y"
```

The core output receives GERP, PhyloP and PhastCons values. SIFT, PolyPhen-2, REVEL, CADD and AlphaMissense values, when available, go to `10_optional_predictor_scores.csv`. They are never silently inserted into the core feature set.

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
- dbNSFP: <https://www.dbnsfp.org/>

