# Source Files and Existing Files That Would Help

## Required for a clean full rebuild

1. **ClinVar GRCh38 VCF (`.vcf.gz`)**
   - Prefer an archived monthly release.
   - If you do not already have one, Stage 01 can download and freeze the current monthly file.

2. **Ensembl VEP REST access**
   - The full workflow uses checkpointed REST batches and verifies the live Ensembl release against the configured expected release.
   - A local release-matched VEP cache is optional as a contingency, not a requirement.

3. **gnomAD gene-constraint table**
   - Must contain a gene, Ensembl transcript, missense constraint and LOEUF-equivalent column.
   - Record its version and checksum.

4. **MyVariant.info API access for dbNSFP annotations**
   - Stage 10 retrieves only the required GRCh38 records, so no local dbNSFP archive or tabix index is required.
   - The pipeline verifies returned coordinates and alleles, records the MyVariant build and API-reported dbNSFP release, and preserves content-addressed response caches.
   - Review the dbNSFP licence before redistributing any derived dataset.

## Optional but useful

5. **Indexed gnomAD exome/genome site VCFs**
   - Used only if direct AC, AN and homozygote counts are enabled as a separate enrichment.
   - The default full workflow leaves this disabled and uses release-verified VEP population frequencies.

6. **Existing UniProt or InterPro caches**
   - Reusable only after checking that the accession, retrieval date/release and response format are recorded.

## Historical files that could save time or improve provenance

If you have any of these, provide them before rerunning the expensive stages:

- the earlier 62,000-row pre-VEP candidate file;
- complete saved VEP JSON/checkpoint batches;
- the 38,662-row rare, transcript-reconciled checkpoint;
- the 35,188-row UniProt-validated checkpoint;
- the UniProt protein cache;
- the InterPro cache;
- the original gnomAD constraint file;
- the original dbNSFP/MyVariant response cache; or
- the exact 30,523-row new cohort before it was merged with the earlier 9,477 rows.

These files will not be trusted automatically. They can be registered only after their schema, release, checksum, completeness and join keys are audited. The old 40,000-row master is not required to rebuild the governed dataset.

