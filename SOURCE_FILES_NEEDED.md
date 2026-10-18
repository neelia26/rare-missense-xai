# Source Files and Existing Files That Would Help

## Required for a clean full rebuild

1. **ClinVar GRCh38 VCF (`.vcf.gz`)**
   - Prefer an archived monthly release.
   - If you do not already have one, Stage 01 can download and freeze the current monthly file.

2. **Ensembl VEP release and matching GRCh38 cache**
   - Needed for the full local run.
   - A small pilot can use the REST route instead.

3. **gnomAD gene-constraint table**
   - Must contain a gene, Ensembl transcript, missense constraint and LOEUF-equivalent column.
   - Record its version and checksum.

4. **dbNSFP variant table plus tabix index**
   - Must contain GRCh38 coordinate/ref/alt fields and the selected conservation columns.
   - Record the exact dbNSFP release and redistribution terms.

## Optional but useful

5. **Indexed gnomAD exome/genome site VCFs**
   - Used only if direct AC, AN and homozygote counts are wanted.
   - They are large; frequencies from pinned VEP can be used without them.

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

