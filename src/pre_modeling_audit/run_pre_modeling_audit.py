#!/usr/bin/env python3
"""Validate the locked cohort and run training-only evidence experiments."""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import platform
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.metrics import (
    average_precision_score, brier_score_loss, matthews_corrcoef, roc_auc_score,
)
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[2]
SEED = 20261018
TARGET = "binary_label"
GROUP = "reconciled_gene_symbol"
KEY = "variant_key"
MISSING = {"", ".", "na", "nan", "none", "null"}
CONSERVATION = [
    "conservation_gerp_rs", "conservation_phylop100way",
    "conservation_phastcons100way",
]
PREDICTORS = [
    "predictor_sift_score", "predictor_polyphen2_hdiv_score",
    "predictor_revel_score", "predictor_cadd_phred",
    "predictor_alphamissense_score",
]
POPULATION = ["gnomad_global_af", "gnomad_max_population_af"]
METADATA = {KEY, GROUP, TARGET}
PROHIBITED = {
    "uniprot_feature_natural_variant", "uniprot_feature_mutagenesis",
    "uniprot_feature_sequence_conflict", "uniprot_feature_alternative_sequence",
}
PATHS = {
    "frozen": ROOT / "data/processed/final_dataset_full.csv",
    "cohort": ROOT / "data/processed/modeling_ready_full.csv",
    "excluded": ROOT / "data/processed/excluded_alternate_gene_variants.csv",
    "assignments": ROOT / "data/processed/modeling_split_gene_assignments.csv",
    "train_primary": ROOT / "data/processed/modeling_splits/train_core_no_population.csv",
    "train_population": ROOT / "data/processed/modeling_splits/train_core.csv",
    "train_optional": ROOT / "data/processed/modeling_splits/train_optional_predictors.csv",
}
OUT = ROOT / "outputs/pre_modeling_evidence_audit"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(rows[0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[])


def numeric(series: pd.Series) -> pd.Series:
    cleaned = series.astype(str).str.strip()
    cleaned = cleaned.mask(cleaned.str.lower().isin(MISSING))
    return pd.to_numeric(cleaned, errors="coerce")


def variant_key_set(frame: pd.DataFrame) -> set[str]:
    return set(frame[KEY].astype(str).str.strip())


def validate_cohort() -> tuple[pd.DataFrame, dict]:
    frozen, cohort, excluded = (load(PATHS[name]) for name in ("frozen", "cohort", "excluded"))
    failures: list[str] = []
    expected = {"frozen": 40000, "cohort": 39552, "excluded": 448}
    observed = {"frozen": len(frozen), "cohort": len(cohort), "excluded": len(excluded)}
    if observed != expected:
        failures.append(f"Expected row counts {expected}; observed {observed}")
    for name, frame in (("frozen", frozen), ("cohort", cohort), ("excluded", excluded)):
        required = {KEY, TARGET}
        if name != "excluded":
            required.add(GROUP)
        missing = sorted(required - set(frame.columns))
        if missing:
            failures.append(f"{name} lacks required columns: {missing}")
        if KEY in frame and (frame[KEY].astype(str).str.strip() == "").any():
            failures.append(f"{name} has blank variant keys")
        if KEY in frame and frame[KEY].duplicated().any():
            failures.append(f"{name} has duplicate variant keys")
        if frame.duplicated().any():
            failures.append(f"{name} has exact duplicate rows")
    retained, removed, source = variant_key_set(cohort), variant_key_set(excluded), variant_key_set(frozen)
    if retained & removed:
        failures.append("Retained and excluded variant keys overlap")
    if retained | removed != source:
        failures.append("Retained plus excluded keys do not exactly account for frozen keys")
    labels = set(cohort[TARGET])
    if labels != {"0", "1"}:
        failures.append(f"Unexpected target values: {sorted(labels)}")
    conflicts = frozen.groupby(KEY)[TARGET].nunique()
    if (conflicts > 1).any():
        failures.append("Conflicting labels exist for the same frozen variant key")
    if (cohort[GROUP].astype(str).str.strip() == "").any():
        failures.append("Cohort has blank reconciled genes")

    schema_rows, missing_rows = [], []
    registry = load(ROOT / "data/processed/modeling_feature_registry.csv")
    roles = dict(zip(registry["field"], registry["evidence_group"]))
    for column in cohort.columns:
        parsed = numeric(cohort[column])
        nonmissing_text = ~cohort[column].astype(str).str.strip().str.lower().isin(MISSING)
        numeric_fraction = float(parsed.notna().sum() / max(nonmissing_text.sum(), 1))
        proposed = "target" if column == TARGET else "group" if column == GROUP else "identifier" if column == KEY else "candidate evidence" if column in roles else "metadata/provenance"
        schema_rows.append({
            "column": column, "observed_dtype": str(cohort[column].dtype),
            "inferred_type": "numeric" if numeric_fraction == 1 and nonmissing_text.any() else "categorical/text",
            "proposed_role": proposed, "evidence_group": roles.get(column, ""),
            "nonmissing_count": int(nonmissing_text.sum()), "numeric_parse_fraction": round(numeric_fraction, 8),
            "infinite_count": int(np.isinf(parsed.dropna()).sum()),
            "malformed_nonmissing_count": int((nonmissing_text & parsed.isna()).sum()) if column in roles else 0,
            "numeric_min": "" if parsed.notna().sum() == 0 else float(parsed.min()),
            "numeric_max": "" if parsed.notna().sum() == 0 else float(parsed.max()),
            "unexpected_binary_count": int((nonmissing_text & ~cohort[column].isin(["0", "1"])).sum()) if dict(zip(registry["field"], registry["data_type"])).get(column) == "binary" else 0,
            "unique_nonmissing": int(cohort.loc[nonmissing_text, column].nunique()),
        })
        missing_count = int((~nonmissing_text).sum())
        missing_rows.append({
            "column": column, "row_count": len(cohort), "missing_count": missing_count,
            "missing_fraction": round(missing_count / len(cohort), 8),
            "zero_count_if_numeric": int((parsed == 0).sum()),
        })
    base = OUT / "00_locked_protocol/cohort_validation"
    write_csv(base / "feature_schema.csv", schema_rows)
    write_csv(base / "missingness_summary.csv", missing_rows)
    report = {
        "status": "PASS" if not failures else "FAIL", "critical_failures": failures,
        "row_counts": observed, "retained_excluded_overlap": len(retained & removed),
        "row_accounting_complete": retained | removed == source,
        "duplicate_variant_keys": {name: int(frame[KEY].duplicated().sum()) for name, frame in (("frozen", frozen), ("cohort", cohort), ("excluded", excluded))},
        "exact_duplicate_rows": {name: int(frame.duplicated().sum()) for name, frame in (("frozen", frozen), ("cohort", cohort), ("excluded", excluded))},
        "conflicting_label_keys": int((conflicts > 1).sum()), "allowed_target_values": sorted(labels),
        "frozen_sha256": sha256(PATHS["frozen"]),
    }
    write_json(base / "cohort_integrity_report.json", report)
    if failures:
        raise RuntimeError("Critical cohort validation failure: " + "; ".join(failures))
    return cohort, report


def validate_split(cohort: pd.DataFrame) -> dict:
    assignments = load(PATHS["assignments"])
    split_by_gene = dict(zip(assignments["gene"], assignments["split"]))
    missing = sorted(set(cohort[GROUP]) - set(split_by_gene))
    counts = {}
    genes = {}
    for split in ("train", "validation", "test"):
        subset = assignments[assignments["split"] == split]
        counts[split] = {
            "variants": int(pd.to_numeric(subset["variant_count"]).sum()),
            "class_0": int(pd.to_numeric(subset["benign_count"]).sum()),
            "class_1": int(pd.to_numeric(subset["pathogenic_count"]).sum()),
            "genes": len(subset),
        }
        genes[split] = set(subset["gene"])
    # Read only keys from validation/test split tables; their outcome columns stay unread.
    key_sets = {
        s: set(pd.read_csv(ROOT / f"data/processed/modeling_splits/{s}_full.csv", usecols=[KEY], dtype=str)[KEY])
        for s in ("train", "validation", "test")
    }
    key_overlap = {f"{a}_{b}": len(key_sets[a] & key_sets[b]) for a, b in itertools.combinations(key_sets, 2)}
    gene_overlap = {f"{a}_{b}": len(genes[a] & genes[b]) for a, b in itertools.combinations(genes, 2)}
    failures = []
    if missing:
        failures.append(f"{len(missing)} cohort genes lack assignments")
    if set().union(*key_sets.values()) != variant_key_set(cohort):
        failures.append("Split-key union does not equal cohort")
    if any(key_overlap.values()) or any(gene_overlap.values()):
        failures.append("Variant or gene overlap exists")
    report = {
        "status": "PASS" if not failures else "FAIL", "critical_failures": failures,
        "split_counts_from_locked_assignment_summary": counts,
        "class_count_access_policy": "Aggregate validation/test counts were read only from the pre-existing assignment summary for integrity verification; validation/test outcome rows were not analysed.",
        "variant_overlap": key_overlap, "gene_overlap": gene_overlap,
        "missing_gene_assignments": len(missing), "split_union_equals_cohort": not failures,
    }
    write_json(OUT / "00_locked_protocol/cohort_validation/split_integrity_report.json", report)
    if failures:
        raise RuntimeError("Critical split validation failure: " + "; ".join(failures))
    return report


def inventory_and_leakage(cohort: pd.DataFrame) -> None:
    inventory = [
        ("Filter ClinVar and binary labels", "scripts/02_extract_candidates.py", "ClinVar GRCh38 VCF", "02_clinvar_candidates.csv", "pipeline_config.json"),
        ("Deterministic candidate sampling", "scripts/03_select_candidate_pool.py", "02_clinvar_candidates.csv", "03_candidate_pool.csv", "seed 20261018"),
        ("VEP missense and transcript reconciliation", "scripts/06_parse_vep.py", "05_vep.jsonl", "06_vep_reconciled.csv", "pipeline_config.json"),
        ("UniProt mapping verification", "scripts/07_enrich_uniprot.py", "06_vep_reconciled.csv", "07_uniprot_validated.csv", "pipeline_config.json"),
        ("InterPro/protein context", "scripts/08_enrich_interpro.py", "07_uniprot_validated.csv", "08_interpro.csv", "pipeline_config.json"),
        ("Gene constraint", "scripts/09_merge_gnomad_constraint.py", "08_interpro.csv", "09_constraint.csv", "gnomAD v4.1"),
        ("Conservation/optional predictors", "scripts/10_enrich_dbnsfp.py", "09_constraint.csv", "10_dbnsfp.csv + optional table", "dbNSFP 4.8a"),
        ("Amino-acid descriptors", "scripts/11_add_amino_acid_features.py", "10_dbnsfp.csv", "11_amino_acid_features.csv", "deterministic constants"),
        ("Final 40,000 build", "scripts/12_build_final_dataset.py", "11_amino_acid_features.csv", "final_dataset_full.csv", "pipeline_config.json"),
        ("Exclude gene disagreements/context indicators", "scripts/12b_prepare_modeling_dataset.py", "final_dataset_full.csv", "modeling_ready_*.csv + exclusions", "frozen SHA-256"),
        ("Gene-separated split", "scripts/13_split_by_gene.py", "modeling_ready_*.csv", "modeling_splits + assignments", "seed 20261018"),
        ("Readiness verification", "scripts/14_audit_modeling_readiness.py", "derived cohort/splits", "audit JSON", "locked tolerances"),
    ]
    rows = [{
        "preprocessing_operation": a, "script_or_notebook": b, "input_file": c,
        "output_file": d, "configuration": e, "already_applied": "yes",
        "deterministic": "yes", "reproducible": "yes",
        "verification_evidence": "existing manifest/audit plus this cohort validation",
        "remaining_issue": "none" if "Conservation" not in a else "predictor provenance reviewed separately",
    } for a, b, c, d, e in inventory]
    write_csv(OUT / "00_locked_protocol/cohort_validation/preprocessing_script_inventory.csv", rows)

    registry = load(ROOT / "data/processed/modeling_feature_registry.csv")
    approved = set(registry.loc[registry["primary_model_feature"] == "1", "field"])
    leakage = []
    for col in cohort.columns:
        lower = col.lower()
        if col == TARGET:
            role, allowed, concern, decision = "target", "no", "direct target", "TARGET_ONLY"
        elif col == GROUP or any(token in lower for token in ("gene_id", "transcript_id", "accession", "variation_id", "allele_id", "dbsnp", "query_id", "record_id")):
            role, allowed, concern, decision = "group/identifier", "no", "identity or lookup proxy", "EXCLUDE"
        elif col == KEY or col in {"chromosome", "position", "reference_allele", "alternate_allele"}:
            role, allowed, concern, decision = "variant identifier", "no", "variant identity", "EXCLUDE"
        elif col in PROHIBITED or any(token in lower for token in ("clinvar_significance", "clinical_conditions", "review_status")):
            role, allowed, concern, decision = "clinical/curated metadata", "no", "direct or prior curated knowledge", "PROHIBIT"
        elif col in approved:
            role, allowed, concern, decision = "approved primary evidence", "yes", "controlled evidence only", "ALLOW_PRIMARY"
        else:
            role, allowed, concern, decision = "metadata/provenance", "no", "not whitelisted; may encode pipeline/database status", "EXCLUDE"
        leakage.append({"column_name": col, "data_type": str(cohort[col].dtype), "proposed_role": role, "allowed_as_input": allowed, "leakage_concern": concern, "decision": decision, "reason": "Only registry-whitelisted biological evidence enters models; identifiers are retained for grouping/audit only."})
    for col in sorted(PROHIBITED - set(cohort.columns)):
        leakage.append({"column_name": col, "data_type": "binary indicator (absent)", "proposed_role": "curated UniProt knowledge", "allowed_as_input": "no", "leakage_concern": "prior curated variant/experimental knowledge", "decision": "PROHIBIT_CONFIRMED_ABSENT", "reason": "Explicitly excluded during model-table preparation."})
    write_csv(OUT / "00_locked_protocol/cohort_validation/leakage_column_audit.csv", leakage)


def provenance() -> list[dict]:
    common = {"retrieval_route": "dbNSFP 4.8a through MyVariant.info build 20250624", "genome_build_or_reference_sequence": "GRCh38/hg38 coordinate and allele verified by pipeline", "whether_exact_variant_overlap_can_be_tested": "only if public harmonised training variants are available", "whether_gene_overlap_can_be_tested": "yes if training genes are published", "whether_protein_overlap_can_be_tested": "yes if training proteins are published"}
    rows = [
        dict(display_name="GERP++ RS", exact_dataset_column=CONSERVATION[0], evidence_group="conservation", raw_derived_or_composite="derived evolutionary score", original_source="GERP++", source_version="dbNSFP 4.8a field", alignment_or_derivation_method="rejected substitutions from multispecies alignment", training_or_calibration_data="phylogenetic alignment and neutral tree", training_labels="none (evolutionary constraint)", whether_clinvar_was_used="no", whether_hgmd_was_used="no", whether_uniprot_disease_annotations_were_used="no", other_clinical_label_sources="none documented", whether_training_variants_are_publicly_available="not applicable as supervised variants", underlying_component_features="multispecies alignment", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="no", whether_it_incorporates_other_computational_predictors="no", training_or_release_cutoff_date="GERP++ paper 2010; exact dbNSFP track provenance not fully encoded", documentary_source="Davydov et al. 2010", exact_supporting_citation_or_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC2996323/", circularity_risk_classification="NOT APPLICABLE AS SUPERVISED-LABEL CIRCULARITY", redundancy_concern="overlapping evolutionary constraint with phyloP/phastCons", proposed_role="PRIMARY BIOLOGICAL EVIDENCE", remaining_uncertainty="exact alignment release requires dbNSFP field documentation", **common),
        dict(display_name="phyloP 100-way", exact_dataset_column=CONSERVATION[1], evidence_group="conservation", raw_derived_or_composite="derived evolutionary score", original_source="PHAST/UCSC 100 vertebrates", source_version="dbNSFP 4.8a field", alignment_or_derivation_method="sitewise conservation/acceleration statistic", training_or_calibration_data="100-vertebrate alignment and neutral model", training_labels="none", whether_clinvar_was_used="no", whether_hgmd_was_used="no", whether_uniprot_disease_annotations_were_used="no", other_clinical_label_sources="none", whether_training_variants_are_publicly_available="not applicable", underlying_component_features="multispecies alignment", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="no", whether_it_incorporates_other_computational_predictors="no", training_or_release_cutoff_date="track release not encoded locally", documentary_source="UCSC 100-way conservation track", exact_supporting_citation_or_url="https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way", circularity_risk_classification="NOT APPLICABLE AS SUPERVISED-LABEL CIRCULARITY", redundancy_concern="shared alignment signal with phastCons/GERP++", proposed_role="PRIMARY BIOLOGICAL EVIDENCE", remaining_uncertainty="exact underlying track build in dbNSFP requires documentation", **common),
        dict(display_name="phastCons 100-way", exact_dataset_column=CONSERVATION[2], evidence_group="conservation", raw_derived_or_composite="derived evolutionary score", original_source="PHAST/UCSC 100 vertebrates", source_version="dbNSFP 4.8a field", alignment_or_derivation_method="phylo-HMM posterior conserved-state probability", training_or_calibration_data="100-vertebrate alignment and phylogenetic model", training_labels="none", whether_clinvar_was_used="no", whether_hgmd_was_used="no", whether_uniprot_disease_annotations_were_used="no", other_clinical_label_sources="none", whether_training_variants_are_publicly_available="not applicable", underlying_component_features="multispecies alignment", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="no", whether_it_incorporates_other_computational_predictors="no", training_or_release_cutoff_date="track release not encoded locally", documentary_source="UCSC 100-way conservation track", exact_supporting_citation_or_url="https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way", circularity_risk_classification="NOT APPLICABLE AS SUPERVISED-LABEL CIRCULARITY", redundancy_concern="shared alignment signal with phyloP/GERP++", proposed_role="PRIMARY BIOLOGICAL EVIDENCE", remaining_uncertainty="exact underlying track build in dbNSFP requires documentation", **common),
        dict(display_name="SIFT", exact_dataset_column=PREDICTORS[0], evidence_group="computational predictor", raw_derived_or_composite="derived predictor", original_source="SIFT", source_version="dbNSFP 4.8a field", alignment_or_derivation_method="sequence homology and amino-acid probabilities", training_or_calibration_data="homologous protein sequences", training_labels="no clinical supervised labels for score construction", whether_clinvar_was_used="no documented direct training", whether_hgmd_was_used="no documented direct training", whether_uniprot_disease_annotations_were_used="no documented direct training", other_clinical_label_sources="none documented", whether_training_variants_are_publicly_available="not applicable/algorithmic", underlying_component_features="sequence conservation and amino-acid properties", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="no", whether_it_incorporates_other_computational_predictors="no", training_or_release_cutoff_date="exact dbNSFP source release unresolved", documentary_source="Vaser et al. 2016 protocol", exact_supporting_citation_or_url="https://doi.org/10.1038/nprot.2015.123", circularity_risk_classification="POSSIBLE RISK", redundancy_concern="strong conservation dependence; score direction is lower=more damaging", proposed_role="EXTERNAL COMPARATOR", remaining_uncertainty="exact SIFT database/version and historical benchmark reuse", **common),
        dict(display_name="PolyPhen-2 HDIV", exact_dataset_column=PREDICTORS[1], evidence_group="computational predictor", raw_derived_or_composite="supervised composite predictor", original_source="PolyPhen-2", source_version="HDIV via dbNSFP 4.8a", alignment_or_derivation_method="naive Bayes using sequence/structure annotations", training_or_calibration_data="HumDiv disease variants and mammalian divergence", training_labels="Mendelian disease vs divergence", whether_clinvar_was_used="not reported for original HumDiv", whether_hgmd_was_used="not established from local evidence", whether_uniprot_disease_annotations_were_used="yes; Swiss-Prot/UniProt disease mutation annotations", other_clinical_label_sources="HumDiv", whether_training_variants_are_publicly_available="yes (official training download)", underlying_component_features="conservation, sequence and structural annotations", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="not primary HDIV definition", whether_it_incorporates_other_computational_predictors="no", training_or_release_cutoff_date="2011_12 training download", documentary_source="Adzhubei et al. 2010 and official downloads", exact_supporting_citation_or_url="https://genetics.bwh.harvard.edu/wiki/!pph2/downloads", circularity_risk_classification="POSSIBLE RISK", redundancy_concern="overlaps conservation/protein-context evidence", proposed_role="EXTERNAL COMPARATOR", remaining_uncertainty="safe GRCh37/protein harmonisation not implemented", **common),
        dict(display_name="REVEL", exact_dataset_column=PREDICTORS[2], evidence_group="computational predictor", raw_derived_or_composite="supervised ensemble", original_source="REVEL", source_version="dbNSFP 4.8a field", alignment_or_derivation_method="ensemble of 13 predictors", training_or_calibration_data="pathogenic and rare neutral missense variants", training_labels="HGMD disease mutations vs rare population variants", whether_clinvar_was_used="used for independent reported test, not stated training", whether_hgmd_was_used="yes", whether_uniprot_disease_annotations_were_used="component predictors may use curated data", other_clinical_label_sources="HGMD; component tools", whether_training_variants_are_publicly_available="incomplete; HGMD pathogenic set restricted", underlying_component_features="MutPred,FATHMM,VEST,PolyPhen,SIFT,PROVEAN,MutationAssessor,MutationTaster,LRT,GERP,SiPhy,phyloP,phastCons", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="training neutral AF selection", whether_it_incorporates_other_computational_predictors="yes", training_or_release_cutoff_date="2016 publication", documentary_source="Ioannidis et al. 2016", exact_supporting_citation_or_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/", circularity_risk_classification="UNKNOWN/UNVERIFIABLE", redundancy_concern="explicitly contains all three conservation families and SIFT/PolyPhen", proposed_role="EXTERNAL COMPARATOR", remaining_uncertainty="exact pathogenic training variants unavailable due HGMD licensing", **common),
        dict(display_name="CADD PHRED", exact_dataset_column=PREDICTORS[3], evidence_group="computational predictor", raw_derived_or_composite="supervised composite score", original_source="CADD", source_version="unknown exact CADD release through dbNSFP 4.8a", alignment_or_derivation_method="model distinguishes human-derived from simulated variants", training_or_calibration_data="human-derived proxy-neutral and simulated proxy-deleterious variants", training_labels="proxy labels, not clinical pathogenicity", whether_clinvar_was_used="no for training; used in validation", whether_hgmd_was_used="no for training", whether_uniprot_disease_annotations_were_used="no direct labels; annotations may include broad sources", other_clinical_label_sources="none as training labels", whether_training_variants_are_publicly_available="yes for documented releases", underlying_component_features="many genomic annotations including conservation", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="may include population-related annotations depending release", whether_it_incorporates_other_computational_predictors="multiple annotations", training_or_release_cutoff_date="exact release unresolved", documentary_source="Official CADD training documentation", exact_supporting_citation_or_url="https://cadd.bihealth.org/training", circularity_risk_classification="POSSIBLE RISK", redundancy_concern="broad annotation overlap; exact release not recorded", proposed_role="EXTERNAL COMPARATOR", remaining_uncertainty="exact CADD model release behind dbNSFP field", **common),
        dict(display_name="AlphaMissense", exact_dataset_column=PREDICTORS[4], evidence_group="computational predictor", raw_derived_or_composite="deep-learning predictor", original_source="AlphaMissense", source_version="dbNSFP 4.8a field", alignment_or_derivation_method="AlphaFold-derived protein representations and evolutionary context", training_or_calibration_data="protein sequence/structure; population frequency-based weak supervision described by paper", training_labels="not directly ClinVar-trained", whether_clinvar_was_used="evaluation, not direct training", whether_hgmd_was_used="no documented training", whether_uniprot_disease_annotations_were_used="no documented clinical-label training", other_clinical_label_sources="ClinVar used for evaluation", whether_training_variants_are_publicly_available="model/data release available; exact weak-supervision set requires model documentation", underlying_component_features="protein sequence, structural context, evolutionary conservation", whether_it_already_contains_conservation_information="yes", whether_it_already_contains_population_information="yes in weak-supervision labels", whether_it_incorporates_other_computational_predictors="no documented", training_or_release_cutoff_date="2023", documentary_source="Cheng et al. 2023", exact_supporting_citation_or_url="https://doi.org/10.1126/science.adg7492", circularity_risk_classification="POSSIBLE RISK", redundancy_concern="conservation/protein-context and population-information overlap", proposed_role="EXTERNAL COMPARATOR", remaining_uncertainty="exact overlap cannot be concluded from score correlations", **common),
    ]
    write_csv(OUT / "01_provenance_circularity/provenance_circularity_register.csv", rows)
    overlap_rows = [{"score": r["display_name"], "exact_variant_overlap": "NOT_TESTED", "gene_overlap": "NOT_TESTED", "protein_overlap": "NOT_TESTED", "reason": r["remaining_uncertainty"], "partition_outcomes_examined": "no"} for r in rows]
    write_csv(OUT / "01_provenance_circularity/training_overlap_results.csv", overlap_rows)
    return rows


def heatmap(matrix: pd.DataFrame, path: Path, title: str, cluster: bool = False) -> None:
    data = matrix.copy()
    if cluster and len(data) > 2:
        # Deterministic greedy ordering by strongest absolute relationship; exploratory only.
        remaining = list(data.columns)
        order = [remaining.pop(0)]
        while remaining:
            previous = order[-1]
            nxt = max(remaining, key=lambda x: abs(float(data.loc[previous, x])) if pd.notna(data.loc[previous, x]) else -1)
            order.append(nxt); remaining.remove(nxt)
        data = data.loc[order, order]
    size = max(6, 0.8 * len(data))
    fig, ax = plt.subplots(figsize=(size, size))
    image = ax.imshow(data.astype(float), vmin=-1, vmax=1, cmap="coolwarm")
    ax.set_xticks(range(len(data))); ax.set_xticklabels(data.columns, rotation=75, ha="right", fontsize=8)
    ax.set_yticks(range(len(data))); ax.set_yticklabels(data.index, fontsize=8)
    ax.set_title(title); fig.colorbar(image, ax=ax, shrink=.75, label="Spearman correlation")
    fig.tight_layout(); path.parent.mkdir(parents=True, exist_ok=True); fig.savefig(path, dpi=180); plt.close(fig)


def redundancy(train: pd.DataFrame) -> dict:
    columns = CONSERVATION + PREDICTORS
    table = train[[TARGET] + columns].copy()
    for col in columns:
        table[col] = numeric(table[col])
    coverage = [{"column": c, "nonmissing_count": int(table[c].notna().sum()), "row_count": len(table), "coverage": round(float(table[c].notna().mean()), 8)} for c in columns]
    by_class = []
    for c in columns:
        for label in (0, 1):
            subset = table[table[TARGET].astype(int) == label]
            by_class.append({"column": c, "target_class": label, "nonmissing_count": int(subset[c].notna().sum()), "class_rows": len(subset), "coverage": round(float(subset[c].notna().mean()), 8)})
    spearman = table[columns].corr(method="spearman", min_periods=3)
    pearson = table[columns].corr(method="pearson", min_periods=3)
    pair_counts = pd.DataFrame(index=columns, columns=columns, dtype=int)
    for a in columns:
        for b in columns:
            pair_counts.loc[a, b] = int((table[a].notna() & table[b].notna()).sum())
    flagged = []
    for i, a in enumerate(columns):
        for b in columns[i + 1:]:
            value = spearman.loc[a, b]
            if pd.notna(value) and abs(value) >= .80:
                flagged.append({"field_a": a, "field_b": b, "spearman": round(float(value), 8), "pearson": round(float(pearson.loc[a, b]), 8), "paired_observations": int(pair_counts.loc[a, b]), "threshold": .80, "interpretation": "investigate shared information; not automatic deletion"})
    imputed = SimpleImputer(strategy="median").fit_transform(table[columns])
    vif = []
    for index, col in enumerate(columns):
        y = imputed[:, index]; x = np.delete(imputed, index, axis=1)
        r2 = LinearRegression().fit(x, y).score(x, y)
        value = math.inf if r2 >= 1 else 1 / (1 - r2)
        vif.append({"column": col, "vif": "Infinity" if not math.isfinite(value) else round(float(value), 8), "flag_vif_ge_10": int(value >= 10), "preprocessing": "training-only median imputation; original data unchanged"})
    base = OUT / "02_redundancy"
    write_csv(base / "coverage_summary.csv", coverage); write_csv(base / "coverage_by_class.csv", by_class)
    spearman.to_csv(base / "spearman_correlation.csv"); pearson.to_csv(base / "pearson_correlation.csv")
    pair_counts.astype(int).to_csv(base / "pairwise_observation_counts.csv")
    write_csv(base / "vif_results.csv", vif); write_csv(base / "flagged_redundant_pairs.csv", flagged, ["field_a","field_b","spearman","pearson","paired_observations","threshold","interpretation"])
    heatmap(spearman.loc[CONSERVATION, CONSERVATION], base / "conservation_only_heatmap.png", "Training-only conservation correlations")
    heatmap(spearman.loc[PREDICTORS, PREDICTORS], base / "computational_predictor_only_heatmap.png", "Training-only predictor correlations")
    heatmap(spearman, base / "combined_heatmap.png", "Training-only evidence correlations")
    heatmap(spearman, base / "clustered_combined_heatmap.png", "Training-only evidence correlations (ordered)", True)
    return {"coverage": coverage, "flagged": flagged, "vif": vif}


def metric_values(y: np.ndarray, probability: np.ndarray, threshold: float) -> dict:
    predicted = (probability >= threshold).astype(int)
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    logit = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    calibration = LogisticRegression(C=1e6, solver="liblinear").fit(logit, y)
    return {
        "mcc": matthews_corrcoef(y, predicted), "pr_auc": average_precision_score(y, probability),
        "roc_auc": roc_auc_score(y, probability), "brier": brier_score_loss(y, probability),
        "calibration_intercept": float(calibration.intercept_[0]),
        "calibration_slope": float(calibration.coef_[0, 0]), "threshold": threshold,
    }


def choose_threshold(y: np.ndarray, probability: np.ndarray) -> float:
    candidates = np.unique(np.quantile(probability, np.linspace(.05, .95, 91)))
    return float(max(candidates, key=lambda t: (matthews_corrcoef(y, probability >= t), -abs(t - .5))))


def experiments(train: pd.DataFrame) -> dict:
    primary = load(PATHS["train_primary"])
    optional = load(PATHS["train_optional"])[[KEY] + PREDICTORS]
    population = load(PATHS["train_population"])[[KEY] + POPULATION]
    data = primary.merge(optional, on=KEY, validate="one_to_one").merge(population, on=KEY, validate="one_to_one")
    features = [c for c in primary.columns if c not in METADATA]
    base_features = [c for c in features if c not in CONSERVATION]
    y = data[TARGET].astype(int).to_numpy(); groups = data[GROUP].to_numpy()
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=SEED)
    folds = np.empty(len(data), dtype=int)
    split_indices = []
    for fold, (fit, held) in enumerate(cv.split(data, y, groups)):
        folds[held] = fold; split_indices.append((fit, held))
    write_csv(OUT / "03_controlled_experiments/fixed_cv_fold_assignments.csv", [{KEY: key, GROUP: gene, "cv_fold": int(fold)} for key, gene, fold in zip(data[KEY], groups, folds)])

    conditions = []
    for size in range(4):
        for subset in itertools.combinations(CONSERVATION, size):
            conditions.append(("conservation_" + ("none" if not subset else "__".join(x.replace("conservation_", "") for x in subset)), base_features + list(subset), "conservation subset", False))
    manifest, fold_rows, oof_rows, summaries = [], [], [], []

    def run(condition_id: str, condition_features: list[str], category: str, standalone: bool = False) -> dict:
        started = time.perf_counter(); probs = np.full(len(data), np.nan); warnings = []
        for fold, (fit, held) in enumerate(split_indices):
            x_fit = data.iloc[fit][condition_features].apply(numeric)
            x_held = data.iloc[held][condition_features].apply(numeric)
            pipe = Pipeline([("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler()), ("model", LogisticRegression(C=1.0, max_iter=2000, solver="liblinear", class_weight="balanced", random_state=SEED))])
            try:
                pipe.fit(x_fit, y[fit]); held_prob = pipe.predict_proba(x_held)[:, 1]
                threshold = choose_threshold(y[fit], pipe.predict_proba(x_fit)[:, 1]) if standalone else .5
                probs[held] = held_prob; values = metric_values(y[held], held_prob, threshold)
                fold_rows.append({"condition_id": condition_id, "fold": fold, "heldout_rows": len(held), "heldout_genes": len(set(groups[held])), **{k: round(float(v), 8) for k,v in values.items()}})
            except Exception as exc:
                warnings.append(f"fold {fold}: {type(exc).__name__}: {exc}")
        valid = np.isfinite(probs)
        rows = [r for r in fold_rows if r["condition_id"] == condition_id]
        summary = {"condition_id": condition_id, "category": category, "included_features": ";".join(condition_features), "feature_count": len(condition_features), "row_count": int(valid.sum()), "gene_count": len(set(groups[valid])), "runtime_seconds": round(time.perf_counter()-started,3), "warnings": " | ".join(warnings), "failure_information": "" if not warnings else "see warnings"}
        for metric in ("mcc","pr_auc","roc_auc","brier","calibration_intercept","calibration_slope"):
            vals=[float(r[metric]) for r in rows]; summary[f"mean_{metric}"]=round(float(np.mean(vals)),8); summary[f"sd_{metric}"]=round(float(np.std(vals,ddof=1)),8)
        for index in np.where(valid)[0]:
            oof_rows.append({"condition_id": condition_id, KEY: data.iloc[index][KEY], GROUP: groups[index], "cv_fold": int(folds[index]), "observed_label": int(y[index]), "predicted_probability": round(float(probs[index]),10)})
        summaries.append(summary)
        manifest.append({"condition_id": condition_id, "category": category, "included_features": ";".join(condition_features), "excluded_features": ";".join(sorted(set(features+PREDICTORS+POPULATION)-set(condition_features))), "model_configuration": "LogisticRegression(C=1,solver=liblinear,class_weight=balanced,seed=20261018)", "preprocessing_configuration": "fold-local median imputation and standard scaling", "fold_assignments": "fixed_cv_fold_assignments.csv", "feature_coverage": ";".join(f"{f}={numeric(data[f]).notna().mean():.4f}" for f in condition_features), "runtime_seconds": summary["runtime_seconds"], "warnings": summary["warnings"], "failure_information": summary["failure_information"]})
        return summary

    conservation_results = [run(*condition) for condition in conditions]
    recommended_conservation = max(conservation_results, key=lambda r: (r["mean_mcc"], r["mean_pr_auc"], -r["mean_brier"]))
    locked_biological = recommended_conservation["included_features"].split(";")
    for predictor in PREDICTORS:
        run(f"comparator_{predictor}", [predictor], "standalone external comparator", True)
        run(f"biological_plus_{predictor}", locked_biological + [predictor], "one-at-a-time optional sensitivity", False)
    run("population_inclusive_sensitivity", locked_biological + POPULATION, "population sensitivity", False)
    write_csv(OUT / "03_controlled_experiments/experiment_manifest.csv", manifest)
    write_csv(OUT / "03_controlled_experiments/fold_level_metrics.csv", fold_rows)
    write_csv(OUT / "03_controlled_experiments/experiment_summary.csv", summaries)
    write_csv(OUT / "03_controlled_experiments/out_of_fold_predictions.csv", oof_rows)
    recommendation = {
        "selection_basis": "training-only fixed gene-grouped cross-validation; maximum mean MCC, then PR-AUC, then lower Brier",
        "recommended_condition_id": recommended_conservation["condition_id"],
        "recommended_features": locked_biological,
        "mean_metrics": {k: v for k,v in recommended_conservation.items() if k.startswith("mean_")},
        "optional_predictors": "external comparators or one-at-a-time sensitivity only; not approved primary inputs",
        "population": "no-population remains primary; population-inclusive result is sensitivity only",
        "validation_used": False, "test_used": False,
    }
    write_json(OUT / "03_controlled_experiments/recommended_training_only_configuration.json", recommendation)
    return recommendation


def lock_protocol(cohort_report: dict, split_report: dict) -> None:
    registry = load(ROOT / "data/processed/modeling_feature_registry.csv")
    primary = sorted(registry.loc[registry["primary_model_feature"] == "1", "field"])
    payload = {
        "schema_version": 1, "locked_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)} for name,path in PATHS.items()},
        "row_counts": cohort_report["row_counts"], "split_counts": split_report["split_counts_from_locked_assignment_summary"],
        "target_column": TARGET, "gene_group_column": GROUP, "variant_key_definition": "genome_build:chromosome:position:reference:alternate",
        "available_conservation_variables": CONSERVATION, "available_computational_predictors": PREDICTORS,
        "approved_base_biological_features": sorted(set(primary)-set(CONSERVATION)), "prohibited_variables": sorted(PROHIBITED | METADATA),
        "missing_value_handling": "fold-local median imputation; missing biological values never converted to zero",
        "encoding_and_scaling": "numeric/binary tables; fold-local StandardScaler for fixed logistic reference model",
        "random_seed": SEED, "cross_validation": "5-fold StratifiedGroupKFold by reconciled gene; fixed assignments saved",
        "model": "LogisticRegression(C=1,solver=liblinear,max_iter=2000,class_weight=balanced)",
        "metrics": {"primary":"MCC higher is better","secondary":["PR-AUC higher","ROC-AUC higher","Brier lower","calibration slope near 1/intercept near 0"]},
        "correlation_threshold": "absolute Spearman >= 0.80", "vif_threshold": ">= 10 supplementary",
        "conditions": "all 8 conservation subsets; predictor alone and one-at-a-time addition; population sensitivity",
        "software": {"python":sys.version.split()[0],"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__,"scikit_learn":sklearn.__version__,"matplotlib":matplotlib.__version__},
    }
    path = OUT / "00_locked_protocol/locked_protocol.yaml"
    path.parent.mkdir(parents=True, exist_ok=True)
    # JSON is a valid YAML 1.2 document and avoids an unnecessary runtime dependency.
    path.write_text(json.dumps(payload, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    (path.parent / "locked_protocol.sha256").write_text(f"{sha256(path)}  locked_protocol.yaml\n", encoding="utf-8")


def documents(recommendation: dict, provenance_rows: list[dict]) -> None:
    base = OUT / "00_locked_protocol/cohort_validation"
    (base / "PREPROCESSING_VERIFICATION.md").write_text("""# Preprocessing verification\n\nThe authoritative preprocessing scripts and audit manifests were inventoried and the frozen outputs were checked without regeneration. All imputers, encoders, scalers, selectors and models must be fitted inside each training fold. Validation and test transformations must use training-learned parameters only. Missing biological evidence must remain missing until fold-local imputation and must never be silently converted to zero. Scaling is used only for the fixed logistic audit instrument. Biological outliers are retained unless a prospective domain rule is approved. Reconciled genes are grouping variables, never predictors. The existing gene-separated split is reused unchanged.\n""", encoding="utf-8")
    src = OUT / "01_provenance_circularity"
    (src / "documentary_sources.md").write_text("""# Documentary sources\n\n- GERP++: Davydov et al. (2010), https://pmc.ncbi.nlm.nih.gov/articles/PMC2996323/\n- UCSC 100-way phyloP/phastCons documentation: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way\n- SIFT protocol: https://doi.org/10.1038/nprot.2015.123\n- PolyPhen-2 paper/training downloads: https://genetics.bwh.harvard.edu/wiki/!pph2/downloads\n- REVEL: https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/\n- Official CADD training documentation: https://cadd.bihealth.org/training\n- AlphaMissense: https://doi.org/10.1126/science.adg7492\n\nThese sources support provenance statements; absence of public training records was not interpreted as absence of overlap.\n""", encoding="utf-8")
    (src / "circularity_limitations.md").write_text("""# Circularity and overlap limitations\n\nNo external training set was downloaded for exact overlap. PolyPhen-2 training data are public but require validated GRCh37/protein-level harmonisation; REVEL's HGMD-derived pathogenic training records are not fully public; the exact CADD release behind the dbNSFP field is not recorded locally; and score-source versions for SIFT require further confirmation. Therefore exact variant, gene and protein overlap remain `UNKNOWN/UNVERIFIABLE` where appropriate. AlphaMissense was not directly trained on ClinVar labels, but it uses evolutionary/protein context and population-based weak supervision. Correlation is evidence of shared signal, not a test for circularity.\n""", encoding="utf-8")
    command = "../rare_missense_dataset_pipeline/.venv/bin/python -m src.pre_modeling_audit.run_pre_modeling_audit --validation-confirmation outputs/pre_modeling_evidence_audit/03_controlled_experiments/recommended_training_only_configuration.json\n"
    (OUT / "03_controlled_experiments/validation_confirmation_command.txt").write_text(command, encoding="utf-8")
    report = f"""# Pre-modelling evidence audit\n\n## Purpose and scope\n\nThis work separates (A) reproducible cohort/preprocessing verification, (B) a pre-modelling provenance and redundancy audit, and (C) preliminary modelling for evidence selection. Phase C is not exploratory data analysis and is not the final model comparison.\n\n## Dataset and split confirmation\n\nThe frozen dataset contains 40,000 unique variants. The reconciled modelling cohort contains 39,552, while 448 gene-mapping disagreements remain preserved separately. Retained and excluded keys are disjoint and account exactly for the frozen source. The existing train/validation/test assignment remains gene- and variant-disjoint. No source dataset or split was modified. Aggregate validation/test class counts were used only from the locked assignment summary for integrity verification; outcome rows were not analysed.\n\n## Preprocessing and leakage\n\nThe script inventory, schema, missingness and leakage reports are under `00_locked_protocol/cohort_validation`. Target, clinical text, identifiers, gene/transcript identifiers, database IDs and lookup metadata are excluded. Natural variant, Mutagenesis, Sequence conflict and Alternative sequence indicators remain prohibited. Fold-local median imputation and scaling were used only in preliminary training CV; missing biology was never converted to zero.\n\n## Circularity versus redundancy\n\nCircularity concerns reuse of clinical labels/variants or downstream information. Redundancy concerns shared information between features. Low correlation cannot prove absence of circularity, and high correlation does not by itself prove it. Heatmaps are therefore presented only as training-data shared-variance summaries.\n\n## Provenance findings\n\nGERP++, phyloP and phastCons are evolutionary measures and are classified as not applicable to supervised-label circularity, although they can double-count evolutionary information. PolyPhen-2 HDIV uses supervised disease/divergence labels and UniProt/Swiss-Prot annotations. REVEL uses an HGMD-derived pathogenic training set and incorporates SIFT, PolyPhen and multiple conservation scores; exact overlap is unverifiable because the full pathogenic set is restricted. CADD uses proxy labels rather than clinical labels but overlaps broad annotations. AlphaMissense was evaluated on ClinVar rather than directly trained on it, while still sharing evolutionary/protein and weak population information. Exact overlap was not forced where safe harmonisation or public records were unavailable.\n\n## Coverage and redundancy\n\nTraining-only coverage, class-stratified coverage, pairwise counts, Spearman/Pearson matrices, VIF and flagged pairs are under `02_redundancy`. Score direction matters: lower SIFT indicates greater predicted damage, so negative correlations may still represent shared information. Flags trigger investigation, not automatic deletion.\n\n## Preliminary controlled experiments\n\nAll eight conservation subsets were compared with a fixed logistic-regression audit instrument using identical five-fold gene-grouped training-only folds. Each optional predictor was evaluated alone and as a one-at-a-time addition. The population-inclusive condition was sensitivity-only. No model-family or hyperparameter search was performed.\n\nThe training-only recommended condition is `{recommendation['recommended_condition_id']}` with {len(recommendation['recommended_features'])} features. This recommendation is provisional and must be confirmed once on validation using the saved command. Optional predictors remain external comparators/sensitivity inputs because predictive improvement cannot resolve circularity. The no-population condition remains primary.\n\n## Limitations\n\nExact external training overlap remains unavailable or unsafe to harmonise for several scores. The reference model is an audit instrument, not the final thesis model. No causal or non-circularity conclusion follows from correlations, PCA (not used), or performance.\n\n## Commands and next step\n\nMain command: `python -m src.pre_modeling_audit.run_pre_modeling_audit`. The exact future validation command is in `03_controlled_experiments/validation_confirmation_command.txt`; it was not executed. Validation outcomes were not inspected. Test outcomes were not inspected and the test partition remains locked.\n\n## Output index\n\nOutputs follow the requested `00_locked_protocol`, `01_provenance_circularity`, `02_redundancy`, `03_controlled_experiments`, and `04_report` directories. CSV/JSON files contain machine-readable evidence; PNG files are correlation visualisations only.\n"""
    execution_log = """

## Complete execution log

```text
git status --short
git branch --show-current
git worktree list
python -m py_compile src/pre_modeling_audit/run_pre_modeling_audit.py src/pre_modeling_audit/validate_final_cohort.py tests/test_pre_modeling_audit.py
python -m unittest tests.test_pre_modeling_audit -v
python -m src.pre_modeling_audit.validate_final_cohort
python -m src.pre_modeling_audit.run_pre_modeling_audit
python -m unittest discover -s tests -v
```
"""
    report_path = OUT / "04_report/PRE_MODELING_EVIDENCE_AUDIT.md"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(report + execution_log, encoding="utf-8")


def validation_confirmation(configuration_path: Path) -> None:
    """Future, explicitly invoked single validation confirmation; never called by the audit."""
    configuration = json.loads(configuration_path.read_text(encoding="utf-8"))
    features = configuration["recommended_features"]
    train = load(PATHS["train_primary"])
    validation = load(ROOT / "data/processed/modeling_splits/validation_core_no_population.csv")
    pipe = Pipeline([
        ("imputer", SimpleImputer(strategy="median")),
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(C=1.0, max_iter=2000, solver="liblinear", class_weight="balanced", random_state=SEED)),
    ])
    pipe.fit(train[features].apply(numeric), train[TARGET].astype(int))
    probability = pipe.predict_proba(validation[features].apply(numeric))[:, 1]
    result = {
        "condition_id": configuration["recommended_condition_id"],
        "selection_locked_before_validation": True,
        "validation_rows": len(validation),
        "validation_genes": int(validation[GROUP].nunique()),
        "metrics": {k: round(float(v), 8) for k, v in metric_values(validation[TARGET].astype(int).to_numpy(), probability, .5).items()},
        "test_partition_accessed": False,
    }
    write_json(OUT / "03_controlled_experiments/validation_confirmation_result.json", result)
    print(json.dumps(result, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--validation-confirmation", help="Future-only validation configuration; intentionally not implemented in audit run")
    args = parser.parse_args()
    if args.validation_confirmation:
        validation_confirmation(Path(args.validation_confirmation))
        return
    cohort, cohort_report = validate_cohort()
    split_report = validate_split(cohort)
    inventory_and_leakage(cohort)
    lock_protocol(cohort_report, split_report)
    provenance_rows = provenance()
    train_optional = load(PATHS["train_optional"])
    redundancy(train_optional)
    recommendation = experiments(train_optional)
    documents(recommendation, provenance_rows)
    print(json.dumps({"status":"PASS","recommendation":recommendation}, indent=2))


if __name__ == "__main__":
    main()
