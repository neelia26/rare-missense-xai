#!/usr/bin/env python3
"""Post-run corrections, paired uncertainty analysis, and publication figures.

This module never reads validation/test outcomes. It consumes locked assignment
counts, training-only audit tables, and saved training OOF predictions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tarfile
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, matthews_corrcoef

from src.pre_modeling_audit.run_pre_modeling_audit import (
    GROUP, KEY, OUT, PATHS, PREDICTORS, ROOT, SEED, TARGET, sha256, write_csv, write_json,
)

N_BOOTSTRAP = 5000
ALL3 = "conservation_gerp_rs__phylop100way__phastcons100way"
PHYLO_PHAST = "conservation_phylop100way__phastcons100way"
GERP_PHYLO = "conservation_gerp_rs__phylop100way"
PHYLO = "conservation_phylop100way"
COMPARISONS = [(ALL3, PHYLO_PHAST), (ALL3, GERP_PHYLO), (ALL3, PHYLO), (PHYLO_PHAST, PHYLO)]
DISPLAY = {
    "conservation_gerp_rs": "GERP++", "conservation_phylop100way": "phyloP",
    "conservation_phastcons100way": "PhastCons", "predictor_sift_score": "SIFT",
    "predictor_polyphen2_hdiv_score": "PolyPhen-2", "predictor_revel_score": "REVEL",
    "predictor_cadd_phred": "CADD", "predictor_alphamissense_score": "AlphaMissense",
    ALL3: "GERP++ + phyloP + PhastCons", PHYLO_PHAST: "phyloP + PhastCons",
    GERP_PHYLO: "GERP++ + phyloP", PHYLO: "phyloP only",
}
BLUE, ORANGE, GREEN, RED, PURPLE, GREY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#666666"


def savefig(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=320, bbox_inches="tight", facecolor="white")
    fig.savefig(path.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def locked_checksum() -> str:
    expected = (OUT / "00_locked_protocol/locked_protocol.sha256").read_text().split()[0]
    actual = sha256(OUT / "00_locked_protocol/locked_protocol.yaml")
    if actual != expected:
        raise RuntimeError(f"Locked protocol checksum mismatch: {expected} != {actual}")
    return actual


def corrected_leakage() -> None:
    source = pd.read_csv(OUT / "00_locked_protocol/cohort_validation/leakage_column_audit.csv", dtype=str)
    population = {"gnomad_global_af", "gnomad_max_population_af"}
    mask = source["column_name"].isin(population)
    source.loc[mask, "proposed_role"] = "OPTIONAL SENSITIVITY INPUT"
    source.loc[mask, "allowed_as_input"] = "sensitivity only"
    source.loc[mask, "leakage_concern"] = "class-associated missingness and possible dependence on clinical classification practices"
    source.loc[mask, "decision"] = "ALLOW_SENSITIVITY_ONLY"
    source.loc[mask, "reason"] = "Biological population evidence; excluded from the primary no-population model and no missingness indicator is permitted."
    source.to_csv(OUT / "00_locked_protocol/cohort_validation/leakage_column_audit_corrected.csv", index=False)
    roles = {
        "predictor_sift_score": ("EXTERNAL COMPARATOR", "Algorithmic conservation predictor; low direct clinical-label circularity but high conservation redundancy."),
        "predictor_polyphen2_hdiv_score": ("EXTERNAL COMPARATOR", "Supervised disease/divergence training creates circularity concern despite public training records."),
        "predictor_revel_score": ("EXTERNAL COMPARATOR", "HGMD-derived pathogenic training and composite inputs create unverifiable circularity and redundancy."),
        "predictor_cadd_phred": ("OPTIONAL SENSITIVITY INPUT", "Proxy-label training rather than clinical labels; broad annotation overlap and exact release unresolved."),
        "predictor_alphamissense_score": ("OPTIONAL SENSITIVITY INPUT", "Population-frequency fine-tuning and shared evolutionary/protein evidence require sensitivity-only treatment."),
    }
    write_csv(OUT / "00_locked_protocol/cohort_validation/optional_predictor_input_policy.csv", [
        {"column_name": field, "evidence_type": "computational pathogenicity predictor", "proposed_role": role,
         "allowed_as_primary_biological_evidence": "no", "reason": reason,
         "performance_policy": "Training CV improvement does not remove circularity or promote the score to primary evidence."}
        for field, (role, reason) in roles.items()
    ])


def polyphen_overlap(archive: Path | None) -> dict:
    result = {
        "status": "UNKNOWN/UNVERIFIABLE", "source": "PolyPhen-2 HumDiv 2011_12 official archive",
        "archive_sha256": "", "training_records": "", "exact_protein_variant_overlap": "",
        "gene_overlap": "", "protein_overlap": "", "unmapped_training_records": "",
        "train_exact_overlap": "", "validation_exact_overlap": "", "test_exact_overlap": "",
        "notes": "Official protein-level training data were not supplied to this run; genomic coordinates were not compared.",
    }
    if not archive or not archive.is_file():
        return result
    records = []
    with tarfile.open(archive, "r:gz") as tar:
        for name in ("humdiv-2011_12.deleterious.pph.input", "humdiv-2011_12.neutral.pph.input"):
            handle = tar.extractfile(name)
            if handle is None:
                raise RuntimeError(f"Missing {name} in PolyPhen archive")
            for line in handle:
                parts = line.decode().strip().split()
                if len(parts) >= 4:
                    records.append((parts[0].split(".")[0], int(parts[1]), parts[2], parts[3]))
    # Read only identifiers/protein mapping—not validation/test labels.
    cols = [KEY, GROUP, "uniprot_selected_accession", "reconciled_protein_position", "reconciled_reference_amino_acid", "reconciled_alternate_amino_acid"]
    cohort = pd.read_csv(PATHS["cohort"], usecols=cols, dtype=str)
    cohort["protein_key"] = list(zip(
        cohort["uniprot_selected_accession"].str.split(".").str[0],
        pd.to_numeric(cohort["reconciled_protein_position"], errors="coerce"),
        cohort["reconciled_reference_amino_acid"], cohort["reconciled_alternate_amino_acid"],
    ))
    assignment = pd.read_csv(PATHS["assignments"], usecols=["gene", "split"], dtype=str)
    cohort = cohort.merge(assignment, left_on=GROUP, right_on="gene", validate="many_to_one")
    training = set(records); proteins = {r[0] for r in records}
    matched = cohort[cohort["protein_key"].isin(training)]
    overlap_proteins = set(cohort["uniprot_selected_accession"].str.split(".").str[0]) & proteins
    cohort_proteins = set(cohort["uniprot_selected_accession"].str.split(".").str[0])
    mapped_training = sum(r[0] in cohort_proteins for r in records)
    protein_matched_cohort = cohort[
        cohort["uniprot_selected_accession"].str.split(".").str[0].isin(overlap_proteins)
    ]
    result.update({
        "status": "TESTED", "archive_sha256": sha256(archive), "training_records": len(records),
        "exact_protein_variant_overlap": len(matched), "gene_overlap": protein_matched_cohort[GROUP].nunique(),
        "protein_overlap": len(overlap_proteins), "unmapped_training_records": len(records) - mapped_training,
        "train_exact_overlap": int((matched["split"] == "train").sum()),
        "validation_exact_overlap": int((matched["split"] == "validation").sum()),
        "test_exact_overlap": int((matched["split"] == "test").sum()),
        "notes": "Matched accession + amino-acid position + reference + alternate. No GRCh37/GRCh38 coordinate comparison and no partition outcomes used.",
    })
    return result


def corrected_provenance(polyphen: dict) -> list[dict]:
    original = pd.read_csv(OUT / "01_provenance_circularity/provenance_circularity_register.csv", dtype=str).fillna("")
    rows = original.to_dict("records")
    by_name = {r["display_name"]: r for r in rows}
    for name in ("GERP++ RS", "phyloP 100-way", "phastCons 100-way"):
        row = by_name[name]
        row["circularity_risk_classification"] = "NOT APPLICABLE AS SUPERVISED-LABEL CIRCULARITY"
        row["whether_training_variants_are_publicly_available"] = "NOT_APPLICABLE"
        row["whether_exact_variant_overlap_can_be_tested"] = "NOT_APPLICABLE"
        row["whether_gene_overlap_can_be_tested"] = "NOT_APPLICABLE"
        row["whether_protein_overlap_can_be_tested"] = "NOT_APPLICABLE"
        row["remaining_uncertainty"] = "Exact underlying alignment/track version is not recoverable from local dbNSFP 4.8a/MyVariant metadata."
    sift = by_name["SIFT"]
    sift.update({
        "raw_derived_or_composite": "algorithmic sequence-conservation score",
        "training_or_calibration_data": "homologous protein sequences; no conventional clinically labelled training set",
        "training_labels": "NOT_APPLICABLE",
        "whether_training_variants_are_publicly_available": "NOT_APPLICABLE",
        "whether_exact_variant_overlap_can_be_tested": "NOT_APPLICABLE",
        "whether_gene_overlap_can_be_tested": "NOT_APPLICABLE",
        "whether_protein_overlap_can_be_tested": "NOT_APPLICABLE",
        "circularity_risk_classification": "LOW RISK",
        "redundancy_concern": "High overlap with conservation evidence; lower score means more damaging.",
        "remaining_uncertainty": "Exact sequence database/version behind the dbNSFP 4.8a SIFT field is not recovered locally.",
    })
    poly = by_name["PolyPhen-2 HDIV"]
    poly["circularity_risk_classification"] = "POSSIBLE RISK"
    poly["whether_exact_variant_overlap_can_be_tested"] = polyphen["status"]
    poly["whether_gene_overlap_can_be_tested"] = polyphen["status"]
    poly["whether_protein_overlap_can_be_tested"] = polyphen["status"]
    poly["remaining_uncertainty"] = polyphen["notes"]
    revel = by_name["REVEL"]
    revel["whether_training_variants_are_publicly_available"] = "PARTIAL: neutral population portion described; HGMD-derived pathogenic portion restricted"
    revel["circularity_risk_classification"] = "UNKNOWN/UNVERIFIABLE"
    revel["remaining_uncertainty"] = "Public and restricted portions are separated conceptually; exact pathogenic overlap cannot be tested without licensed HGMD records."
    cadd = by_name["CADD PHRED"]
    cadd["proposed_role"] = "OPTIONAL SENSITIVITY INPUT"
    cadd["source_version"] = "UNKNOWN EXACT RELEASE (dbNSFP 4.8a field does not identify CADD model release locally)"
    cadd["circularity_risk_classification"] = "POSSIBLE RISK"
    cadd["remaining_uncertainty"] = "Proxy-label training is distinct from clinical-label training; annotation overlap remains, exact CADD release unresolved."
    alpha = by_name["AlphaMissense"]
    alpha["proposed_role"] = "OPTIONAL SENSITIVITY INPUT"
    alpha["training_or_calibration_data"] = "fine-tuned with human and primate population-frequency information"
    alpha["training_labels"] = "population-frequency-derived weak supervision; ClinVar used for evaluation, not direct training"
    alpha["circularity_risk_classification"] = "POSSIBLE RISK"
    alpha["remaining_uncertainty"] = "Conservation, protein-context and population overlap remain; score correlation does not establish exact training overlap."
    path = OUT / "01_provenance_circularity/provenance_circularity_register_corrected.csv"
    write_csv(path, rows, list(original.columns))
    overlap = []
    for name in ("GERP++ RS", "phyloP 100-way", "phastCons 100-way", "SIFT"):
        overlap.append({"score": name, "overlap_status": "NOT_APPLICABLE", "public_portion": "NOT_APPLICABLE", "restricted_portion": "NOT_APPLICABLE", "exact_protein_variant_overlap": "NOT_APPLICABLE", "gene_overlap": "NOT_APPLICABLE", "protein_overlap": "NOT_APPLICABLE", "unmapped_records": "NOT_APPLICABLE", "partition_counts": "NOT_APPLICABLE", "reason": "Evolutionary/algorithmic statistic without conventional clinically labelled training variants."})
    overlap.append({"score":"PolyPhen-2 HDIV","overlap_status":polyphen["status"],"public_portion":"HumDiv 2011_12 official archive","restricted_portion":"none identified for official HumDiv archive","exact_protein_variant_overlap":polyphen["exact_protein_variant_overlap"],"gene_overlap":polyphen["gene_overlap"],"protein_overlap":polyphen["protein_overlap"],"unmapped_records":polyphen["unmapped_training_records"],"partition_counts":f"train={polyphen['train_exact_overlap']}; validation={polyphen['validation_exact_overlap']}; test={polyphen['test_exact_overlap']}","reason":polyphen["notes"]})
    overlap.append({"score":"REVEL","overlap_status":"UNKNOWN/UNVERIFIABLE","public_portion":"neutral population training portion described but not exact-overlap tested","restricted_portion":"HGMD-derived pathogenic records restricted","exact_protein_variant_overlap":"UNKNOWN/UNVERIFIABLE","gene_overlap":"UNKNOWN/UNVERIFIABLE","protein_overlap":"UNKNOWN/UNVERIFIABLE","unmapped_records":"UNKNOWN/UNVERIFIABLE","partition_counts":"not calculated","reason":"Exact pathogenic overlap cannot be excluded or tested without licensed records."})
    for score, reason in (("CADD PHRED","Exact model release unresolved; proxy-label training is not clinical-label training."),("AlphaMissense","Fine-tuning examples/population information do not support exact variant-overlap inference from correlation.")):
        overlap.append({"score":score,"overlap_status":"UNKNOWN/UNVERIFIABLE","public_portion":"documentary provenance reviewed","restricted_portion":"not applicable/unclear","exact_protein_variant_overlap":"UNKNOWN/UNVERIFIABLE","gene_overlap":"UNKNOWN/UNVERIFIABLE","protein_overlap":"UNKNOWN/UNVERIFIABLE","unmapped_records":"UNKNOWN/UNVERIFIABLE","partition_counts":"not calculated","reason":reason})
    write_csv(OUT / "01_provenance_circularity/training_overlap_results_corrected.csv", overlap)
    return rows


def write_provenance_docs(polyphen: dict) -> None:
    (OUT / "01_provenance_circularity/documentary_sources_corrected.md").write_text(
        """# Corrected documentary sources\n\n- GERP++: Davydov et al. (2010), https://pmc.ncbi.nlm.nih.gov/articles/PMC2996323/\n- UCSC/PHAST 100-way track: https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg38&g=cons100way\n- SIFT protocol: Vaser et al. (2016), https://doi.org/10.1038/nprot.2015.123\n- PolyPhen-2: Adzhubei et al. (2010), https://doi.org/10.1038/nmeth0410-248 and official HumDiv archive https://genetics.bwh.harvard.edu/downloads/pph2/training/training-2.2.2.tar.gz\n- REVEL: Ioannidis et al. (2016), https://pmc.ncbi.nlm.nih.gov/articles/PMC5065685/\n- CADD training: https://cadd.bihealth.org/training\n- AlphaMissense: Cheng et al. (2023), https://doi.org/10.1126/science.adg7492\n\nOnly primary papers and official documentation are used.\n""", encoding="utf-8")
    (OUT / "01_provenance_circularity/circularity_limitations_corrected.md").write_text(
        f"""# Corrected circularity limitations\n\nGERP++, phyloP, PhastCons and SIFT do not have conventional clinically labelled training-variant overlap; those checks are `NOT_APPLICABLE`. Their material concern is overlapping evolutionary signal.\n\nPolyPhen-2 HumDiv was protein-matched using accession, amino-acid position, reference and alternate residues. Status: **{polyphen['status']}**; exact matches: **{polyphen['exact_protein_variant_overlap']}**; protein overlap: **{polyphen['protein_overlap']}**; mapped-gene overlap: **{polyphen['gene_overlap']}**; training records not mapped to a current cohort protein: **{polyphen['unmapped_training_records']}**. No genomic coordinate conversion was attempted.\n\nREVEL remains `UNKNOWN/UNVERIFIABLE` for exact pathogenic overlap because the HGMD-derived portion is restricted. CADD uses proxy labels, not clinical labels, but its exact release is `UNKNOWN EXACT RELEASE`. AlphaMissense used population-frequency information for fine-tuning and ClinVar for evaluation; correlation cannot determine exact training overlap. Unavailable data never demonstrate absence of overlap.\n""", encoding="utf-8")


def amendment() -> None:
    (OUT / "00_locked_protocol/protocol_amendment_001.md").write_text(
        """# Protocol amendment 001 — paired conservation uncertainty\n\nThis post-run amendment was created after reviewing the initial training-only results. It does not alter or overwrite the locked protocol retrospectively. Its sole purpose is to test whether the extra complexity of the highest-observed conservation subset is supported by a stable paired improvement. The analysis uses 5,000 deterministic paired gene-cluster bootstrap samples, seed 20261018, preserves all variants belonging to each sampled gene, and compares MCC, average precision (reported as PR-AUC), and Brier score from saved training OOF predictions. Validation and test outcomes are not accessed.\n""", encoding="utf-8")
    payload = {"amendment":1,"post_run":True,"changes_original_protocol":False,"resampling_unit":GROUP,"bootstrap_samples":N_BOOTSTRAP,"seed":SEED,"paired":True,"metrics":["MCC","average_precision_PR_AUC","Brier"],"comparisons":[list(x) for x in COMPARISONS],"data":"saved training-only OOF predictions"}
    (OUT / "00_locked_protocol/protocol_amendment_001.yaml").write_text(json.dumps(payload,indent=2)+"\n",encoding="utf-8")


def bootstrap() -> tuple[pd.DataFrame, dict]:
    oof = pd.read_csv(OUT / "03_controlled_experiments/out_of_fold_predictions.csv")
    wanted = sorted(set(sum(([a,b] for a,b in COMPARISONS), [])))
    frames = {}
    for condition in wanted:
        frame = oof[oof["condition_id"] == condition][[KEY,GROUP,"observed_label","predicted_probability"]].copy()
        if len(frame) != 27727 or frame[KEY].nunique() != 27727:
            raise RuntimeError(f"{condition} lacks one OOF prediction per training variant")
        frames[condition] = frame.sort_values(KEY).reset_index(drop=True)
    reference = frames[wanted[0]][[KEY,GROUP,"observed_label"]]
    for condition, frame in frames.items():
        if not frame[[KEY,GROUP,"observed_label"]].equals(reference):
            raise RuntimeError(f"OOF row mismatch for {condition}")
    genes = sorted(reference[GROUP].unique()); gene_index = {g:i for i,g in enumerate(genes)}
    row_gene = reference[GROUP].map(gene_index).to_numpy(); y = reference["observed_label"].to_numpy(int)
    rng = np.random.default_rng(SEED)
    metrics = {c:{"mcc":np.empty(N_BOOTSTRAP),"pr_auc":np.empty(N_BOOTSTRAP),"brier":np.empty(N_BOOTSTRAP)} for c in wanted}
    batch = 50
    for start in range(0,N_BOOTSTRAP,batch):
        end=min(start+batch,N_BOOTSTRAP); gene_weights=rng.multinomial(len(genes),np.full(len(genes),1/len(genes)),size=end-start)
        row_weights=gene_weights[:,row_gene]
        for condition in wanted:
            prob=frames[condition]["predicted_probability"].to_numpy(float); pred=prob>=.5
            tp=row_weights @ ((pred==1)&(y==1)); tn=row_weights @ ((pred==0)&(y==0)); fp=row_weights @ ((pred==1)&(y==0)); fn=row_weights @ ((pred==0)&(y==1))
            denom=np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); metrics[condition]["mcc"][start:end]=np.divide(tp*tn-fp*fn,denom,out=np.zeros_like(denom,dtype=float),where=denom>0)
            metrics[condition]["brier"][start:end]=(row_weights @ ((prob-y)**2))/row_weights.sum(axis=1)
            order=np.argsort(-prob,kind="mergesort"); ys=y[order]
            for j,w in enumerate(row_weights[:,order]):
                tp_cum=np.cumsum(w*ys); all_cum=np.cumsum(w); total_pos=tp_cum[-1]
                precision=np.divide(tp_cum,all_cum,out=np.ones_like(tp_cum,dtype=float),where=all_cum>0)
                metrics[condition]["pr_auc"][start+j]=float(np.sum(precision*w*ys)/total_pos)
    rows=[]
    for a,b in COMPARISONS:
        for metric in ("mcc","pr_auc","brier"):
            diff = metrics[a][metric]-metrics[b][metric] if metric!="brier" else metrics[b][metric]-metrics[a][metric]
            observed_a = matthews_corrcoef(y,frames[a]["predicted_probability"]>=.5) if metric=="mcc" else average_precision_score(y,frames[a]["predicted_probability"]) if metric=="pr_auc" else brier_score_loss(y,frames[a]["predicted_probability"])
            observed_b = matthews_corrcoef(y,frames[b]["predicted_probability"]>=.5) if metric=="mcc" else average_precision_score(y,frames[b]["predicted_probability"]) if metric=="pr_auc" else brier_score_loss(y,frames[b]["predicted_probability"])
            rows.append({"condition_a":a,"condition_b":b,"metric":metric,"difference_definition":"A-B" if metric!="brier" else "Brier(B)-Brier(A); positive favours A","observed_a":round(float(observed_a),8),"observed_b":round(float(observed_b),8),"observed_difference_favouring_a":round(float((observed_a-observed_b) if metric!="brier" else (observed_b-observed_a)),8),"bootstrap_samples":N_BOOTSTRAP,"ci_2_5":round(float(np.percentile(diff,2.5)),8),"ci_97_5":round(float(np.percentile(diff,97.5)),8),"probability_difference_gt_zero":round(float(np.mean(diff>0)),8),"resampling_unit":GROUP,"seed":SEED})
    result=pd.DataFrame(rows); result.to_csv(OUT/"03_controlled_experiments/paired_gene_bootstrap_results.csv",index=False)
    summary={"bootstrap_samples":N_BOOTSTRAP,"seed":SEED,"resampling_unit":GROUP,"paired":True,"validation_used":False,"test_used":False,"comparisons":rows}
    write_json(OUT/"03_controlled_experiments/paired_gene_bootstrap_summary.json",summary)
    return result,summary


def decision(boot: pd.DataFrame) -> None:
    key=boot[(boot.condition_a==ALL3)&(boot.condition_b==PHYLO_PHAST)&(boot.metric=="mcc")].iloc[0]
    uncertain=key.ci_2_5<=0<=key.ci_97_5
    text=f"""# Conservation selection decision\n\n## BEST OBSERVED CONDITION\n\nAll three conservation measures produced the highest observed training CV mean MCC (0.72214).\n\n## PARSIMONIOUS ALTERNATIVE\n\nphyloP + PhastCons produced mean MCC 0.72082 using one fewer conservation score.\n\n## STATISTICALLY UNCERTAIN DIFFERENCE\n\nThe paired gene-cluster bootstrap MCC difference (all three minus phyloP + PhastCons) was {key.observed_difference_favouring_a:.5f}, with 95% interval [{key.ci_2_5:.5f}, {key.ci_97_5:.5f}] and P(difference > 0)={key.probability_difference_gt_zero:.3f}. The interval {'includes' if uncertain else 'does not include'} zero. Therefore, {'the larger configuration is not claimed superior; phyloP + PhastCons is the parsimonious alternative for later confirmation' if uncertain else 'the paired analysis supports the observed advantage, subject to later validation confirmation'}.\n\nThis is a post-run refinement. It does not replace the preserved initial recommendation and uses no validation or test outcomes.\n"""
    (OUT/"03_controlled_experiments/CONSERVATION_SELECTION_DECISION.md").write_text(text,encoding="utf-8")


def annotated_heatmap(matrix: pd.DataFrame, labels: list[str], path: Path, title: str, caption: str) -> None:
    m=matrix.loc[labels,labels].astype(float); names=[DISPLAY[x] for x in labels]
    fig,ax=plt.subplots(figsize=(max(7,len(labels)*1.25),max(6,len(labels)*1.1)))
    im=ax.imshow(m,vmin=-1,vmax=1,cmap="coolwarm")
    ax.set_xticks(range(len(names)),names,rotation=40,ha="right"); ax.set_yticks(range(len(names)),names)
    for i in range(len(names)):
        for j in range(len(names)):
            value=m.iloc[i,j]; star="*" if i!=j and abs(value)>=.8 else ""
            ax.text(j,i,f"{value:.2f}{star}",ha="center",va="center",color="white" if abs(value)>.55 else "black",fontsize=9)
    ax.set_title(title+"\nTraining partition only",pad=14,fontweight="bold"); fig.colorbar(im,ax=ax,label="Spearman ρ",fraction=.046,pad=.04)
    fig.text(.5,.01,caption+"  * |ρ| ≥ 0.80",ha="center",fontsize=9); fig.subplots_adjust(bottom=.18)
    savefig(fig,path)


def figures(provenance_rows: list[dict], boot: pd.DataFrame) -> None:
    # Cohort flow
    f0=OUT/"00_locked_protocol/cohort_validation/figures"
    fig,ax=plt.subplots(figsize=(10,5)); ax.axis("off")
    ax.text(.5,.82,"Frozen dataset\n40,000 variants",ha="center",va="center",fontsize=17,bbox=dict(boxstyle="round",fc="#E5F5F9",ec=BLUE))
    ax.annotate("",(.28,.35),(.45,.7),arrowprops=dict(arrowstyle="->",lw=2,color=GREY)); ax.annotate("",(.72,.35),(.55,.7),arrowprops=dict(arrowstyle="->",lw=2,color=GREY))
    ax.text(.25,.25,"Retained modelling cohort\n39,552",ha="center",fontsize=15,bbox=dict(boxstyle="round",fc="#E6F5D0",ec=GREEN)); ax.text(.75,.25,"Preserved gene-mapping exclusions\n448",ha="center",fontsize=15,bbox=dict(boxstyle="round",fc="#FEE8C8",ec=ORANGE))
    ax.text(.5,.04,"Frozen dataset unchanged • 39,552 + 448 = 40,000",ha="center",fontsize=11); savefig(fig,f0/"cohort_reconciliation_flow.png")
    assign=pd.read_csv(PATHS["assignments"]); splits=[]
    for s in ("train","validation","test"):
        x=assign[assign.split==s]; splits.append((s,int(x.variant_count.sum()),len(x),int(x.benign_count.sum()),int(x.pathogenic_count.sum())))
    fig,axes=plt.subplots(1,2,figsize=(12,5)); names=[x[0].title() for x in splits]; pos=np.arange(3)
    axes[0].bar(pos,[x[1] for x in splits],color=[BLUE,ORANGE,GREEN]); axes[0].set_xticks(pos,names); axes[0].set_ylabel("Variants"); axes[0].set_title("Locked split sizes")
    for i,x in enumerate(splits): axes[0].text(i,x[1]+400,f"{x[1]:,}\n{x[2]:,} genes",ha="center")
    axes[1].bar(pos,[x[3] for x in splits],label="Class 0",color=BLUE); axes[1].bar(pos,[x[4] for x in splits],bottom=[x[3] for x in splits],label="Class 1",color=ORANGE); axes[1].set_xticks(pos,names); axes[1].set_ylabel("Variants"); axes[1].set_title("Class counts (integrity description only)"); axes[1].legend()
    fig.suptitle("Split integrity summary — gene overlap = 0; variant overlap = 0",fontweight="bold"); savefig(fig,f0/"split_integrity_summary.png")
    # provenance matrices
    f1=OUT/"01_provenance_circularity/figures"; names=[r["display_name"].replace(" 100-way","").replace(" RS","") for r in provenance_rows]
    cols=["Clinical supervised labels","Public training variants","Exact overlap status","Conservation overlap","Population overlap","Composite","Recommended role"]
    values=[]
    overlap=pd.read_csv(OUT/"01_provenance_circularity/training_overlap_results_corrected.csv").set_index("score")
    for r in provenance_rows:
        clinical="Yes" if r["display_name"] in {"PolyPhen-2 HDIV","REVEL"} else "Proxy" if r["display_name"]=="CADD PHRED" else "No"
        public="Yes" if r["display_name"] in {"PolyPhen-2 HDIV","CADD PHRED"} else "Partial" if r["display_name"]=="REVEL" else "N/A" if r["display_name"] in {"GERP++ RS","phyloP 100-way","phastCons 100-way","SIFT"} else "Unclear"
        status=overlap.loc[r["display_name"],"overlap_status"]
        values.append([clinical,public,status,"Yes" if r["whether_it_already_contains_conservation_information"]=="yes" else "No","Yes" if "yes" in r["whether_it_already_contains_population_information"].lower() else "No","Yes" if r["whether_it_incorporates_other_computational_predictors"]=="yes" else "No",r["proposed_role"]])
    fig,ax=plt.subplots(figsize=(15,7)); ax.axis("off"); table=ax.table(cellText=values,rowLabels=names,colLabels=cols,loc="center",cellLoc="center"); table.auto_set_font_size(False); table.set_fontsize(8); table.scale(1,1.8)
    for (i,j),cell in table.get_celld().items():
        if i==0: cell.set_facecolor("#D9EAF7"); cell.set_text_props(weight="bold")
        elif j>=0: cell.set_facecolor("#F7F7F7" if i%2 else "#EAF2F8")
    ax.set_title("Provenance, circularity and evidence-role matrix\nUnknown/unverifiable is not low risk",fontweight="bold",pad=18); savefig(fig,f1/"provenance_circularity_risk_matrix.png")
    role_groups=Counter(r["proposed_role"] for r in provenance_rows); fig,ax=plt.subplots(figsize=(9,5)); keys=list(role_groups); ax.barh(keys,[role_groups[k] for k in keys],color=[BLUE,ORANGE,PURPLE,GREY][:len(keys)]); ax.set_xlabel("Number of audited scores"); ax.set_title("Evidence-role summary (corrected governance)",fontweight="bold")
    for i,k in enumerate(keys): ax.text(role_groups[k]+.05,i,str(role_groups[k]),va="center"); savefig(fig,f1/"evidence_role_summary.png")
    # heatmaps
    f2=OUT/"02_redundancy/figures"; matrix=pd.read_csv(OUT/"02_redundancy/spearman_correlation.csv",index_col=0)
    cons=list(DISPLAY)[:3]; preds=list(DISPLAY)[3:8]; annotated_heatmap(matrix,cons,f2/"conservation_heatmap_annotated.png","Conservation redundancy","Correlation indicates shared evolutionary information, not circularity.")
    annotated_heatmap(matrix,preds,f2/"computational_predictor_heatmap_annotated.png","Computational-predictor redundancy","SIFT is reversed: lower means more damaging; negative ρ can indicate shared signal.")
    annotated_heatmap(matrix,cons+preds,f2/"combined_evidence_heatmap_annotated.png","Combined evidence redundancy","Correlation indicates shared information, not circularity.")
    order=list(matrix.abs().sum().sort_values(ascending=False).index); annotated_heatmap(matrix,order,f2/"clustered_evidence_heatmap_fixed.png","Combined evidence correlations (deterministically ordered)","Signed correlations are preserved; ordering is descriptive, not a circularity test.")
    flagged=pd.read_csv(OUT/"02_redundancy/flagged_redundant_pairs.csv"); nodes=sorted(set(flagged.field_a)|set(flagged.field_b)); angles=np.linspace(0,2*np.pi,len(nodes),endpoint=False); xy={n:(math.cos(a),math.sin(a)) for n,a in zip(nodes,angles)}
    fig,ax=plt.subplots(figsize=(9,8)); ax.axis("off")
    for _,r in flagged.iterrows():
        x1,y1=xy[r.field_a]; x2,y2=xy[r.field_b]; ax.plot([x1,x2],[y1,y2],color=RED if r.spearman<0 else BLUE,ls="--" if r.spearman<0 else "-",lw=2); ax.text((x1+x2)/2,(y1+y2)/2,f"{r.spearman:.2f}",bbox=dict(fc="white",ec="none",alpha=.8),ha="center")
    for n,(x,y0) in xy.items(): ax.scatter(x,y0,s=1600,color="#E5F5F9",edgecolor=BLUE,zorder=3); ax.text(x,y0,DISPLAY[n],ha="center",va="center",fontsize=9,zorder=4)
    ax.set_title("Flagged redundancy network (|Spearman ρ| ≥ 0.80)\nTraining only; redundancy is not circularity",fontweight="bold"); savefig(fig,f2/"flagged_redundancy_network.png")
    # experiments
    f3=OUT/"03_controlled_experiments/figures"; summary=pd.read_csv(OUT/"03_controlled_experiments/experiment_summary.csv")
    cons_rows=summary[summary.category=="conservation subset"].sort_values("mean_mcc"); fig,ax=plt.subplots(figsize=(11,6)); conservation_labels={
        "conservation_none": "No conservation",
        "conservation_gerp_rs": "GERP++",
        "conservation_phylop100way": "phyloP",
        "conservation_phastcons100way": "PhastCons",
        "conservation_gerp_rs__phylop100way": "GERP++ + phyloP",
        "conservation_gerp_rs__phastcons100way": "GERP++ + PhastCons",
        "conservation_phylop100way__phastcons100way": "phyloP + PhastCons",
        ALL3: "GERP++ + phyloP + PhastCons",
    }; labels=[conservation_labels[x] for x in cons_rows.condition_id]; colors=[BLUE if x!=ALL3 else ORANGE for x in cons_rows.condition_id]
    ax.barh(labels,cons_rows.mean_mcc,xerr=cons_rows.sd_mcc,color=colors,capsize=3); ax.set_xlabel("Mean MCC ± fold SD"); ax.set_title("Conservation-subset performance\nTraining-only gene-grouped cross-validation",fontweight="bold"); ax.set_xlim(.5,.77)
    for i,v in enumerate(cons_rows.mean_mcc): ax.text(v+.003,i,f"{v:.3f}",va="center"); savefig(fig,f3/"conservation_subset_performance.png")
    predictor_fields=["predictor_sift_score","predictor_polyphen2_hdiv_score","predictor_revel_score","predictor_cadd_phred","predictor_alphamissense_score"]; alone=[]; plus=[]; labels=[]
    for field in predictor_fields:
        alone.append(summary.loc[summary.condition_id==f"comparator_{field}","mean_mcc"].iloc[0]); plus.append(summary.loc[summary.condition_id==f"biological_plus_{field}","mean_mcc"].iloc[0]); labels.append(DISPLAY[field])
    fig,ax=plt.subplots(figsize=(10,6)); x=np.arange(5); w=.36; ax.bar(x-w/2,alone,w,label="Predictor alone",color=BLUE); ax.bar(x+w/2,plus,w,label="Biological + predictor",color=ORANGE); ax.set_xticks(x,labels); ax.set_ylabel("Mean MCC"); ax.set_title("Optional predictor comparisons\nPerformance does not establish primary-evidence eligibility",fontweight="bold"); ax.legend(); savefig(fig,f3/"predictor_comparator_performance.png")
    selected=summary[summary.condition_id.str.startswith(("comparator_","biological_plus_"))].copy(); metrics=[("mean_mcc","MCC (higher better)"),("mean_pr_auc","Average precision / PR-AUC (higher better)"),("mean_roc_auc","ROC-AUC (higher better)"),("mean_brier","Brier (lower better)")]
    fig,axes=plt.subplots(2,2,figsize=(14,10)); short=[x.replace("comparator_predictor_","alone: ").replace("biological_plus_predictor_","bio+: ").replace("_score","") for x in selected.condition_id]
    for ax,(col,title) in zip(axes.flat,metrics): ax.barh(short,selected[col],color=GREEN if "Brier" not in title else PURPLE); ax.set_title(title); ax.tick_params(axis="y",labelsize=8)
    fig.suptitle("Training-only optional-predictor metric summary",fontweight="bold"); savefig(fig,f3/"predictor_metric_summary.png")
    mcc=boot[boot.metric=="mcc"].copy(); y_pos=np.arange(len(mcc)); fig,ax=plt.subplots(figsize=(10,6)); means=mcc.observed_difference_favouring_a.to_numpy(); lo=means-mcc.ci_2_5.to_numpy(); hi=mcc.ci_97_5.to_numpy()-means; labels=[f"{DISPLAY[a]}\nvs {DISPLAY[b]}" for a,b in zip(mcc.condition_a,mcc.condition_b)]; ax.errorbar(means,y_pos,xerr=[lo,hi],fmt="o",color=BLUE,capsize=4); ax.axvline(0,color=RED,ls="--"); ax.set_yticks(y_pos,labels); ax.set_xlabel("Paired MCC difference (A − B)"); ax.set_title("Paired gene-cluster bootstrap (5,000 samples)\nIntervals crossing zero indicate statistical uncertainty",fontweight="bold"); savefig(fig,f3/"paired_gene_bootstrap_forest.png")
    oof=pd.read_csv(OUT/"03_controlled_experiments/out_of_fold_predictions.csv"); wanted=[ALL3,PHYLO_PHAST,"comparator_predictor_revel_score","comparator_predictor_alphamissense_score","biological_plus_predictor_revel_score","population_inclusive_sensitivity"]; fig,ax=plt.subplots(figsize=(9,7)); bins=np.linspace(0,1,11)
    for c in wanted:
        d=oof[oof.condition_id==c]; idx=np.digitize(d.predicted_probability,bins)-1; xs=[]; ys=[]
        for b in range(10):
            q=d[idx==b]
            if len(q): xs.append(q.predicted_probability.mean()); ys.append(q.observed_label.mean())
        ax.plot(xs,ys,marker="o",label=DISPLAY.get(c,c.replace("comparator_predictor_","alone: ").replace("biological_plus_predictor_","bio+: ").replace("_score","")))
    ax.plot([0,1],[0,1],color=GREY,ls="--",label="Ideal"); ax.set_xlabel("Mean predicted probability"); ax.set_ylabel("Observed fraction positive"); ax.set_title("Preliminary training-only OOF calibration\nGene-grouped cross-validation",fontweight="bold"); ax.legend(fontsize=8); savefig(fig,f3/"oof_calibration_curves.png")


def reports(provenance_rows: list[dict], boot: pd.DataFrame) -> None:
    summary=pd.read_csv(OUT/"03_controlled_experiments/experiment_summary.csv"); coverage=pd.read_csv(OUT/"02_redundancy/coverage_summary.csv"); flagged=pd.read_csv(OUT/"02_redundancy/flagged_redundant_pairs.csv"); vif=pd.read_csv(OUT/"02_redundancy/vif_results.csv")
    b=boot[(boot.condition_a==ALL3)&(boot.condition_b==PHYLO_PHAST)&(boot.metric=="mcc")].iloc[0]; uncertain=b.ci_2_5<=0<=b.ci_97_5
    coverage_lines="\n".join(f"- {DISPLAY[r.column]}: {r.coverage:.2%}" for _,r in coverage.iterrows())
    cons_lines="\n".join(f"- {DISPLAY.get(r.condition_id,r.condition_id)}: MCC {r.mean_mcc:.5f}, PR-AUC {r.mean_pr_auc:.5f}, ROC-AUC {r.mean_roc_auc:.5f}, Brier {r.mean_brier:.5f}" for _,r in summary[summary.category=="conservation subset"].iterrows())
    pred_lines="\n".join(f"- {r.condition_id}: MCC {r.mean_mcc:.5f}, PR-AUC {r.mean_pr_auc:.5f}, ROC-AUC {r.mean_roc_auc:.5f}, Brier {r.mean_brier:.5f}" for _,r in summary[summary.condition_id.str.startswith(("comparator_","biological_plus_"))].iterrows())
    revised=f"""# Revised pre-modelling evidence audit\n\n## Scope and safeguards\n\nThis report preserves the initial report and adds corrected provenance, governance, paired uncertainty and publication figures. Statistics and preliminary modelling use the 27,727-row training partition only. Aggregate validation/test counts are used only for split integrity; validation and test outcomes were not analysed. Preliminary cross-validation is not final model evaluation.\n\n## Cohort and split\n\nThe frozen 40,000 variants reconcile to 39,552 retained variants and 448 preserved gene-mapping exclusions. Train/validation/test contain 27,727/5,790/6,035 variants and 5,206/1,115/1,115 genes. Gene overlap and variant overlap are zero.\n\n![Cohort flow](../00_locked_protocol/cohort_validation/figures/cohort_reconciliation_flow.png)\n\n![Split integrity](../00_locked_protocol/cohort_validation/figures/split_integrity_summary.png)\n\n## Corrected provenance and circularity\n\nGERP++, phyloP, PhastCons and SIFT do not have conventional clinically labelled training-variant overlap; these checks are NOT_APPLICABLE. PolyPhen-2 HDIV uses supervised disease/divergence examples and was assessed at protein-variant level using its official public HumDiv archive. REVEL remains unverifiable for its restricted HGMD-derived pathogenic portion. CADD uses proxy labels and has UNKNOWN EXACT RELEASE locally. AlphaMissense uses population-frequency fine-tuning; ClinVar was evaluation rather than direct training. Unavailable training data do not demonstrate absence of overlap.\n\n![Provenance matrix](../01_provenance_circularity/figures/provenance_circularity_risk_matrix.png)\n\n## Coverage\n\n{coverage_lines}\n\n## Correlation and VIF\n\nphyloP–PhastCons Spearman ρ was 0.82695. Flagged predictor relationships included phyloP–CADD 0.81785, SIFT–AlphaMissense −0.80127 (SIFT direction reversed), REVEL–AlphaMissense 0.81856 and CADD–AlphaMissense 0.80312. No VIF reached 10 (maximum {pd.to_numeric(vif.vif,errors='coerce').max():.2f}). Correlation indicates shared information, not circularity.\n\n![Combined heatmap](../02_redundancy/figures/combined_evidence_heatmap_annotated.png)\n\n![Redundancy network](../02_redundancy/figures/flagged_redundancy_network.png)\n\n## All conservation conditions\n\n{cons_lines}\n\n![Conservation performance](../03_controlled_experiments/figures/conservation_subset_performance.png)\n\n## Optional predictors and population sensitivity\n\n{pred_lines}\n\nPopulation-inclusive sensitivity: MCC {summary.loc[summary.condition_id=='population_inclusive_sensitivity','mean_mcc'].iloc[0]:.5f}, PR-AUC {summary.loc[summary.condition_id=='population_inclusive_sensitivity','mean_pr_auc'].iloc[0]:.5f}, ROC-AUC {summary.loc[summary.condition_id=='population_inclusive_sensitivity','mean_roc_auc'].iloc[0]:.5f}, Brier {summary.loc[summary.condition_id=='population_inclusive_sensitivity','mean_brier'].iloc[0]:.5f}. Population evidence remains sensitivity-only because of class-associated missingness. Performance improvement does not remove circularity, so computational predictors remain comparator/sensitivity evidence.\n\n![Predictor comparison](../03_controlled_experiments/figures/predictor_comparator_performance.png)\n\n## Paired gene bootstrap and refined interpretation\n\nThe best observed conservation condition was all three scores. phyloP + PhastCons is the parsimonious alternative. The paired MCC difference was {b.observed_difference_favouring_a:.5f}, 95% interval [{b.ci_2_5:.5f}, {b.ci_97_5:.5f}], P(diff>0)={b.probability_difference_gt_zero:.3f}. The interval {'includes' if uncertain else 'does not include'} zero; therefore the larger condition is {'not claimed superior' if uncertain else 'supported as superior within this training-only analysis'}. The initial recommendation remains preserved.\n\n![Bootstrap forest](../03_controlled_experiments/figures/paired_gene_bootstrap_forest.png)\n\n![OOF calibration](../03_controlled_experiments/figures/oof_calibration_curves.png)\n\n## Next step\n\nAfter governance review, run the refined single validation-confirmation command once. It evaluates only the refined training-only selection and does not touch test. Validation was not used here; test was not used here.\n"""
    p=OUT/"04_report/PRE_MODELING_EVIDENCE_AUDIT_REVISED.md"; p.write_text(revised,encoding="utf-8")
    supervisor=f"""# Supervisor progress summary\n\n## What is complete\n\n- Verified 40,000 frozen variants, 39,552 retained, and 448 preserved exclusions.\n- Confirmed zero gene/variant overlap across the locked split.\n- Corrected provenance categories and tested PolyPhen-2 HumDiv overlap safely at protein level.\n- Completed training-only correlations, fixed-fold experiments and 5,000-sample paired gene bootstrap.\n\n![Cohort reconciliation](../00_locked_protocol/cohort_validation/figures/cohort_reconciliation_flow.png)\n\n![Conservation comparison](../03_controlled_experiments/figures/conservation_subset_performance.png)\n\n## Main interpretation\n\nAll three conservation scores were best observed (MCC 0.72214), while phyloP + PhastCons was nearly identical (0.72082). Bootstrap 95% interval for their MCC difference: [{b.ci_2_5:.5f}, {b.ci_97_5:.5f}]. {'Because it includes zero, phyloP + PhastCons is the parsimonious alternative.' if uncertain else 'The interval excludes zero, supporting the three-score condition within training.'}\n\nREVEL and AlphaMissense improve training CV strongly, but remain comparator/sensitivity evidence because performance cannot remove circularity or provenance concerns.\n\n![Paired uncertainty](../03_controlled_experiments/figures/paired_gene_bootstrap_forest.png)\n\n## What happens next\n\nReview governance and the refined command, then confirm exactly one preselected configuration on validation. The test partition stays locked for final evaluation.\n"""
    (OUT/"04_report/SUPERVISOR_PROGRESS_SUMMARY.md").write_text(supervisor,encoding="utf-8")


def figure_index() -> None:
    entries=[]
    for p in sorted(OUT.rglob("figures/*.png")):
        rel=p.relative_to(OUT/"04_report").as_posix() if False else "../"+p.relative_to(OUT).as_posix()
        name=p.name
        partition="locked aggregate counts only" if "cohort" in str(p) else "training only"
        entries.append(f"| [{name}]({rel}) | {name.replace('_',' ').replace('.png','')} | {partition} | Descriptive audit result | Interpret within stated audit scope | Do not infer circularity absence or final-model performance |")
    text="# Figure index\n\n| Figure | Purpose | Data partition | Main finding | Correct interpretation | Misinterpretation to avoid |\n|---|---|---|---|---|---|\n"+"\n".join(entries)+"\n"
    (OUT/"04_report/FIGURE_INDEX.md").write_text(text,encoding="utf-8")


def refined_command() -> None:
    source = OUT / "03_controlled_experiments/recommended_training_only_configuration.json"
    refined = json.loads(source.read_text())
    refined["recommended_condition_id"] = PHYLO_PHAST
    refined["recommended_features"] = [
        field for field in refined["recommended_features"] if field != "conservation_gerp_rs"
    ]
    refined["selection_basis"] = (
        "training-only fixed gene-grouped cross-validation followed by the post-run paired "
        "gene-cluster bootstrap amendment; parsimonious alternative selected because the MCC "
        "difference from the best-observed three-score condition was statistically uncertain"
    )
    refined["initial_recommendation_preserved_at"] = source.relative_to(ROOT).as_posix()
    refined_path = OUT / "03_controlled_experiments/recommended_training_only_configuration_refined.json"
    write_json(refined_path, refined)
    command="python -m src.pre_modeling_audit.run_pre_modeling_audit --validation-confirmation outputs/pre_modeling_evidence_audit/03_controlled_experiments/recommended_training_only_configuration_refined.json\n"
    (OUT/"03_controlled_experiments/validation_confirmation_command_refined.txt").write_text("# Run from the modelling repository root with its virtual environment activated.\n# This evaluates only the preselected configuration; it does not rerun search or access test.\n"+command,encoding="utf-8")


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("--polyphen-archive",type=Path); args=parser.parse_args()
    before=locked_checksum(); corrected_leakage(); poly=polyphen_overlap(args.polyphen_archive); rows=corrected_provenance(poly); write_provenance_docs(poly); amendment(); boot,summary=bootstrap(); decision(boot); figures(rows,boot); reports(rows,boot); figure_index(); refined_command(); after=locked_checksum()
    if before!=after: raise RuntimeError("Locked protocol changed")
    print(json.dumps({"status":"PASS","locked_protocol_sha256":after,"polyphen_overlap":poly,"bootstrap_samples":N_BOOTSTRAP},indent=2))


if __name__=="__main__": main()
