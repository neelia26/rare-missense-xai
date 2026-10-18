#!/usr/bin/env python3
"""List pipeline status or run one numbered stage with transparent arguments."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STAGES = {
    "00": ("scripts/00_check_environment.py", "audits/00_environment_report.json"),
    "01": ("scripts/01_fetch_clinvar.py", "data/raw/clinvar_grch38.vcf.gz"),
    "02": ("scripts/02_extract_candidates.py", "data/interim/02_clinvar_candidates.csv"),
    "03": ("scripts/03_select_candidate_pool.py", "data/interim/03_candidate_pool.csv"),
    "04": ("scripts/04_prepare_vep_input.py", "data/raw/vep/04_candidates.vcf"),
    "05": ("scripts/05_run_vep.py", "data/raw/vep/05_vep.jsonl"),
    "06": ("scripts/06_parse_vep.py", "data/interim/06_vep_reconciled.csv"),
    "07": ("scripts/07_enrich_uniprot.py", "data/interim/07_uniprot_validated.csv"),
    "08": ("scripts/08_enrich_interpro.py", "data/interim/08_interpro.csv"),
    "08b": ("scripts/08b_enrich_gnomad_counts.py", "data/interim/08b_gnomad_counts.csv"),
    "09": ("scripts/09_merge_gnomad_constraint.py", "data/interim/09_constraint.csv"),
    "10": ("scripts/10_enrich_dbnsfp.py", "data/interim/10_dbnsfp.csv"),
    "11": ("scripts/11_add_amino_acid_features.py", "data/interim/11_amino_acid_features.csv"),
    "12": ("scripts/12_build_final_dataset.py", "data/processed/final_dataset_full.csv"),
    "13": ("scripts/13_split_by_gene.py", "data/processed/split_gene_assignments.csv"),
    "14": ("scripts/14_audit_final_dataset.py", "audits/14_final_dataset_audit.json"),
}


def show_status() -> None:
    print("Stage  Status   Expected output")
    for stage, (_, output) in STAGES.items():
        exists = (ROOT / output).exists()
        print(f"{stage:>4}   {'DONE' if exists else 'PENDING':<7}  {output}")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("status")
    run = subparsers.add_parser("run")
    run.add_argument("stage", choices=list(STAGES))
    run.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    if args.action == "status":
        show_status()
        return
    script = ROOT / STAGES[args.stage][0]
    extra = args.arguments
    if extra[:1] == ["--"]:
        extra = extra[1:]
    command = [sys.executable, str(script), *extra]
    print("Running:", " ".join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()

