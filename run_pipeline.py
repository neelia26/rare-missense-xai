#!/usr/bin/env python3
"""Inspect, run one stage, or execute the governed full dataset rebuild."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))

from src.pipeline_utils import file_record, load_config, sha256_file, utc_now, write_json


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


def project_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def show_status() -> None:
    print("Stage  Status   Expected output")
    for stage, (_, output) in STAGES.items():
        exists = (ROOT / output).exists()
        print(f"{stage:>4}   {'DONE' if exists else 'PENDING':<7}  {output}")
    default_state = ROOT / "logs/full_pipeline_state.json"
    if default_state.exists():
        state = json.loads(default_state.read_text(encoding="utf-8"))
        print(f"\nFull-run state: {state.get('status', 'UNKNOWN')} ({default_state})")


def validate_full_config(config: dict[str, Any]) -> Path:
    if config.get("vep", {}).get("mode") != "rest":
        raise ValueError("The automated API workflow requires vep.mode='rest'")
    sampling = config.get("sampling", {})
    final = int(sampling.get("final_per_class", 0))
    if final < 1:
        raise ValueError("sampling.final_per_class must be positive")
    if int(sampling.get("candidate_benign", 0)) < final:
        raise ValueError("candidate_benign cannot be smaller than final_per_class")
    if int(sampling.get("candidate_pathogenic", 0)) < final:
        raise ValueError("candidate_pathogenic cannot be smaller than final_per_class")
    expected_dbnsfp = str(config.get("dbnsfp", {}).get("expected_release", ""))
    if not expected_dbnsfp or expected_dbnsfp == "set_after_downloading":
        raise ValueError("dbnsfp.expected_release must be an explicit API-reported release")
    gnomad = config.get("gnomad", {})
    if not gnomad.get("source_release"):
        raise ValueError("gnomad.source_release is required")
    constraint_path = project_path(str(gnomad.get("constraint_file", "")))
    if not constraint_path.is_file():
        raise FileNotFoundError(
            f"Configured gnomAD constraint file does not exist: {constraint_path}"
        )
    if gnomad.get("direct_counts_enabled"):
        if not gnomad.get("exomes_vcf_template") and not gnomad.get("genomes_vcf_template"):
            raise ValueError(
                "Enable direct counts only with an exome or genome VCF template"
            )
    return constraint_path


def full_stage_arguments(
    config_path: Path, config: dict[str, Any], constraint_path: Path
) -> list[tuple[str, list[str]]]:
    cfg = str(config_path)
    plan: list[tuple[str, list[str]]] = [
        ("00", ["--config", cfg]),
        ("01", ["--config", cfg]),
        ("02", ["--config", cfg]),
        ("03", ["--config", cfg]),
        ("04", []),
        ("05", ["--mode", "rest", "--config", cfg]),
        ("06", ["--config", cfg]),
        ("07", ["--config", cfg]),
        ("08", ["--config", cfg]),
    ]
    gnomad = config["gnomad"]
    stage09_input = ROOT / "data/interim/08_interpro.csv"
    if gnomad.get("direct_counts_enabled"):
        arguments = [
            "--source-release",
            str(gnomad["source_release"]),
            "--rarity-threshold",
            str(config["rarity"]["maximum_global_af"]),
        ]
        if gnomad.get("exomes_vcf_template"):
            arguments.extend(["--exomes-vcf", str(gnomad["exomes_vcf_template"])])
        if gnomad.get("genomes_vcf_template"):
            arguments.extend(["--genomes-vcf", str(gnomad["genomes_vcf_template"])])
        plan.append(("08b", arguments))
        stage09_input = ROOT / "data/interim/08b_gnomad_counts.csv"
    plan.extend(
        [
            (
                "09",
                [
                    "--input",
                    str(stage09_input),
                    "--constraint-file",
                    str(constraint_path),
                    "--source-release",
                    str(gnomad["source_release"]),
                ],
            ),
            ("10", ["--config", cfg]),
            ("11", []),
            ("12", ["--config", cfg]),
            ("13", ["--config", cfg]),
            ("14", []),
        ]
    )
    return plan


def frozen_clinvar_record() -> dict[str, Any] | None:
    source = ROOT / "data/raw/clinvar_grch38.vcf.gz"
    provenance_path = ROOT / "data/raw/clinvar_grch38.provenance.json"
    if not source.exists() and not provenance_path.exists():
        return None
    if not source.exists() or not provenance_path.exists():
        raise RuntimeError(
            "ClinVar frozen file and provenance JSON must either both exist or both be absent"
        )
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    expected = str(provenance.get("file", {}).get("sha256", ""))
    actual = sha256_file(source)
    if not expected or actual != expected:
        raise RuntimeError(
            "Existing ClinVar VCF does not match its frozen provenance SHA-256"
        )
    return {
        "reused_frozen_source": True,
        "output": file_record(source),
        "provenance": file_record(provenance_path),
    }


def state_entry_is_valid(entry: dict[str, Any], expected_output: Path) -> bool:
    if not expected_output.exists():
        return False
    expected_hash = str(entry.get("output", {}).get("sha256", ""))
    return bool(expected_hash) and sha256_file(expected_output) == expected_hash


def run_full(config_value: str, resume: bool, dry_run: bool) -> None:
    config_path = project_path(config_value)
    config = load_config(config_path)
    constraint_path = validate_full_config(config)
    plan = full_stage_arguments(config_path, config, constraint_path)
    state_path = project_path(
        str(config.get("automation", {}).get("state_file", "logs/full_pipeline_state.json"))
    )
    config_sha256 = sha256_file(config_path)
    frozen = frozen_clinvar_record()

    print("Governed full-run plan")
    print(f"Config: {config_path}")
    print(
        "Candidates: "
        f"{int(config['sampling']['candidate_benign']):,} benign-class + "
        f"{int(config['sampling']['candidate_pathogenic']):,} pathogenic-class"
    )
    print(f"Final target: {int(config['sampling']['final_per_class']):,} per class")
    print(
        "Direct gnomAD counts: "
        + ("enabled" if config["gnomad"].get("direct_counts_enabled") else "disabled (optional)")
    )
    for stage, arguments in plan:
        if stage == "01" and frozen is not None:
            print("  01  REUSE frozen ClinVar file after SHA-256 verification")
            continue
        script = ROOT / STAGES[stage][0]
        print(f"  {stage:>3}  {shlex.join([sys.executable, str(script), *arguments])}")
    if dry_run:
        print("Dry run only: nothing was executed or overwritten.")
        return

    if state_path.exists():
        if not resume:
            raise SystemExit(
                f"A full-run state already exists at {state_path}. Use --resume to continue it."
            )
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("config_sha256") != config_sha256:
            raise RuntimeError(
                "The configuration changed after this full run started. Archive the state and "
                "begin a reviewed new run instead of mixing configurations."
            )
    else:
        state = {
            "schema_version": 1,
            "status": "RUNNING",
            "created_at_utc": utc_now(),
            "updated_at_utc": utc_now(),
            "config": file_record(config_path),
            "config_sha256": config_sha256,
            "completed_stages": {},
        }
        write_json(state_path, state)

    downstream_dirty = False
    for stage, arguments in plan:
        expected_output = ROOT / STAGES[stage][1]
        previous = state["completed_stages"].get(stage)
        if stage == "01" and frozen is not None:
            state["completed_stages"][stage] = {
                "completed_at_utc": utc_now(),
                "mode": "reused_frozen_source",
                **frozen,
            }
            state["updated_at_utc"] = utc_now()
            write_json(state_path, state)
            print("Stage 01: verified and reused the frozen ClinVar release")
            continue
        if (
            resume
            and not downstream_dirty
            and isinstance(previous, dict)
            and state_entry_is_valid(previous, expected_output)
        ):
            print(f"Stage {stage}: verified previous output; skipping")
            continue

        downstream_dirty = True
        script = ROOT / STAGES[stage][0]
        command = [sys.executable, str(script), *arguments]
        print(f"\nRunning stage {stage}: {shlex.join(command)}", flush=True)
        subprocess.run(command, cwd=ROOT, check=True)
        if not expected_output.exists():
            raise RuntimeError(
                f"Stage {stage} completed without its expected output: {expected_output}"
            )
        state["completed_stages"][stage] = {
            "completed_at_utc": utc_now(),
            "command": command,
            "output": file_record(expected_output),
        }
        state["updated_at_utc"] = utc_now()
        write_json(state_path, state)

    state["status"] = "COMPLETE"
    state["completed_at_utc"] = utc_now()
    state["updated_at_utc"] = utc_now()
    write_json(state_path, state)
    print(f"\nFull dataset pipeline completed. State: {state_path}")


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="action", required=True)
    subparsers.add_parser("status")

    run = subparsers.add_parser("run")
    run.add_argument("stage", choices=list(STAGES))
    run.add_argument("arguments", nargs=argparse.REMAINDER)

    run_full_parser = subparsers.add_parser(
        "run-full", help="Execute the governed full rebuild with resumable state"
    )
    run_full_parser.add_argument(
        "--config", default="config/pipeline_config.json"
    )
    run_full_parser.add_argument(
        "--resume",
        action="store_true",
        help="Resume verified stages from logs/full_pipeline_state.json",
    )
    run_full_parser.add_argument(
        "--dry-run", action="store_true", help="Print and validate the complete plan only"
    )

    args = parser.parse_args()
    if args.action == "status":
        show_status()
        return
    if args.action == "run-full":
        run_full(args.config, args.resume, args.dry_run)
        return

    script = ROOT / STAGES[args.stage][0]
    extra = args.arguments
    if extra[:1] == ["--"]:
        extra = extra[1:]
    command = [sys.executable, str(script), *extra]
    print("Running:", shlex.join(command), flush=True)
    subprocess.run(command, cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
