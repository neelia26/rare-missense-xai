#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class TailIntegrationTest(unittest.TestCase):
    def test_full_run_dry_run_builds_complete_plan(self):
        with tempfile.TemporaryDirectory() as temporary:
            temporary_path = Path(temporary)
            constraint = temporary_path / "constraint.tsv"
            constraint.write_text("gene\ttranscript\tmis.z_score\tlof.oe_ci.upper\n", encoding="utf-8")
            config = json.loads(
                (ROOT / "config/pipeline_config.json").read_text(encoding="utf-8")
            )
            config["gnomad"]["constraint_file"] = str(constraint)
            config_path = temporary_path / "pipeline_config.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            result = subprocess.run(
                [
                    sys.executable,
                    "run_pipeline.py",
                    "run-full",
                    "--config",
                    str(config_path),
                    "--dry-run",
                ],
                cwd=ROOT,
                check=True,
                capture_output=True,
                text=True,
            )
            self.assertIn("Governed full-run plan", result.stdout)
            self.assertIn("Dry run only", result.stdout)
            self.assertNotIn("08b_enrich_gnomad_counts.py", result.stdout)

    def test_build_split_and_audit(self):
        with tempfile.TemporaryDirectory() as temporary:
            clone = Path(temporary) / "pipeline"
            shutil.copytree(
                ROOT,
                clone,
                ignore=shutil.ignore_patterns(".venv", "__pycache__", "*.pyc"),
            )
            input_path = clone / "synthetic_stage11.csv"
            fields = [
                "variant_key", "binary_label", "strict_missense_eligible",
                "uniprot_validation_status", "reconciled_gene_symbol",
                "gnomad_rarity_status", "gnomad_global_af", "gnomad_max_population_af",
                "conservation_gerp_rs", "conservation_phylop100way",
                "conservation_phastcons100way", "gnomad_missense_z", "gnomad_loeuf",
                "aa_grantham_score", "aa_composition_delta", "aa_polarity_delta",
                "aa_volume_delta", "aa_hydropathy_delta", "aa_charge_delta",
                "aa_charge_changed", "aa_polarity_class_changed",
                "aa_relative_protein_position", "interpro_overlap_count",
                "interpro_overlap_types", "uniprot_overlapping_feature_types",
            ]
            with input_path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                for gene_number in range(20):
                    for label in (0, 1):
                        writer.writerow(
                            {
                                "variant_key": f"GRCh38:1:{1000 + gene_number * 2 + label}:A:G",
                                "binary_label": label,
                                "strict_missense_eligible": 1,
                                "uniprot_validation_status": "MATCH",
                                "reconciled_gene_symbol": f"GENE{gene_number:02d}",
                                "gnomad_rarity_status": "RARE",
                                "gnomad_global_af": 0.0001 + gene_number / 1000000,
                                "gnomad_max_population_af": 0.0002 + gene_number / 1000000,
                                "conservation_gerp_rs": gene_number / 10,
                                "conservation_phylop100way": gene_number / 20,
                                "conservation_phastcons100way": gene_number / 100,
                                "gnomad_missense_z": gene_number / 5,
                                "gnomad_loeuf": 0.2 + gene_number / 100,
                                "aa_grantham_score": 50 + gene_number,
                                "aa_composition_delta": 0.1,
                                "aa_polarity_delta": 0.2,
                                "aa_volume_delta": 1.0,
                                "aa_hydropathy_delta": -0.5,
                                "aa_charge_delta": label,
                                "aa_charge_changed": label,
                                "aa_polarity_class_changed": gene_number % 2,
                                "aa_relative_protein_position": 0.1 + gene_number / 100,
                                "interpro_overlap_count": gene_number % 3,
                                "interpro_overlap_types": "domain" if gene_number % 2 else "",
                                "uniprot_overlapping_feature_types": "Region",
                            }
                        )
            subprocess.run(
                [
                    sys.executable,
                    "scripts/12_build_final_dataset.py",
                    "--input",
                    str(input_path),
                    "--target-per-class",
                    "15",
                    "--allow-smaller",
                ],
                cwd=clone,
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sys.executable, "scripts/13_split_by_gene.py"],
                cwd=clone,
                check=True,
                capture_output=True,
                text=True,
            )
            subprocess.run(
                [sys.executable, "scripts/14_audit_final_dataset.py"],
                cwd=clone,
                check=True,
                capture_output=True,
                text=True,
            )
            report = json.loads(
                (clone / "audits/14_final_dataset_audit.json").read_text(encoding="utf-8")
            )
            self.assertEqual(report["status"], "PASS")
            self.assertFalse(any(report["split_gene_overlap"].values()))
            self.assertEqual(report["dataset"]["class_counts"], {"0": 15, "1": 15})


if __name__ == "__main__":
    unittest.main()
