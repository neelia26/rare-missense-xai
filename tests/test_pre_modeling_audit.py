from __future__ import annotations

import tempfile
import unittest
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from src.pre_modeling_audit import run_pre_modeling_audit as audit


def load_preparation_script():
    path = audit.ROOT / "scripts/12b_prepare_modeling_dataset.py"
    spec = importlib.util.spec_from_file_location("prepare_modeling_dataset", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class PreModelingAuditTests(unittest.TestCase):
    def test_registry_generator_keeps_population_sensitivity_only(self):
        module = load_preparation_script()
        rows = module.registry_rows(
            ["gnomad_global_af", "conservation_gerp_rs"],
            {"conservation_gerp_rs"},
        )
        by_field = {row["field"]: row for row in rows}
        self.assertEqual(by_field["gnomad_global_af"]["primary_model_feature"], 0)
        self.assertEqual(by_field["conservation_gerp_rs"]["primary_model_feature"], 1)

    def test_variant_accounting_matches_locked_cohort(self):
        cohort, report = audit.validate_cohort()
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["row_counts"], {"frozen": 40000, "cohort": 39552, "excluded": 448})
        self.assertEqual(len(cohort), 39552)

    def test_existing_split_is_disjoint_and_complete(self):
        cohort = audit.load(audit.PATHS["cohort"])
        report = audit.validate_split(cohort)
        self.assertEqual(report["status"], "PASS")
        self.assertTrue(report["split_union_equals_cohort"])
        self.assertFalse(any(report["gene_overlap"].values()))
        self.assertFalse(any(report["variant_overlap"].values()))

    def test_missing_values_are_not_zero(self):
        values = pd.Series(["", ".", "NA", "0", "1.5"])
        parsed = audit.numeric(values)
        self.assertTrue(parsed.iloc[:3].isna().all())
        self.assertEqual(parsed.iloc[3], 0)

    def test_fixed_gene_folds_do_not_overlap(self):
        data = audit.load(audit.PATHS["train_primary"])
        y = data[audit.TARGET].astype(int).to_numpy()
        groups = data[audit.GROUP].to_numpy()
        splitter = audit.StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=audit.SEED)
        for fit, held in splitter.split(data, y, groups):
            self.assertFalse(set(groups[fit]) & set(groups[held]))

    def test_threshold_is_chosen_from_training_probabilities(self):
        threshold = audit.choose_threshold(np.array([0, 0, 1, 1]), np.array([.1, .2, .8, .9]))
        self.assertGreater(threshold, .2)
        self.assertLessEqual(threshold, .8)

    def test_validation_entry_point_is_locked(self):
        source = Path(audit.__file__).read_text()
        self.assertIn("selection_locked_before_validation", source)
        self.assertNotIn("test_core_no_population.csv", source)

    def test_every_condition_has_one_oof_prediction_per_training_variant(self):
        path = audit.OUT / "03_controlled_experiments/out_of_fold_predictions.csv"
        rows = pd.read_csv(path, usecols=["condition_id", audit.KEY])
        counts = rows.groupby("condition_id").agg(rows=(audit.KEY, "size"), unique=(audit.KEY, "nunique"))
        self.assertTrue((counts["rows"] == 27727).all())
        self.assertTrue((counts["unique"] == 27727).all())

    def test_saved_cv_fold_assignments_are_gene_disjoint_and_deterministic(self):
        path = audit.OUT / "03_controlled_experiments/fixed_cv_fold_assignments.csv"
        rows = pd.read_csv(path)
        self.assertEqual(len(rows), 27727)
        self.assertEqual(rows.groupby(audit.GROUP)["cv_fold"].nunique().max(), 1)
        self.assertEqual(set(rows["cv_fold"]), {0, 1, 2, 3, 4})

    def test_locked_protocol_checksum_matches(self):
        path = audit.OUT / "00_locked_protocol/locked_protocol.yaml"
        expected = (audit.OUT / "00_locked_protocol/locked_protocol.sha256").read_text().split()[0]
        self.assertEqual(audit.sha256(path), expected)

    def test_metrics_use_probability_and_stated_threshold(self):
        y = np.array([0, 0, 1, 1])
        probability = np.array([0.1, 0.6, 0.7, 0.9])
        result = audit.metric_values(y, probability, 0.65)
        self.assertEqual(result["threshold"], 0.65)
        self.assertEqual(result["mcc"], 1.0)
        self.assertAlmostEqual(result["brier"], np.mean((probability - y) ** 2))
        self.assertAlmostEqual(result["pr_auc"], 1.0)

    def test_threshold_selection_does_not_accept_heldout_labels(self):
        import inspect
        parameters = list(inspect.signature(audit.choose_threshold).parameters)
        self.assertEqual(parameters, ["y", "probability"])

    def test_corrected_governance_separates_roles(self):
        path = audit.OUT / "01_provenance_circularity/provenance_circularity_register_corrected.csv"
        rows = pd.read_csv(path).set_index("display_name")
        self.assertEqual(rows.loc["CADD PHRED", "proposed_role"], "OPTIONAL SENSITIVITY INPUT")
        self.assertEqual(rows.loc["AlphaMissense", "proposed_role"], "OPTIONAL SENSITIVITY INPUT")
        self.assertEqual(rows.loc["SIFT", "proposed_role"], "EXTERNAL COMPARATOR")

    def test_paired_gene_bootstrap_is_complete_and_deterministic(self):
        path = audit.OUT / "03_controlled_experiments/paired_gene_bootstrap_results.csv"
        rows = pd.read_csv(path)
        self.assertEqual(len(rows), 12)
        self.assertEqual(set(rows["bootstrap_samples"]), {5000})
        self.assertEqual(set(rows["seed"]), {20261018})
        self.assertEqual(set(rows["resampling_unit"]), {"reconciled_gene_symbol"})

    def test_refined_selection_is_single_and_training_only(self):
        path = audit.OUT / "03_controlled_experiments/recommended_training_only_configuration_refined.json"
        selected = json.loads(path.read_text())
        self.assertEqual(selected["recommended_condition_id"], "conservation_phylop100way__phastcons100way")
        self.assertNotIn("conservation_gerp_rs", selected["recommended_features"])
        self.assertFalse(selected["validation_used"])
        self.assertFalse(selected["test_used"])

    def test_required_figures_exist_as_nonempty_png(self):
        figure_names = {
            "cohort_reconciliation_flow", "split_integrity_summary",
            "provenance_circularity_risk_matrix", "evidence_role_summary",
            "conservation_heatmap_annotated", "computational_predictor_heatmap_annotated",
            "combined_evidence_heatmap_annotated", "clustered_evidence_heatmap_fixed",
            "flagged_redundancy_network", "conservation_subset_performance",
            "predictor_comparator_performance", "predictor_metric_summary",
            "paired_gene_bootstrap_forest", "oof_calibration_curves",
        }
        figures = list(audit.OUT.rglob("figures/*"))
        found = {(path.stem, path.suffix) for path in figures if path.stat().st_size > 0}
        for name in figure_names:
            self.assertIn((name, ".png"), found)

    def test_markdown_image_paths_resolve(self):
        report_dir = audit.OUT / "04_report"
        for name in ("PRE_MODELING_EVIDENCE_AUDIT_REVISED.md", "SUPERVISOR_PROGRESS_SUMMARY.md", "FIGURE_INDEX.md"):
            path = report_dir / name
            links = re.findall(r"!?(?:\[[^]]*\])\(([^)]+)\)", path.read_text())
            self.assertTrue(links, name)
            for link in links:
                self.assertTrue((path.parent / link).resolve().is_file(), f"{name}: {link}")


if __name__ == "__main__":
    unittest.main()
