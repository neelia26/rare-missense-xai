from __future__ import annotations

import math
import unittest
from pathlib import Path

import numpy as np

from src.pre_modeling_audit import dominance_hierarchical_partitioning as dominance
from src.pre_modeling_audit import run_pre_modeling_audit as audit


class DominanceHierarchicalPartitioningTests(unittest.TestCase):
    def test_exact_shapley_weighting_and_efficiency(self):
        features = ["a", "b", "c"]
        values = {
            frozenset(): 0.10,
            frozenset(["a"]): 0.20,
            frozenset(["b"]): 0.18,
            frozenset(["c"]): 0.14,
            frozenset(["a", "b"]): 0.27,
            frozenset(["a", "c"]): 0.23,
            frozenset(["b", "c"]): 0.22,
            frozenset(features): 0.30,
        }
        result = dominance.shapley_contributions(values, features)
        self.assertAlmostEqual(sum(result.values()), 0.20, places=12)
        expected_a = (0.20 - 0.10) / 3 + ((0.27 - 0.18) + (0.23 - 0.14)) / 6 + (0.30 - 0.22) / 3
        self.assertAlmostEqual(result["a"], expected_a, places=12)

    def test_normalised_shares_sum_to_one(self):
        shares = dominance.normalised_shares({"a": 0.12, "b": 0.08}, 0.20)
        self.assertAlmostEqual(sum(value for value in shares.values() if value is not None), 1.0, places=12)

    def test_negative_contributions_are_preserved(self):
        values = {frozenset(): 0.5, frozenset(["a"]): 0.4, frozenset(["b"]): 0.7, frozenset(["a", "b"]): 0.6}
        result = dominance.shapley_contributions(values, ["a", "b"])
        self.assertLess(result["a"], 0)
        self.assertGreater(result["b"], 0)

    def test_zero_and_negative_gain_are_not_normalised(self):
        self.assertEqual(dominance.normalised_shares({"a": 0.0}, 0.0), {"a": None})
        self.assertEqual(dominance.normalised_shares({"a": -0.1}, -0.1), {"a": None})

    def test_gene_cluster_bootstrap_is_deterministic_and_paired(self):
        old = dominance.N_BOOTSTRAP
        dominance.N_BOOTSTRAP = 50
        try:
            y = np.array([0, 0, 1, 1, 0, 1])
            groups = np.array(["g1", "g1", "g2", "g2", "g3", "g3"])
            probabilities = {
                frozenset(): np.array([.1, .2, .6, .8, .3, .7]),
                frozenset(["x"]): np.array([.1, .2, .8, .9, .2, .8]),
            }
            values1, phi1 = dominance.bootstrap_block(y, groups, ["x"], probabilities)
            values2, phi2 = dominance.bootstrap_block(y, groups, ["x"], probabilities)
            np.testing.assert_array_equal(values1[frozenset()], values2[frozenset()])
            np.testing.assert_array_equal(phi1["x"], phi2["x"])
            np.testing.assert_allclose(phi1["x"], values1[frozenset(["x"])] - values1[frozenset()])
        finally:
            dominance.N_BOOTSTRAP = old

    def test_stage_has_no_validation_or_test_data_path(self):
        source = Path(dominance.__file__).read_text(encoding="utf-8")
        self.assertNotIn("validation_core", source)
        self.assertNotIn("test_core", source)
        self.assertNotIn("validation_confirmation(", source)

    def test_frozen_and_split_checksums_match_locked_protocol(self):
        locked = __import__("json").loads((audit.OUT / "00_locked_protocol/locked_protocol.yaml").read_text())
        self.assertEqual(audit.sha256(audit.PATHS["frozen"]), locked["inputs"]["frozen"]["sha256"])
        self.assertEqual(audit.sha256(audit.PATHS["assignments"]), locked["inputs"]["assignments"]["sha256"])

    def test_locked_split_remains_gene_and_variant_disjoint(self):
        report = __import__("json").loads((audit.OUT / "00_locked_protocol/cohort_validation/split_integrity_report.json").read_text())
        self.assertEqual(report["status"], "PASS")
        self.assertTrue(report["split_union_equals_cohort"])
        self.assertFalse(any(report["gene_overlap"].values()))
        self.assertFalse(any(report["variant_overlap"].values()))

    def test_completed_outputs_cover_every_governed_subset(self):
        rows = __import__("pandas").read_csv(dominance.STAGE / "subset_performance.csv")
        self.assertEqual(rows.groupby("block").size().to_dict(), {
            "primary_conservation": 8,
            "primary_amino_acid": 4,
            "primary_gene_constraint": 4,
            "population_inclusive_sensitivity": 4,
            "computational_predictor_sensitivity": 32,
        })
        self.assertTrue((rows["row_count"] == 27727).all())

    def test_completed_shapley_efficiency_is_exact(self):
        rows = __import__("pandas").read_csv(dominance.STAGE / "dominance_contributions.csv")
        self.assertTrue((rows["efficiency_residual"].abs() <= dominance.TOLERANCE).all())
        for _, group in rows.groupby("block"):
            if group["normalisation_defined"].all():
                self.assertAlmostEqual(group["normalised_share"].sum(), 1.0, places=9)

    def test_run_metadata_confirms_no_heldout_access(self):
        metadata = __import__("json").loads((dominance.STAGE / "run_metadata.json").read_text())
        self.assertFalse(metadata["validation_accessed"])
        self.assertFalse(metadata["test_accessed"])
        self.assertFalse(metadata["authoritative_recommendation_modified"])
        self.assertFalse(metadata["validation_command_modified"])
        self.assertEqual(metadata["analysis_specification_sha256"], audit.sha256(dominance.SPEC))

    def test_required_stage_figures_are_nonempty(self):
        dominance_names = {
            "conservation_dominance_contributions",
            "amino_acid_dominance_contributions", "constraint_dominance_contributions",
            "population_sensitivity_dominance", "predictor_dominance_contributions",
        }
        subset_names = {
            "conservation_subset_performance", "subset_performance_vs_complexity", "amino_acid_subset_performance",
            "gene_constraint_subset_performance", "population_subset_performance",
            "computational_predictor_subset_performance",
        }
        for name in dominance_names:
            path = dominance.STAGE / "figures" / f"{name}.png"
            self.assertTrue(path.is_file() and path.stat().st_size > 0, str(path))
        for name in subset_names:
            path = dominance.STAGE / "figures" / "subset_combinations" / f"{name}.png"
            self.assertTrue(path.is_file() and path.stat().st_size > 0, str(path))

    def test_conservation_mcc_reconciliation_is_exact(self):
        rows = __import__("pandas").read_csv(dominance.STAGE / "conservation_mcc_reconciliation.csv")
        self.assertEqual(len(rows), 8)
        self.assertTrue(np.allclose(rows["previous_reported_mean_fold_mcc"], rows["recomputed_unweighted_mean_fold_mcc"], atol=5e-9))
        for column in [name for name in rows if name.endswith("_identical")]:
            self.assertTrue(rows[column].all(), column)
        self.assertTrue(rows["exact_cause"].str.contains("aggregation method only").all())

    def test_requested_predictor_reduction_comparison_exists(self):
        rows = __import__("pandas").read_csv(dominance.STAGE / "paired_subset_comparisons.csv")
        selected = rows[rows["comparison_id"] == "all_five_vs_revel_cadd_alphamissense"]
        self.assertEqual(len(selected), 1)
        row = selected.iloc[0]
        self.assertEqual(row["bootstrap_samples"], 5000)
        self.assertLess(row["ci_2_5"], 0)
        self.assertGreater(row["ci_97_5"], 0)

    def test_addendum_report_has_required_title(self):
        report = (dominance.STAGE / "DOMINANCE_HIERARCHICAL_PARTITIONING_REPORT.md").read_text()
        self.assertTrue(report.startswith("# Rare Missense-SNV Pre-Modelling Evidence Audit:\n## Dominance and Hierarchical-Partitioning Addendum"))


if __name__ == "__main__":
    unittest.main()
