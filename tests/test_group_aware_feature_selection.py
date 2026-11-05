from __future__ import annotations

import inspect
import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

from src.pre_modeling_audit import group_aware_feature_selection as selection


class GroupAwareFeatureSelectionTests(unittest.TestCase):
    def test_training_only_guard_rejects_heldout_and_unknown_paths(self):
        with self.assertRaises(RuntimeError):
            selection.guarded_training_load(selection.ROOT / "data/processed/modeling_splits/validation_core_no_population.csv")
        with self.assertRaises(RuntimeError):
            selection.guarded_training_load(selection.ROOT / "data/processed/modeling_splits/test_core_no_population.csv")

    def test_all_expected_subsets_are_enumerated(self):
        self.assertEqual(len(selection.subsets(selection.CONSERVATION)), 8)
        self.assertEqual(len(selection.subsets(selection.AA)), 4)
        self.assertEqual(len(selection.subsets(selection.CONSTRAINT)), 4)
        self.assertEqual(len(selection.subsets(selection.POPULATION)), 4)
        self.assertEqual(len(selection.subsets(selection.PREDICTORS)), 32)

    def test_mrmr_synthetic_relevance_and_determinism(self):
        rng = np.random.default_rng(7)
        y = np.tile([0, 1], 100)
        frame = pd.DataFrame({"relevant": y + rng.normal(0, .02, len(y)), "noise": rng.normal(size=len(y))})
        first = selection.mrmr_ranking(frame, y, ["relevant", "noise"])
        second = selection.mrmr_ranking(frame, y, ["relevant", "noise"])
        self.assertEqual(first, second)
        self.assertEqual(first[0]["feature"], "relevant")

    def test_inner_gene_splits_are_isolated(self):
        groups = np.repeat([f"g{i}" for i in range(40)], 4)
        y = np.tile([0, 0, 1, 1], 40)
        splitter = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=selection.SEED + 100)
        for fit, held in splitter.split(np.arange(len(y)), y, groups):
            self.assertFalse(set(groups[fit]) & set(groups[held]))

    def test_preprocessing_statistics_are_fit_only_on_fit_rows(self):
        model = selection.pipeline()
        x_fit = pd.DataFrame({"x": [0.0, np.nan, 2.0, 4.0], "z": [1.0, 2.0, 3.0, 4.0]})
        y_fit = np.array([0, 1, 0, 1])
        model.fit(x_fit, y_fit)
        self.assertAlmostEqual(model.named_steps["imputer"].statistics_[0], 2.0)
        model.predict_proba(pd.DataFrame({"x": [1000.0], "z": [1000.0]}))
        self.assertAlmostEqual(model.named_steps["imputer"].statistics_[0], 2.0)

    def test_subset_selection_is_deterministic_and_prefers_parsimony_on_tie(self):
        scores = {frozenset(): .5, frozenset(["a"]): .6, frozenset(["a", "b"]): .6}
        self.assertEqual(selection.choose_subset(scores), frozenset(["a"]))
        self.assertEqual(selection.choose_subset(scores), selection.choose_subset(scores))

    def test_dominance_normalisation_is_not_final_selection_rule(self):
        decision_section = inspect.getsource(selection.synthesise_decisions)
        self.assertNotIn("normalised_share", decision_section)
        self.assertNotIn("dominance_contribution", decision_section)

    def test_runner_contains_no_heldout_data_filename(self):
        source = inspect.getsource(selection)
        self.assertNotIn("validation_core", source)
        self.assertNotIn("test_core", source)

    def test_completed_fixed_subset_outputs_have_all_conditions(self):
        performance = pd.read_csv(selection.STAGE / "nested_wrapper_subset_performance.csv")
        fixed = performance[performance.analysis_type.eq("fixed_subset_outer_oof")]
        expected = {
            "conservation": 8,
            "amino_acid_properties": 4,
            "gene_constraint": 4,
            "population_frequency": 4,
            "computational_predictors": 32,
        }
        self.assertEqual(fixed.groupby("block").size().to_dict(), expected)
        self.assertEqual(fixed.groupby("block").subset_key.nunique().to_dict(), expected)

    def test_completed_mrmr_outputs_cover_every_fold_and_candidate(self):
        rankings = pd.read_csv(selection.STAGE / "mrmr_fold_rankings.csv")
        for block in selection.block_definitions([]):
            rows = rankings[rankings.block.eq(block["block"])]
            self.assertEqual(set(rows.outer_fold), {0, 1, 2, 3, 4})
            for _, fold_rows in rows.groupby("outer_fold"):
                self.assertEqual(set(fold_rows.feature), set(block["candidates"]))
                self.assertEqual(set(fold_rows["rank"]), set(range(1, len(block["candidates"]) + 1)))

    def test_completed_run_metadata_confirms_training_only(self):
        metadata = json.loads((selection.STAGE / "run_metadata.json").read_text())
        serialized = json.dumps(metadata).lower()
        self.assertNotIn("validation_core", serialized)
        self.assertNotIn("test_core", serialized)
        self.assertFalse(metadata["validation_accessed"])
        self.assertFalse(metadata["test_accessed"])
        self.assertEqual(metadata["training_rows"], 27727)

    def test_completed_figures_exist_as_png(self):
        names = [
            "mrmr_selection_frequency", "subset_performance_confidence_intervals",
            "subset_selection_stability", "performance_vs_subset_complexity",
            "block_level_decision_summary",
        ]
        for name in names:
            path = selection.STAGE / "figures" / f"{name}.png"
            self.assertTrue(path.is_file(), path)
            self.assertGreater(path.stat().st_size, 1000, path)

    def test_analysis_specification_checksum_is_locked(self):
        self.assertEqual(selection.verify_specification(), selection.SPEC_HASH.read_text().split()[0])


if __name__ == "__main__":
    unittest.main()
