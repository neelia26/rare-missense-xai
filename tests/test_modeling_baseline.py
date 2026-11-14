import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from src.modeling.common import FrozenInputError, load_frozen_training_data
from src.modeling.logistic_baseline import (
    calibration_intercept_slope,
    classification_metrics,
    select_mcc_threshold,
)


ROOT = Path(__file__).resolve().parents[1]


class FrozenPrimaryInputTests(unittest.TestCase):
    def test_frozen_primary_input_resolves_to_47_features(self):
        data = load_frozen_training_data(ROOT)
        self.assertEqual(len(data.features), 47)
        self.assertEqual(
            {name: len(values) for name, values in data.evidence_groups.items()},
            {
                "conservation": 2,
                "gene_constraint": 2,
                "amino_acid_properties": 8,
                "protein_context": 35,
            },
        )

    def test_primary_input_excludes_governed_sensitivity_features(self):
        data = load_frozen_training_data(ROOT)
        forbidden = {
            "conservation_gerp_rs",
            "aa_polarity_delta",
            "gnomad_global_af",
            "gnomad_max_population_af",
        }
        self.assertTrue(forbidden.isdisjoint(data.features))
        self.assertFalse(any(name.startswith("predictor_") for name in data.features))

    def test_hash_change_fails_before_training_rows_are_used(self):
        source = ROOT / "config/grouped_additive_primary_input.json"
        config = json.loads(source.read_text())
        config["dataset"]["training_sha256"] = "0" * 64
        with tempfile.TemporaryDirectory() as directory:
            changed = Path(directory) / "changed.json"
            changed.write_text(json.dumps(config))
            with self.assertRaises(FrozenInputError):
                load_frozen_training_data(ROOT, input_config_path=changed)


class BaselineMetricTests(unittest.TestCase):
    def test_metrics_use_all_four_confusion_matrix_cells(self):
        labels = np.array([0, 0, 1, 1])
        probabilities = np.array([0.1, 0.7, 0.8, 0.2])
        metrics = classification_metrics(labels, probabilities, 0.5)
        self.assertEqual(
            (metrics["tn"], metrics["fp"], metrics["fn"], metrics["tp"]),
            (1, 1, 1, 1),
        )
        self.assertAlmostEqual(metrics["mcc"], 0.0)

    def test_threshold_selection_is_deterministic(self):
        labels = np.array([0, 0, 0, 1, 1, 1])
        probabilities = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 0.95])
        first = select_mcc_threshold(labels, probabilities)
        second = select_mcc_threshold(labels, probabilities)
        self.assertEqual(first, second)
        self.assertGreater(first[1], 0.0)

    def test_calibration_fit_returns_finite_values(self):
        labels = np.array([0, 0, 0, 1, 1, 1])
        probabilities = np.array([0.05, 0.2, 0.4, 0.6, 0.8, 0.95])
        intercept, slope = calibration_intercept_slope(labels, probabilities)
        self.assertTrue(np.isfinite(intercept))
        self.assertTrue(np.isfinite(slope))
        self.assertGreater(slope, 0)


if __name__ == "__main__":
    unittest.main()
