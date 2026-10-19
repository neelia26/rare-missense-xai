#!/usr/bin/env python3

from __future__ import annotations

import csv
import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_audit():
    path = ROOT / "scripts" / "15_audit_circularity.py"
    spec = importlib.util.spec_from_file_location("audit_circularity", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class LabelGuard(dict):
    def get(self, key, default=None):
        if key == "binary_label":
            raise AssertionError("test label was read")
        return super().get(key, default)


class CircularityAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = load_audit()

    def test_governed_missing_markers(self):
        for value in ("", ".", "NA", "NaN", "None", "null", None):
            self.assertTrue(self.audit.is_missing(value))
        self.assertFalse(self.audit.is_missing("0"))

    def test_structural_test_check_never_reads_labels(self):
        rows = [
            LabelGuard(
                variant_key="v1",
                reconciled_gene_symbol="G1",
                binary_label="DO_NOT_READ",
            )
        ]
        result = self.audit.structural_rows(
            rows, ("variant_key", "reconciled_gene_symbol", "binary_label")
        )
        self.assertEqual(result["rows"], 1)

    def test_label_dependent_outputs_are_marked_training_only(self):
        rows = [
            {"binary_label": "0", "x": ""},
            {"binary_label": "1", "x": "1"},
        ]
        result = self.audit.missingness_rows(rows, ["x"])
        self.assertEqual(result[0]["analysis_partition"], "train")

    def test_metadata_is_never_a_predictor(self):
        self.assertEqual(
            self.audit.METADATA,
            {"variant_key", "reconciled_gene_symbol", "binary_label"},
        )

    def test_prohibited_uniprot_features_are_enumerated(self):
        self.assertIn("uniprot_feature_natural_variant", self.audit.PROHIBITED_UNIPROT)
        self.assertEqual(len(self.audit.PROHIBITED_UNIPROT), 4)

    def test_optional_predictors_are_not_primary(self):
        with (ROOT / "data/processed/evidence_lineage_registry.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))
        optional = [r for r in rows if r["evidence_group"] == "optional_computational_predictor"]
        self.assertEqual(len(optional), 5)
        self.assertTrue(all(r["approved_primary"] == "0" for r in optional))

    def test_no_population_governance_is_mandatory(self):
        self.assertIn("no-population versus population-inclusive", self.audit.MANDATORY_ABLATIONS)

    def test_pairwise_complete_numeric_correlation(self):
        rows = [
            {"a": "1", "b": "2"},
            {"a": "2", "b": "4"},
            {"a": "", "b": "6"},
            {"a": "3", "b": "6"},
        ]
        result = self.audit.numeric_correlations(rows, ["a", "b"])[0]
        self.assertEqual(result["paired_observations"], 3)
        self.assertEqual(result["pearson"], 1.0)

    def test_constant_and_all_missing_fields_are_safe(self):
        rows = [{"a": "1", "b": ""}, {"a": "1", "b": ""}, {"a": "1", "b": ""}]
        numeric = self.audit.numeric_correlations(rows, ["a", "b"])[0]
        self.assertEqual(numeric["pearson"], "")
        self.assertNotEqual(numeric["pearson_status"], "VALID")
        value, paired, status = self.audit.corrected_cramers_v(
            ["x", "x", "x"], ["", "", ""]
        )
        self.assertIsNone(value)
        self.assertEqual(paired, 0)
        self.assertEqual(status, "INSUFFICIENT_PAIRED_OBSERVATIONS")

    def test_lineage_and_primary_registry_are_consistent(self):
        with (ROOT / "data/processed/modeling_feature_registry.csv").open(newline="") as f:
            primary = {r["field"] for r in csv.DictReader(f) if r["primary_model_feature"] == "1"}
        with (ROOT / "data/processed/evidence_lineage_registry.csv").open(newline="") as f:
            lineage = {r["field"] for r in csv.DictReader(f) if r["approved_primary"] == "1"}
        self.assertEqual(primary, lineage)

    def test_analysis_functions_are_deterministic(self):
        rows = [
            {"binary_label": "0", "a": "1", "b": "2", "c": "x", "d": "y"},
            {"binary_label": "1", "a": "2", "b": "4", "c": "y", "d": "x"},
            {"binary_label": "0", "a": "3", "b": "6", "c": "x", "d": "y"},
        ]
        first = (
            self.audit.missingness_rows(rows, ["a", "b"]),
            self.audit.numeric_correlations(rows, ["a", "b"]),
            self.audit.categorical_associations(rows, ["c", "d"]),
        )
        second = (
            self.audit.missingness_rows(rows, ["a", "b"]),
            self.audit.numeric_correlations(rows, ["a", "b"]),
            self.audit.categorical_associations(rows, ["c", "d"]),
        )
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
