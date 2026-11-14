import unittest

from src.modeling.evidence_report import agreement_statement, interpret_contribution


class EvidenceReportTests(unittest.TestCase):
    def setUp(self):
        self.thresholds = {
            "conservation": {"weak_boundary": 0.2, "strong_boundary": 1.0},
            "gene_constraint": {"weak_boundary": 0.2, "strong_boundary": 1.0},
        }

    def test_signed_contribution_interpretation(self):
        self.assertEqual(interpret_contribution(1.2, 0.2, 1.0), ("strongly", "pathogenicity"))
        self.assertEqual(interpret_contribution(-0.5, 0.2, 1.0), ("moderately", "benignity"))

    def test_conflicting_groups_are_reported(self):
        text = agreement_statement(
            {"conservation": 1.2, "gene_constraint": -0.7}, self.thresholds
        )
        self.assertIn("conflicting evidence", text)
        self.assertIn("Conservation", text)
        self.assertIn("Gene constraint", text)

    def test_weak_groups_do_not_create_false_conflict(self):
        text = agreement_statement(
            {"conservation": 0.1, "gene_constraint": -0.1}, self.thresholds
        )
        self.assertIn("weak", text)


if __name__ == "__main__":
    unittest.main()
