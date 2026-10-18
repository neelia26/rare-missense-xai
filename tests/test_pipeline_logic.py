#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_script(filename: str):
    path = ROOT / "scripts" / filename
    spec = importlib.util.spec_from_file_location(filename.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class PipelineLogicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.extract = load_script("02_extract_candidates.py")
        cls.vep = load_script("06_parse_vep.py")
        cls.aa = load_script("11_add_amino_acid_features.py")
        cls.build = load_script("12_build_final_dataset.py")
        cls.split = load_script("13_split_by_gene.py")

    def test_clinvar_row_extraction(self):
        fields = [
            "1", "100", "123", "A", "G", ".", ".",
            "CLNVC=single_nucleotide_variant;MC=SO:0001583|missense_variant;"
            "CLNSIG=Pathogenic;CLNREVSTAT=criteria_provided,_single_submitter;"
            "GENEINFO=GENE1:1;ALLELEID=10;RS=99",
        ]
        row, status = self.extract.row_from_vcf(fields, {"Pathogenic": 1}, 1, "20261001")
        self.assertEqual(status, "eligible")
        self.assertEqual(row["variant_key"], "GRCh38:1:100:A:G")
        self.assertEqual(row["binary_label"], 1)

    def test_transcript_prefers_mane_for_matching_gene(self):
        row = {"gene_symbol": "GENE1", "all_gene_info": "GENE1:1"}
        result = {
            "transcript_consequences": [
                {
                    "gene_symbol": "GENE1", "transcript_id": "ENST2",
                    "consequence_terms": ["missense_variant"], "amino_acids": "A/G",
                    "protein_start": 5, "protein_end": 5, "canonical": 1,
                },
                {
                    "gene_symbol": "GENE1", "transcript_id": "ENST1",
                    "consequence_terms": ["missense_variant"], "amino_acids": "A/G",
                    "protein_start": 5, "protein_end": 5, "mane_select": "NM_1.1",
                },
            ]
        }
        selected, basis = self.vep.choose_transcript(row, result)
        self.assertEqual(selected["transcript_id"], "ENST1")
        self.assertIn("MANE_SELECT", basis)

    def test_frequency_not_reported_is_not_zero(self):
        row = {"reference_allele": "A", "alternate_allele": "G"}
        values = self.vep.population_values(row, {"colocated_variants": []}, 0.01)
        self.assertEqual(values["population_lookup_status"], "NOT_REPORTED")
        self.assertEqual(values["gnomad_global_af"], "")

    def test_amino_acid_features(self):
        values = self.aa.values(
            {
                "reconciled_reference_amino_acid": "A",
                "reconciled_alternate_amino_acid": "G",
                "reconciled_protein_position": "5",
                "uniprot_protein_length": "100",
            }
        )
        self.assertEqual(values["aa_feature_status"], "COMPLETE")
        self.assertGreater(values["aa_grantham_score"], 0)
        self.assertAlmostEqual(values["aa_relative_protein_position"], 0.05)

    def test_final_exclusion_rules(self):
        row = {
            "binary_label": "1",
            "strict_missense_eligible": "1",
            "uniprot_validation_status": "MATCH",
            "reconciled_gene_symbol": "GENE1",
            "gnomad_rarity_status": "NOT_RARE",
        }
        self.assertEqual(self.build.exclusion_reason(row, True), "NOT_RARE")

    def test_gene_assignment_is_disjoint(self):
        rows = []
        for gene_number in range(20):
            for label in ("0", "1"):
                rows.append(
                    {
                        "variant_key": f"v{gene_number}_{label}",
                        "reconciled_gene_symbol": f"G{gene_number}",
                        "binary_label": label,
                    }
                )
        groups = self.split.group_statistics(rows)
        assignment = self.split.assign_groups(
            groups, {"train": 0.7, "validation": 0.15, "test": 0.15}, 42
        )
        self.assertEqual(set(assignment), set(groups))
        self.assertEqual(set(assignment.values()), {"train", "validation", "test"})


if __name__ == "__main__":
    unittest.main()
