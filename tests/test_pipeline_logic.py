#!/usr/bin/env python3

from __future__ import annotations

import importlib.util
import json
import tempfile
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


def load_root_script(filename: str):
    path = ROOT / filename
    spec = importlib.util.spec_from_file_location(filename.replace(".py", ""), path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


class PipelineLogicTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.extract = load_script("02_extract_candidates.py")
        cls.vep_runner = load_script("05_run_vep.py")
        cls.vep = load_script("06_parse_vep.py")
        cls.aa = load_script("11_add_amino_acid_features.py")
        cls.build = load_script("12_build_final_dataset.py")
        cls.split = load_script("13_split_by_gene.py")
        cls.runner = load_root_script("run_pipeline.py")

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

    def test_vep_checkpoint_fingerprint_depends_on_variants(self):
        common = (
            "https://rest.ensembl.org/vep/homo_sapiens/region",
            {"mane": 1},
            "GRCh38",
            "116",
        )
        first = self.vep_runner.request_fingerprint(*common, ["1 10 v1 A G . . ."])
        second = self.vep_runner.request_fingerprint(*common, ["1 11 v2 C T . . ."])
        self.assertNotEqual(first, second)
        path = self.vep_runner.checkpoint_path(Path("cache"), 0, first)
        self.assertIn(first[:16], path.name)

    def test_vep_release_mismatch_stops_the_run(self):
        self.assertEqual(
            self.vep_runner.validate_api_release({"release": 116}, "116"), "116"
        )
        with self.assertRaises(RuntimeError):
            self.vep_runner.validate_api_release({"release": 117}, "116")

    def test_old_number_only_vep_checkpoint_is_not_reused(self):
        variants = ["1 10 v1 A G . . ."]
        fingerprint = self.vep_runner.request_fingerprint(
            "https://rest.ensembl.org/vep/homo_sapiens/region",
            {"mane": 1},
            "GRCh38",
            "116",
            variants,
        )
        with tempfile.TemporaryDirectory() as temporary:
            old_checkpoint = Path(temporary) / "batch_00000.json"
            old_checkpoint.write_text(json.dumps([{"input": variants[0]}]), encoding="utf-8")
            loaded = self.vep_runner.load_checkpoint(
                old_checkpoint, fingerprint, variants, "116"
            )
        self.assertIsNone(loaded)

    def test_full_plan_keeps_direct_gnomad_counts_optional(self):
        config = {
            "rarity": {"maximum_global_af": 0.01},
            "gnomad": {
                "source_release": "gnomAD v4.1",
                "direct_counts_enabled": False,
            },
        }
        plan = self.runner.full_stage_arguments(
            Path("config/pipeline_config.json"), config, Path("constraint.tsv")
        )
        stages = [stage for stage, _ in plan]
        self.assertNotIn("08b", stages)
        stage09 = dict(plan)["09"]
        self.assertTrue(any(value.endswith("08_interpro.csv") for value in stage09))


if __name__ == "__main__":
    unittest.main()
