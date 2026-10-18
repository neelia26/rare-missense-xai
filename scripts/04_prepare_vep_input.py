#!/usr/bin/env python3
"""Convert the frozen candidate pool into a VCF accepted by Ensembl VEP."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import PROJECT_ROOT, append_stage_manifest, read_csv, require_columns


def chromosome_sort_key(value: str) -> tuple[int, str]:
    chromosome = value.removeprefix("chr")
    if chromosome.isdigit():
        return int(chromosome), ""
    return {"X": 23, "Y": 24, "MT": 25, "M": 25}.get(chromosome, 99), chromosome


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/03_candidate_pool.csv"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/raw/vep/04_candidates.vcf"))
    args = parser.parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)
    rows = read_csv(input_path)
    require_columns(
        rows,
        ["variant_key", "chromosome", "position", "reference_allele", "alternate_allele"],
        "04_prepare_vep_input",
    )
    rows.sort(key=lambda row: (chromosome_sort_key(row["chromosome"]), int(row["position"])))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=GRCh38\n")
        handle.write('##INFO=<ID=SOURCE_KEY,Number=1,Type=String,Description="Pipeline variant key">\n')
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for row in rows:
            handle.write(
                "\t".join(
                    [
                        row["chromosome"].removeprefix("chr"),
                        row["position"],
                        row["variant_key"],
                        row["reference_allele"],
                        row["alternate_allele"],
                        ".",
                        "PASS",
                        f"SOURCE_KEY={row['variant_key']}",
                    ]
                )
                + "\n"
            )
    manifest = append_stage_manifest(
        "04_prepare_vep_input",
        [input_path],
        [output_path],
        {"rows": len(rows), "format": "VCFv4.2", "assembly": "GRCh38"},
    )
    print(json.dumps({"rows": len(rows), "output": str(output_path), "manifest": str(manifest)}, indent=2))


if __name__ == "__main__":
    main()

