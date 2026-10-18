#!/usr/bin/env python3
"""Extract conservation and separate optional predictor values from a tabix dbNSFP file."""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    load_config,
    normalize_chromosome,
    parse_float,
    read_csv,
    write_csv,
    write_json,
)


CONSERVATION_OUTPUTS = {
    "GERP++_RS": "conservation_gerp_rs",
    "phyloP100way_vertebrate": "conservation_phylop100way",
    "phastCons100way_vertebrate": "conservation_phastcons100way",
}


def find_header(tabix_file: Any) -> list[str]:
    headers = list(tabix_file.header)
    if not headers:
        raise ValueError("dbNSFP file has no tabix header")
    raw = headers[-1].lstrip("#")
    return raw.split("\t")


def first_present(columns: list[str], candidates: list[str]) -> str:
    for candidate in candidates:
        if candidate in columns:
            return candidate
    raise ValueError(f"dbNSFP header lacks all expected columns: {candidates}")


def numeric_max(values: list[str]) -> str:
    parsed: list[float] = []
    for value in values:
        for item in value.replace(",", ";").split(";"):
            number = parse_float(item)
            if number is not None:
                parsed.append(number)
    return "" if not parsed else str(max(parsed))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/09_constraint.csv"))
    parser.add_argument("--dbnsfp-file", required=True, help="bgzip-compressed, tabix-indexed dbNSFP variant file")
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/10_dbnsfp.csv"))
    parser.add_argument(
        "--optional-predictor-output",
        default=str(PROJECT_ROOT / "data/interim/10_optional_predictor_scores.csv"),
    )
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    args = parser.parse_args()
    try:
        import pysam  # type: ignore
    except ImportError as error:
        raise RuntimeError("Install requirements.txt; pysam is required for dbNSFP tabix lookup") from error

    config = load_config(args.config)
    input_path = Path(args.input)
    dbnsfp_path = Path(args.dbnsfp_file).expanduser().resolve()
    output_path = Path(args.output)
    optional_path = Path(args.optional_predictor_output)
    rows = read_csv(input_path)
    table = pysam.TabixFile(str(dbnsfp_path))
    columns = find_header(table)
    column_index = {name: index for index, name in enumerate(columns)}
    chromosome_column = first_present(columns, ["chr", "#chr", "chrom", "chromosome"])
    position_column = first_present(columns, ["pos(1-based)", "pos", "position"])
    reference_column = first_present(columns, ["ref", "reference"])
    alternate_column = first_present(columns, ["alt", "alternate"])

    conservation_fields = config["dbnsfp"]["conservation_fields"]
    predictor_fields = config["dbnsfp"]["optional_predictor_fields"]
    missing_requested = [field for field in conservation_fields if field not in column_index]
    if missing_requested:
        raise ValueError(f"Requested conservation columns missing from dbNSFP: {missing_requested}")

    core_rows: list[dict[str, Any]] = []
    optional_rows: list[dict[str, Any]] = []
    for row in rows:
        chromosome = normalize_chromosome(row["chromosome"])
        position = int(row["position"])
        contig = chromosome if chromosome in table.contigs else f"chr{chromosome}"
        matches: list[dict[str, str]] = []
        if contig in table.contigs:
            try:
                lines = table.fetch(contig, position - 1, position)
            except ValueError:
                lines = []
            for line in lines:
                values = line.split("\t")
                if len(values) != len(columns):
                    continue
                item = dict(zip(columns, values))
                if (
                    normalize_chromosome(item[chromosome_column]) == chromosome
                    and int(item[position_column]) == position
                    and item[reference_column].upper() == row["reference_allele"].upper()
                    and item[alternate_column].upper() == row["alternate_allele"].upper()
                ):
                    matches.append(item)
        combined = dict(row)
        combined["dbnsfp_release"] = args.source_release
        combined["dbnsfp_lookup_status"] = "FOUND" if matches else "NOT_FOUND"
        for source, destination in CONSERVATION_OUTPUTS.items():
            combined[destination] = numeric_max([item.get(source, "") for item in matches])
        core_rows.append(combined)

        optional: dict[str, Any] = {
            "variant_key": row["variant_key"],
            "dbnsfp_release": args.source_release,
            "dbnsfp_lookup_status": combined["dbnsfp_lookup_status"],
        }
        for field in predictor_fields:
            destination = "predictor_" + field.lower().replace("++", "pp").replace("-", "_")
            optional[destination] = (
                numeric_max([item.get(field, "") for item in matches]) if field in column_index else ""
            )
        optional_rows.append(optional)

    core_added = ["dbnsfp_release", "dbnsfp_lookup_status", *CONSERVATION_OUTPUTS.values()]
    write_csv(output_path, core_rows, list(rows[0]) + [x for x in core_added if x not in rows[0]])
    write_csv(optional_path, optional_rows, list(optional_rows[0]))
    audit = {
        **audit_basic_rows(core_rows),
        "source_release": args.source_release,
        "lookup_status_counts": dict(Counter(row["dbnsfp_lookup_status"] for row in core_rows)),
        "conservation_reduction": "Maximum valid numeric value when dbNSFP contains multiple transcript-aligned values.",
        "predictor_policy": (
            "Composite predictor scores are written to a separate file and are not included in the core evidence table."
        ),
        "predictor_columns_available": [field for field in predictor_fields if field in column_index],
    }
    audit_path = PROJECT_ROOT / "audits/10_dbnsfp_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "10_enrich_dbnsfp",
        [input_path, dbnsfp_path],
        [output_path, optional_path, audit_path],
        {"source_release": args.source_release, "conservation_fields": conservation_fields},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

