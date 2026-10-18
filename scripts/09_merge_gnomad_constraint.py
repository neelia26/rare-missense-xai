#!/usr/bin/env python3
"""Merge a versioned gnomAD gene-constraint table without guessing ambiguous genes."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    open_text_auto,
    read_csv,
    write_csv,
    write_json,
)


ADDED_COLUMNS = [
    "gnomad_constraint_release",
    "gnomad_constraint_match_status",
    "gnomad_constraint_match_basis",
    "gnomad_constraint_gene",
    "gnomad_constraint_transcript",
    "gnomad_missense_z",
    "gnomad_loeuf",
]


def detect_delimiter(path: Path) -> str:
    with open_text_auto(path) as handle:
        sample = "".join(handle.readline() for _ in range(5))
    return "\t" if sample.count("\t") >= sample.count(",") else ","


def first_column(fieldnames: list[str], requested: str, candidates: list[str]) -> str:
    if requested:
        if requested not in fieldnames:
            raise ValueError(f"Requested column {requested!r} was not found")
        return requested
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"None of these columns were found: {', '.join(candidates)}")


def clean_transcript(value: str) -> str:
    return value.strip().split(".", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/08_interpro.csv"))
    parser.add_argument("--constraint-file", required=True)
    parser.add_argument("--source-release", required=True, help="For example: gnomAD v4.1")
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/09_constraint.csv"))
    parser.add_argument("--gene-column", default="")
    parser.add_argument("--transcript-column", default="")
    parser.add_argument("--missense-z-column", default="")
    parser.add_argument("--loeuf-column", default="")
    args = parser.parse_args()
    input_path = Path(args.input)
    constraint_path = Path(args.constraint_file).expanduser().resolve()
    output_path = Path(args.output)
    rows = read_csv(input_path)

    delimiter = detect_delimiter(constraint_path)
    with open_text_auto(constraint_path) as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        gene_column = first_column(fieldnames, args.gene_column, ["gene", "gene_symbol", "symbol"])
        transcript_column = first_column(
            fieldnames, args.transcript_column, ["transcript", "transcript_id", "canonical_transcript"]
        )
        missense_column = first_column(
            fieldnames, args.missense_z_column, ["mis.z_score", "mis_z", "mis_z_score", "missense_z"]
        )
        loeuf_column = first_column(
            fieldnames, args.loeuf_column, ["lof.oe_ci.upper", "oe_lof_upper", "LOEUF", "loeuf"]
        )
        constraint_rows = list(reader)

    by_transcript: dict[str, dict[str, str]] = {}
    by_gene: dict[str, list[dict[str, str]]] = defaultdict(list)
    for item in constraint_rows:
        transcript = clean_transcript(item.get(transcript_column, ""))
        gene = item.get(gene_column, "").strip()
        if transcript:
            by_transcript[transcript] = item
        if gene:
            by_gene[gene].append(item)

    enriched = []
    for row in rows:
        transcript = clean_transcript(row.get("reconciled_transcript_id", ""))
        gene = row.get("reconciled_gene_symbol", "") or row.get("gene_symbol", "")
        item = by_transcript.get(transcript)
        basis = "TRANSCRIPT"
        status = "MATCHED"
        if item is None:
            gene_candidates = by_gene.get(gene, [])
            if len(gene_candidates) == 1:
                item = gene_candidates[0]
                basis = "UNIQUE_GENE_FALLBACK"
            elif len(gene_candidates) > 1:
                basis = ""
                status = "AMBIGUOUS_GENE"
            else:
                basis = ""
                status = "NOT_FOUND"
        values = {column: "" for column in ADDED_COLUMNS}
        values.update(
            {
                "gnomad_constraint_release": args.source_release,
                "gnomad_constraint_match_status": status,
                "gnomad_constraint_match_basis": basis,
            }
        )
        if item is not None:
            values.update(
                {
                    "gnomad_constraint_gene": item.get(gene_column, ""),
                    "gnomad_constraint_transcript": item.get(transcript_column, ""),
                    "gnomad_missense_z": item.get(missense_column, ""),
                    "gnomad_loeuf": item.get(loeuf_column, ""),
                }
            )
        combined = dict(row)
        combined.update(values)
        enriched.append(combined)

    fieldnames_out = list(rows[0]) + [column for column in ADDED_COLUMNS if column not in rows[0]]
    write_csv(output_path, enriched, fieldnames_out)
    audit = {
        **audit_basic_rows(enriched),
        "source_release": args.source_release,
        "source_columns": {
            "gene": gene_column,
            "transcript": transcript_column,
            "missense_z": missense_column,
            "loeuf": loeuf_column,
        },
        "match_status_counts": dict(Counter(row["gnomad_constraint_match_status"] for row in enriched)),
        "rule": "Match the reconciled Ensembl transcript first; use gene fallback only when the source has one row for that gene.",
    }
    audit_path = PROJECT_ROOT / "audits/09_constraint_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "09_merge_gnomad_constraint",
        [input_path, constraint_path],
        [output_path, audit_path],
        {"source_release": args.source_release, **audit["source_columns"]},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

