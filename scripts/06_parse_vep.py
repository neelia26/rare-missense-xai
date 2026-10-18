#!/usr/bin/env python3
"""Select one auditable coding transcript and extract population evidence from VEP JSON."""

from __future__ import annotations

import argparse
import json
import re
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
    parse_float,
    read_csv,
    split_multi_value,
    variant_key,
    write_csv,
    write_json,
)


ADDED_COLUMNS = [
    "vep_result_status",
    "transcript_reconciliation_status",
    "transcript_selection_basis",
    "reconciled_gene_symbol",
    "reconciled_gene_id",
    "reconciled_transcript_id",
    "reconciled_mane_select",
    "reconciled_canonical_flag",
    "reconciled_consequence_terms",
    "reconciled_hgvsc",
    "reconciled_hgvsp",
    "reconciled_cds_position",
    "reconciled_codon_position",
    "reconciled_codon_change",
    "reconciled_amino_acid_change",
    "reconciled_reference_amino_acid",
    "reconciled_alternate_amino_acid",
    "reconciled_protein_position",
    "ensembl_protein_id",
    "uniprot_swissprot_accessions",
    "uniprot_trembl_accessions",
    "vep_domains",
    "strict_missense_eligible",
    "gnomad_exome_af",
    "gnomad_genome_af",
    "gnomad_global_af",
    "gnomad_max_population_af",
    "gnomad_frequency_fields",
    "gnomad_rarity_threshold",
    "gnomad_rarity_status",
    "population_lookup_status",
    "population_lookup_reason",
]


def input_key(result: dict[str, Any]) -> str | None:
    raw = str(result.get("input", ""))
    fields = raw.split()
    if len(fields) >= 5:
        chromosome, position, possible_id, reference, alternate = fields[:5]
        if possible_id.startswith("GRCh38:"):
            return possible_id
        try:
            return variant_key(chromosome, position, reference, alternate)
        except (TypeError, ValueError):
            pass
    identifier = str(result.get("id", ""))
    if identifier.startswith("GRCh38:"):
        return identifier
    return None


def genes(row: dict[str, str]) -> list[str]:
    values: list[str] = []
    for item in row.get("all_gene_info", "").split("|"):
        symbol = item.split(":", 1)[0].strip()
        if symbol and symbol not in values:
            values.append(symbol)
    if row.get("gene_symbol") and row["gene_symbol"] not in values:
        values.insert(0, row["gene_symbol"])
    return values


def is_missense(transcript: dict[str, Any]) -> bool:
    terms = transcript.get("consequence_terms") or []
    amino_acids = str(transcript.get("amino_acids", ""))
    parts = amino_acids.split("/", 1)
    return (
        "missense_variant" in terms
        and len(parts) == 2
        and len(parts[0]) == 1
        and len(parts[1]) == 1
        and bool(transcript.get("protein_start"))
    )


def transcript_rank(transcript: dict[str, Any], gene_order: list[str]) -> tuple:
    symbol = str(transcript.get("gene_symbol", ""))
    gene_rank = gene_order.index(symbol) if symbol in gene_order else 999
    mane_rank = 0 if transcript.get("mane_select") else 1
    canonical_rank = 0 if transcript.get("canonical") in {1, "1"} else 1
    coding_rank = 0 if transcript.get("biotype") == "protein_coding" else 1
    protein_length = transcript.get("protein_end") or transcript.get("protein_start") or 0
    return (
        gene_rank,
        mane_rank,
        canonical_rank,
        coding_rank,
        -int(protein_length),
        str(transcript.get("transcript_id", "")),
    )


def choose_transcript(row: dict[str, str], result: dict[str, Any]) -> tuple[dict | None, str]:
    gene_order = genes(row)
    candidates = [
        transcript
        for transcript in result.get("transcript_consequences", [])
        if is_missense(transcript) and transcript.get("gene_symbol") in gene_order
    ]
    if not candidates:
        return None, "NO_GENE_MATCHED_MISSENSE_TRANSCRIPT"
    selected = min(candidates, key=lambda item: transcript_rank(item, gene_order))
    basis = []
    basis.append("PRIMARY_CLINVAR_GENE" if selected.get("gene_symbol") == gene_order[0] else "ALTERNATE_CLINVAR_GENE")
    if selected.get("mane_select"):
        basis.append("MANE_SELECT")
    elif selected.get("canonical") in {1, "1"}:
        basis.append("CANONICAL")
    else:
        basis.append("DETERMINISTIC_RANK")
    return selected, ":".join(basis)


def rendered_domains(value: Any) -> str:
    rendered: list[str] = []
    for item in value if isinstance(value, list) else []:
        if isinstance(item, dict):
            rendered.append(f"{item.get('db', '')}:{item.get('name', '')}".strip(":"))
    return ";".join(item for item in rendered if item)


def accessions(value: Any) -> str:
    if isinstance(value, list):
        return ";".join(str(item) for item in value)
    return str(value or "")


def transcript_values(row: dict[str, str], result: dict[str, Any] | None) -> dict[str, Any]:
    values = {column: "" for column in ADDED_COLUMNS}
    values["strict_missense_eligible"] = 0
    if result is None:
        values["vep_result_status"] = "MISSING"
        values["transcript_reconciliation_status"] = "VEP_RESULT_MISSING"
        return values
    values["vep_result_status"] = "FOUND"
    transcript, basis = choose_transcript(row, result)
    if transcript is None:
        values["transcript_reconciliation_status"] = basis
        return values
    amino = str(transcript.get("amino_acids", ""))
    reference_aa, alternate_aa = amino.split("/", 1)
    cds_position = transcript.get("cds_start", "")
    codon_position = ((int(cds_position) - 1) % 3) + 1 if str(cds_position).isdigit() else ""
    values.update(
        {
            "transcript_reconciliation_status": "RECONCILED_MISSENSE",
            "transcript_selection_basis": basis,
            "reconciled_gene_symbol": transcript.get("gene_symbol", ""),
            "reconciled_gene_id": transcript.get("gene_id", ""),
            "reconciled_transcript_id": transcript.get("transcript_id", ""),
            "reconciled_mane_select": transcript.get("mane_select", ""),
            "reconciled_canonical_flag": transcript.get("canonical", ""),
            "reconciled_consequence_terms": ";".join(transcript.get("consequence_terms", [])),
            "reconciled_hgvsc": transcript.get("hgvsc", ""),
            "reconciled_hgvsp": transcript.get("hgvsp", ""),
            "reconciled_cds_position": cds_position,
            "reconciled_codon_position": codon_position,
            "reconciled_codon_change": transcript.get("codons", ""),
            "reconciled_amino_acid_change": amino,
            "reconciled_reference_amino_acid": reference_aa,
            "reconciled_alternate_amino_acid": alternate_aa,
            "reconciled_protein_position": transcript.get("protein_start", ""),
            "ensembl_protein_id": transcript.get("protein_id", ""),
            "uniprot_swissprot_accessions": accessions(transcript.get("swissprot")),
            "uniprot_trembl_accessions": accessions(transcript.get("trembl")),
            "vep_domains": rendered_domains(transcript.get("domains")),
            "strict_missense_eligible": 1,
        }
    )
    return values


def matching_frequencies(result: dict[str, Any], reference: str, alternate: str) -> dict[str, float] | None:
    merged: dict[str, float] = {}
    for colocated in result.get("colocated_variants", []) or []:
        alleles = str(colocated.get("allele_string", "")).split("/")
        if not alleles or alleles[0].upper() != reference.upper() or alternate.upper() not in [a.upper() for a in alleles[1:]]:
            continue
        frequency_object = colocated.get("frequencies", {}) or {}
        alternate_values = frequency_object.get(alternate) or frequency_object.get(alternate.upper())
        if not isinstance(alternate_values, dict):
            continue
        for name, raw in alternate_values.items():
            value = parse_float(raw)
            if value is not None and ("gnomad" in name.lower()):
                merged[name.lower()] = max(value, merged.get(name.lower(), value))
    return merged or None


def population_values(
    row: dict[str, str], result: dict[str, Any] | None, threshold: float
) -> dict[str, Any]:
    values = {
        "gnomad_exome_af": "",
        "gnomad_genome_af": "",
        "gnomad_global_af": "",
        "gnomad_max_population_af": "",
        "gnomad_frequency_fields": "",
        "gnomad_rarity_threshold": threshold,
        "gnomad_rarity_status": "UNKNOWN",
        "population_lookup_status": "",
        "population_lookup_reason": "",
    }
    if result is None:
        values["population_lookup_status"] = "LOOKUP_FAILED"
        values["population_lookup_reason"] = "No VEP result matched the variant key"
        return values
    frequencies = matching_frequencies(result, row["reference_allele"], row["alternate_allele"])
    if frequencies is None:
        values["population_lookup_status"] = "NOT_REPORTED"
        values["population_lookup_reason"] = (
            "VEP succeeded but returned no matching gnomAD allele frequency; AF remains missing"
        )
        values["gnomad_rarity_status"] = "NOT_REPORTED"
        return values
    exome = frequencies.get("gnomade")
    genome = frequencies.get("gnomadg")
    global_values = [value for value in (exome, genome) if value is not None]
    population_values_found = [
        value
        for name, value in frequencies.items()
        if re.match(r"gnomad[eg]_", name)
    ]
    global_af = max(global_values) if global_values else None
    population_max = max(population_values_found) if population_values_found else global_af
    conservative = [value for value in (global_af, population_max) if value is not None]
    if not conservative:
        values["population_lookup_status"] = "FIELDS_EMPTY"
        values["population_lookup_reason"] = "A frequency object existed without usable gnomAD values"
        return values
    maximum = max(conservative)
    values.update(
        {
            "gnomad_exome_af": "" if exome is None else exome,
            "gnomad_genome_af": "" if genome is None else genome,
            "gnomad_global_af": "" if global_af is None else global_af,
            "gnomad_max_population_af": "" if population_max is None else population_max,
            "gnomad_frequency_fields": json.dumps(frequencies, sort_keys=True),
            "gnomad_rarity_status": "RARE" if maximum <= threshold else "NOT_RARE",
            "population_lookup_status": "FOUND",
        }
    )
    return values


def load_jsonl(path: Path) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            value = json.loads(line)
            if isinstance(value, list):
                objects = value
            else:
                objects = [value]
            for result in objects:
                key = input_key(result)
                if key:
                    results[key] = result
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/03_candidate_pool.csv"))
    parser.add_argument("--vep-jsonl", default=str(PROJECT_ROOT / "data/raw/vep/05_vep.jsonl"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/06_vep_reconciled.csv"))
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    args = parser.parse_args()
    config = load_config(args.config)
    threshold = float(config["rarity"]["maximum_global_af"])
    input_path = Path(args.input)
    vep_path = Path(args.vep_jsonl)
    output_path = Path(args.output)
    rows = read_csv(input_path)
    results = load_jsonl(vep_path)
    enriched = []
    for row in rows:
        result = results.get(row["variant_key"])
        combined = dict(row)
        combined.update(transcript_values(row, result))
        combined.update(population_values(row, result, threshold))
        enriched.append(combined)
    fieldnames = list(rows[0]) + [column for column in ADDED_COLUMNS if column not in rows[0]]
    write_csv(output_path, enriched, fieldnames)
    audit = {
        **audit_basic_rows(enriched),
        "vep_results_by_key": len(results),
        "transcript_status_counts": dict(Counter(row["transcript_reconciliation_status"] for row in enriched)),
        "population_lookup_counts": dict(Counter(row["population_lookup_status"] for row in enriched)),
        "rarity_status_counts": dict(Counter(row["gnomad_rarity_status"] for row in enriched)),
        "rarity_threshold": threshold,
        "not_reported_policy": (
            "Kept distinct from AF=0. Eligibility is decided later from configuration; lookup "
            "failures remain different from successfully queried but unreported variants."
        ),
    }
    audit_path = PROJECT_ROOT / "audits/06_vep_reconciliation_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "06_parse_vep",
        [input_path, vep_path],
        [output_path, audit_path],
        {"rarity_threshold": threshold, "vep_release": config["vep"]["release"]},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
