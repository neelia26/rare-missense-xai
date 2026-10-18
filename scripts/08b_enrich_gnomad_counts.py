#!/usr/bin/env python3
"""Optionally add direct gnomAD AC/AN/homozygote counts from local indexed site VCFs."""

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
    normalize_chromosome,
    parse_float,
    read_csv,
    write_csv,
    write_json,
)


def allele_value(value: Any, index: int) -> Any:
    if isinstance(value, (tuple, list)):
        return value[index] if index < len(value) else ""
    return value if index == 0 else ""


def query(vcf: Any, row: dict[str, str]) -> dict[str, Any]:
    chromosome = normalize_chromosome(row["chromosome"])
    contig = chromosome if chromosome in vcf.header.contigs else f"chr{chromosome}"
    if contig not in vcf.header.contigs:
        return {"status": "CONTIG_NOT_FOUND"}
    position = int(row["position"])
    try:
        records = vcf.fetch(contig, position - 1, position)
    except (ValueError, OSError):
        return {"status": "LOOKUP_FAILED"}
    for record in records:
        if record.pos != position or record.ref.upper() != row["reference_allele"].upper():
            continue
        alternates = list(record.alts or [])
        for index, alternate in enumerate(alternates):
            if alternate.upper() != row["alternate_allele"].upper():
                continue
            return {
                "status": "FOUND",
                "af": allele_value(record.info.get("AF", ""), index),
                "ac": allele_value(record.info.get("AC", ""), index),
                "an": record.info.get("AN", ""),
                "nhomalt": allele_value(
                    record.info.get("nhomalt", record.info.get("nhom_alt", "")), index
                ),
            }
    return {"status": "NOT_REPORTED"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/08_interpro.csv"))
    parser.add_argument("--exomes-vcf")
    parser.add_argument("--genomes-vcf")
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--rarity-threshold", type=float, default=0.01)
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/08b_gnomad_counts.csv"))
    args = parser.parse_args()
    if not args.exomes_vcf and not args.genomes_vcf:
        raise SystemExit("Provide --exomes-vcf, --genomes-vcf, or both")
    try:
        import pysam  # type: ignore
    except ImportError as error:
        raise RuntimeError("Install requirements.txt; pysam is required for indexed VCF lookup") from error

    input_path = Path(args.input)
    output_path = Path(args.output)
    rows = read_csv(input_path)
    sources: dict[str, Any] = {}
    source_paths: list[Path] = []
    for name, value in (("exome", args.exomes_vcf), ("genome", args.genomes_vcf)):
        if value:
            path = Path(value).expanduser().resolve()
            sources[name] = pysam.VariantFile(str(path))
            source_paths.append(path)

    added = ["gnomad_direct_release"]
    for name in ("exome", "genome"):
        added.extend(
            [
                f"gnomad_direct_{name}_status",
                f"gnomad_direct_{name}_af",
                f"gnomad_direct_{name}_ac",
                f"gnomad_direct_{name}_an",
                f"gnomad_direct_{name}_nhomalt",
            ]
        )
    added.extend(["gnomad_direct_max_af", "gnomad_direct_rarity_status"])

    enriched = []
    for row in rows:
        values: dict[str, Any] = {column: "" for column in added}
        values["gnomad_direct_release"] = args.source_release
        results: dict[str, dict[str, Any]] = {}
        for name in ("exome", "genome"):
            result = query(sources[name], row) if name in sources else {"status": "FILE_NOT_PROVIDED"}
            results[name] = result
            values[f"gnomad_direct_{name}_status"] = result["status"]
            for field in ("af", "ac", "an", "nhomalt"):
                values[f"gnomad_direct_{name}_{field}"] = result.get(field, "")
        frequencies = [
            parse_float(result.get("af"))
            for result in results.values()
            if result.get("status") == "FOUND"
        ]
        frequencies = [value for value in frequencies if value is not None]
        if frequencies:
            maximum = max(frequencies)
            rarity = "RARE" if maximum <= args.rarity_threshold else "NOT_RARE"
            values["gnomad_direct_max_af"] = maximum
        elif any(result["status"] in {"LOOKUP_FAILED", "CONTIG_NOT_FOUND"} for result in results.values()):
            rarity = "LOOKUP_FAILED"
        else:
            rarity = "NOT_REPORTED"
        values["gnomad_direct_rarity_status"] = rarity
        combined = dict(row)
        combined.update(values)
        enriched.append(combined)

    write_csv(output_path, enriched, list(rows[0]) + [column for column in added if column not in rows[0]])
    audit = {
        **audit_basic_rows(enriched),
        "source_release": args.source_release,
        "rarity_threshold": args.rarity_threshold,
        "direct_rarity_status_counts": dict(Counter(row["gnomad_direct_rarity_status"] for row in enriched)),
        "important_note": (
            "NOT_REPORTED is not converted to AF=0. Coverage-aware interpretation requires an "
            "appropriate gnomAD coverage resource and remains separate from lookup failure."
        ),
    }
    audit_path = PROJECT_ROOT / "audits/08b_gnomad_counts_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "08b_enrich_gnomad_counts",
        [input_path, *source_paths],
        [output_path, audit_path],
        {"source_release": args.source_release, "rarity_threshold": args.rarity_threshold},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

