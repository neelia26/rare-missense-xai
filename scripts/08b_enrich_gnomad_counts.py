#!/usr/bin/env python3
"""Optionally add direct gnomAD AC/AN/homozygote counts from indexed site VCFs.

Each of --exomes-vcf / --genomes-vcf may be:
  * a local indexed .vcf.bgz (index .tbi/.csi next to it),
  * a remote https/gs/s3 URL of an indexed .vcf.bgz (index expected at <url>.tbi),
  * a per-chromosome template containing {chrom}, e.g.
        /data/gnomad.exomes.v4.1.sites.{chrom}.vcf.bgz
    where {chrom} is replaced by the contig name with a "chr" prefix (chr1, chrX, ...).
    Only the chromosomes your variants fall on are ever opened.
"""

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


REMOTE_PREFIXES = ("http://", "https://", "gs://", "s3://", "ftp://")


def is_remote(value: str) -> bool:
    return value.lower().startswith(REMOTE_PREFIXES)


def allele_value(value: Any, index: int) -> Any:
    if isinstance(value, (tuple, list)):
        return value[index] if index < len(value) else ""
    return value if index == 0 else ""


def info_value(record: Any, *keys: str) -> Any:
    """First INFO field present in the record; pysam raises ValueError for keys the header lacks."""
    for key in keys:
        try:
            value = record.info.get(key)
        except (ValueError, KeyError):
            continue
        if value is not None:
            return value
    return ""


def bare_chromosome(value: str) -> str:
    chromosome = normalize_chromosome(value)
    return chromosome[3:] if chromosome.lower().startswith("chr") else chromosome


class VcfSource:
    """Lazily opens one indexed VCF, or one per chromosome when given a {chrom} template."""

    def __init__(self, name: str, template: str, pysam: Any) -> None:
        self.name = name
        self.template = template
        self.pysam = pysam
        self.per_chromosome = "{chrom}" in template
        self.handles: dict[str, Any] = {}
        self.errors: dict[str, str] = {}

    def location(self, chrom: str) -> str:
        value = self.template.replace("{chrom}", chrom)
        return value if is_remote(value) else str(Path(value).expanduser().resolve())

    def get(self, chrom: str) -> Any:
        key = chrom if self.per_chromosome else ""
        if key in self.handles:
            return self.handles[key]
        location = self.location(chrom)
        try:
            if is_remote(location):
                handle = self.pysam.VariantFile(location, index_filename=location + ".tbi")
            else:
                handle = self.pysam.VariantFile(location)
        except (OSError, ValueError) as error:
            self.errors[key or location] = f"{location}: {error}"
            handle = None
        self.handles[key] = handle  # failures are cached so they are not retried per row
        return handle


def query(vcf: Any, row: dict[str, str], bare: str) -> dict[str, Any]:
    contig = next((c for c in (bare, f"chr{bare}") if c in vcf.header.contigs), None)
    if contig is None:
        return {"status": "CONTIG_NOT_FOUND"}
    position = int(row["position"])
    try:
        records = list(vcf.fetch(contig, position - 1, position))
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
                "af": allele_value(info_value(record, "AF"), index),
                "ac": allele_value(info_value(record, "AC"), index),
                "an": info_value(record, "AN"),
                "nhomalt": allele_value(info_value(record, "nhomalt", "nhom_alt"), index),
            }
    return {"status": "NOT_REPORTED"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/08_interpro.csv"))
    parser.add_argument("--exomes-vcf", help="indexed VCF path/URL, or a template containing {chrom}")
    parser.add_argument("--genomes-vcf", help="indexed VCF path/URL, or a template containing {chrom}")
    parser.add_argument("--source-release", required=True)
    parser.add_argument("--rarity-threshold", type=float, default=0.01)
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/08b_gnomad_counts.csv"))
    parser.add_argument(
        "--smoke-test",
        type=int,
        default=0,
        metavar="N",
        help="query only the first N rows, print status counts and errors, write nothing",
    )
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
    sources: dict[str, VcfSource] = {}
    source_paths: list[Path] = []
    for name, value in (("exome", args.exomes_vcf), ("genome", args.genomes_vcf)):
        if value:
            sources[name] = VcfSource(name, value, pysam)
            if "{chrom}" not in value and not is_remote(value):
                source_paths.append(Path(value).expanduser().resolve())

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

    target_rows = rows[: args.smoke_test] if args.smoke_test else rows
    enriched = []
    for row in target_rows:
        values: dict[str, Any] = {column: "" for column in added}
        values["gnomad_direct_release"] = args.source_release
        bare = bare_chromosome(row["chromosome"])
        results: dict[str, dict[str, Any]] = {}
        for name in ("exome", "genome"):
            if name not in sources:
                result: dict[str, Any] = {"status": "FILE_NOT_PROVIDED"}
            else:
                vcf = sources[name].get(f"chr{bare}")
                result = query(vcf, row, bare) if vcf is not None else {"status": "LOOKUP_FAILED"}
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

    exome_counts = dict(Counter(row["gnomad_direct_exome_status"] for row in enriched))
    genome_counts = dict(Counter(row["gnomad_direct_genome_status"] for row in enriched))
    rarity_counts = dict(Counter(row["gnomad_direct_rarity_status"] for row in enriched))
    open_errors = {name: dict(list(source.errors.items())[:5]) for name, source in sources.items() if source.errors}

    if args.smoke_test:
        print(json.dumps(
            {
                "smoke_test_rows": len(enriched),
                "exome_status_counts": exome_counts,
                "genome_status_counts": genome_counts,
                "rarity_status_counts": rarity_counts,
                "open_errors": open_errors,
            },
            indent=2,
        ))
        print("Smoke test only: nothing was written.")
        return

    # If every lookup failed, keep the table for inspection under a .FAILED name so the
    # stage does not look DONE in `run_pipeline.py status`.
    all_failed = bool(enriched) and rarity_counts.get("LOOKUP_FAILED", 0) == len(enriched)
    written_path = output_path.with_name(output_path.stem + ".FAILED.csv") if all_failed else output_path
    write_csv(written_path, enriched, list(rows[0]) + [column for column in added if column not in rows[0]])
    audit = {
        **audit_basic_rows(enriched),
        "source_release": args.source_release,
        "rarity_threshold": args.rarity_threshold,
        "exomes_vcf": args.exomes_vcf or "",
        "genomes_vcf": args.genomes_vcf or "",
        "direct_rarity_status_counts": rarity_counts,
        "exome_status_counts": exome_counts,
        "genome_status_counts": genome_counts,
        "open_errors": open_errors,
        "important_note": (
            "NOT_REPORTED is not converted to AF=0. Coverage-aware interpretation requires an "
            "appropriate gnomAD coverage resource and remains separate from lookup failure."
        ),
    }
    audit_path = PROJECT_ROOT / "audits/08b_gnomad_counts_audit.json"
    write_json(audit_path, audit)
    print(json.dumps(audit, indent=2))
    if all_failed:
        sys.exit(
            "Every gnomAD lookup failed; check contig names, genome build, paths/URLs and index files. "
            f"See {audit_path}; stage not marked complete."
        )
    append_stage_manifest(
        "08b_enrich_gnomad_counts",
        [input_path, *source_paths],
        [output_path, audit_path],
        {
            "source_release": args.source_release,
            "rarity_threshold": args.rarity_threshold,
            "exomes_vcf": args.exomes_vcf or "",
            "genomes_vcf": args.genomes_vcf or "",
        },
    )


if __name__ == "__main__":
    main()