#!/usr/bin/env python3
"""Retrieve InterPro entries for validated proteins and test residue overlap."""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
import time
import urllib.error
import urllib.request
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    load_config,
    parse_int,
    read_csv,
    write_csv,
    write_json,
)


ADDED_COLUMNS = [
    "interpro_lookup_status",
    "interpro_release",
    "interpro_entry_count",
    "interpro_overlap_count",
    "interpro_overlap_accessions",
    "interpro_overlap_names",
    "interpro_overlap_types",
]


def safe_filename(accession: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", accession)


def fetch_page(url: str, retries: int) -> tuple[dict[str, Any] | None, str, str]:
    last_error = ""
    for attempt in range(retries):
        request = urllib.request.Request(
            url,
            headers={"Accept": "application/json", "User-Agent": "rare-missense-rebuild/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                release = response.headers.get("InterPro-Version", "")
                return json.load(response), release, ""
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return {"results": [], "next": None}, "", ""
            last_error = f"HTTP {error.code}: {error.reason}"
            if error.code not in {429, 500, 502, 503, 504}:
                break
            time.sleep(min(60.0, 2**attempt + random.random()))
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = str(error)
            time.sleep(min(60.0, 2**attempt + random.random()))
    return None, "", last_error


def fetch_accession(accession: str, retries: int) -> dict[str, Any]:
    url = (
        "https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/"
        f"{accession}/?page_size=200"
    )
    results: list[dict[str, Any]] = []
    release = ""
    pages = 0
    while url:
        pages += 1
        if pages > 100:
            return {"status": "ERROR", "error": "pagination exceeded 100 pages", "results": []}
        page, page_release, error = fetch_page(url, retries)
        if page is None:
            return {"status": "ERROR", "error": error, "results": []}
        release = release or page_release
        results.extend(page.get("results", []) or [])
        url = str(page.get("next") or "")
    return {"status": "FOUND", "release": release, "pages": pages, "results": results}


def fragments(entry: dict[str, Any]) -> list[tuple[int, int]]:
    output: list[tuple[int, int]] = []
    for protein in entry.get("proteins", []) or []:
        for location in protein.get("entry_protein_locations", []) or []:
            for fragment in location.get("fragments", []) or []:
                start = parse_int(fragment.get("start"))
                end = parse_int(fragment.get("end"))
                if start is not None and end is not None:
                    output.append((start, end))
    return output


def values_for_row(row: dict[str, str], record: dict[str, Any] | None) -> dict[str, Any]:
    values = {column: "" for column in ADDED_COLUMNS}
    if row.get("uniprot_validation_status") != "MATCH":
        values["interpro_lookup_status"] = "SKIPPED_UNVALIDATED_PROTEIN"
        return values
    if not record or record.get("status") != "FOUND":
        values["interpro_lookup_status"] = str((record or {}).get("status", "ERROR"))
        return values
    position = parse_int(row.get("reconciled_protein_position"))
    overlaps: list[dict[str, str]] = []
    for entry in record.get("results", []):
        if position is None or not any(start <= position <= end for start, end in fragments(entry)):
            continue
        metadata = entry.get("metadata", {}) or {}
        overlaps.append(
            {
                "accession": str(metadata.get("accession", "")),
                "name": str(metadata.get("name", "")),
                "type": str(metadata.get("type", "")),
            }
        )
    values.update(
        {
            "interpro_lookup_status": "FOUND",
            "interpro_release": record.get("release", ""),
            "interpro_entry_count": len(record.get("results", [])),
            "interpro_overlap_count": len(overlaps),
            "interpro_overlap_accessions": ";".join(item["accession"] for item in overlaps if item["accession"]),
            "interpro_overlap_names": ";".join(item["name"] for item in overlaps if item["name"]),
            "interpro_overlap_types": ";".join(item["type"] for item in overlaps if item["type"]),
        }
    )
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/07_uniprot_validated.csv"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/08_interpro.csv"))
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    args = parser.parse_args()
    config = load_config(args.config)
    settings = config["interpro"]
    input_path = Path(args.input)
    output_path = Path(args.output)
    cache_dir = PROJECT_ROOT / "data/raw/interpro/cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    rows = read_csv(input_path)
    accessions = sorted(
        {
            row["uniprot_selected_accession"]
            for row in rows
            if row.get("uniprot_validation_status") == "MATCH"
            and row.get("uniprot_selected_accession")
        }
    )
    records: dict[str, dict[str, Any]] = {}
    pending: list[str] = []
    for accession in accessions:
        checkpoint = cache_dir / f"{safe_filename(accession)}.json"
        if checkpoint.exists():
            records[accession] = json.loads(checkpoint.read_text(encoding="utf-8"))
        else:
            pending.append(accession)
    with ThreadPoolExecutor(max_workers=int(settings["workers"])) as executor:
        futures = {
            executor.submit(fetch_accession, accession, int(settings["max_retries"])): accession
            for accession in pending
        }
        for complete, future in enumerate(as_completed(futures), start=1):
            accession = futures[future]
            record = future.result()
            records[accession] = record
            write_json(cache_dir / f"{safe_filename(accession)}.json", record)
            print(f"InterPro {complete}/{len(pending)} completed", flush=True)
    enriched = []
    for row in rows:
        combined = dict(row)
        combined.update(values_for_row(row, records.get(row.get("uniprot_selected_accession", ""))))
        enriched.append(combined)
    fieldnames = list(rows[0]) + [column for column in ADDED_COLUMNS if column not in rows[0]]
    write_csv(output_path, enriched, fieldnames)
    audit = {
        **audit_basic_rows(enriched),
        "unique_validated_accessions": len(accessions),
        "cache_records_reused": len(accessions) - len(pending),
        "records_requested": len(pending),
        "lookup_status_counts": dict(Counter(row["interpro_lookup_status"] for row in enriched)),
        "variants_with_interpro_overlap": sum(int(row.get("interpro_overlap_count") or 0) > 0 for row in enriched),
    }
    audit_path = PROJECT_ROOT / "audits/08_interpro_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "08_enrich_interpro",
        [input_path],
        [output_path, audit_path],
        {"api": "InterPro protein entry endpoint", "workers": settings["workers"]},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
