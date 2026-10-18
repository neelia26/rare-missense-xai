#!/usr/bin/env python3
"""Retrieve UniProt records once per accession and validate the reference residue."""

from __future__ import annotations

import argparse
import json
import random
import re
import sys
import time
import urllib.error
import urllib.parse
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
    split_multi_value,
    write_csv,
    write_json,
)


ADDED_COLUMNS = [
    "uniprot_selected_accession",
    "uniprot_accession_source",
    "uniprot_lookup_status",
    "uniprot_primary_accession",
    "uniprot_reviewed",
    "uniprot_release",
    "uniprot_entry_version",
    "uniprot_sequence_version",
    "uniprot_protein_name",
    "uniprot_gene_name",
    "uniprot_protein_length",
    "uniprot_reference_residue",
    "uniprot_validation_status",
    "uniprot_overlapping_feature_types",
    "uniprot_overlapping_feature_descriptions",
]

# Only these outcomes are stable enough to cache. Transient failures (ERROR) and
# empty/stub bodies (EMPTY_RECORD) are retried on the next run.
CACHEABLE_STATUSES = {"FOUND", "NOT_FOUND"}


def safe_filename(accession: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]", "_", accession)


def bare_accession(value: str) -> str:
    """Strip a trailing version suffix (e.g. P01130.272 -> P01130).

    VEP's SWISSPROT/TREMBL fields append the entry version. UniProt's REST API treats
    the dotted form as a version lookup and returns a stub with no sequence.
    """
    return re.sub(r"\.\d+$", "", value.strip())


def select_accession(row: dict[str, str], allow_trembl: bool) -> tuple[str, str]:
    swissprot = split_multi_value(row.get("uniprot_swissprot_accessions"))
    if swissprot:
        return bare_accession(swissprot[0]), "SWISSPROT"
    trembl = split_multi_value(row.get("uniprot_trembl_accessions"))
    if allow_trembl and trembl:
        return bare_accession(trembl[0]), "TREMBL_FALLBACK"
    return "", "NO_ELIGIBLE_ACCESSION"


def retry_delay(error: urllib.error.HTTPError, attempt: int) -> float:
    try:
        retry_after = float(error.headers.get("Retry-After", 0) or 0)
    except (TypeError, ValueError):
        retry_after = 0.0
    if retry_after > 0:
        return min(60.0, retry_after)
    return min(60.0, 2**attempt + random.random())


def is_usable_record(record: dict[str, Any] | None) -> bool:
    """True for a record that is safe to reuse from cache."""
    if not record:
        return False
    status = record.get("status")
    if status == "NOT_FOUND":
        return True
    if status != "FOUND":
        return False
    payload = record.get("payload") or {}
    return bool(payload.get("primaryAccession")) and bool((payload.get("sequence") or {}).get("value"))


def fetch_record(accession: str, retries: int) -> dict[str, Any]:
    url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(accession)}.json"
    last_error = ""
    for attempt in range(retries):
        request = urllib.request.Request(
            url,
            headers={"Accept": "application/json", "User-Agent": "rare-missense-rebuild/1.0"},
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                release = response.headers.get("x-uniprot-release", "")
                payload = json.load(response)
            sequence = (payload.get("sequence") or {}).get("value")
            if not payload.get("primaryAccession") or not sequence:
                return {"status": "EMPTY_RECORD", "release": release, "payload": {}}
            return {"status": "FOUND", "release": release, "payload": payload}
        except urllib.error.HTTPError as error:
            if error.code == 404:
                return {"status": "NOT_FOUND", "release": "", "payload": {}}
            last_error = f"HTTP {error.code}: {error.reason}"
            if error.code not in {429, 500, 502, 503, 504}:
                break
            time.sleep(retry_delay(error, attempt))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as error:
            last_error = str(error)
            time.sleep(min(60.0, 2**attempt + random.random()))
    return {"status": "ERROR", "release": "", "error": last_error, "payload": {}}


def position_value(value: Any) -> int | None:
    if isinstance(value, dict):
        value = value.get("value")
    return parse_int(value)


def overlap_features(payload: dict[str, Any], position: int) -> tuple[str, str]:
    types: list[str] = []
    descriptions: list[str] = []
    for feature in payload.get("features", []) or []:
        location = feature.get("location", {}) or {}
        start = position_value(location.get("start"))
        end = position_value(location.get("end"))
        if start is None or end is None or not (start <= position <= end):
            continue
        feature_type = str(feature.get("type", ""))
        description = str(feature.get("description", ""))
        if feature_type and feature_type not in types:
            types.append(feature_type)
        if description and description not in descriptions:
            descriptions.append(description)
    return ";".join(types), ";".join(descriptions)


def protein_name(payload: dict[str, Any]) -> str:
    description = payload.get("proteinDescription", {}) or {}
    recommended = description.get("recommendedName", {}) or {}
    full = recommended.get("fullName", {}) or {}
    if full.get("value"):
        return str(full["value"])
    submission = description.get("submissionNames", []) or []
    if submission:
        return str((submission[0].get("fullName") or {}).get("value", ""))
    return ""


def gene_name(payload: dict[str, Any]) -> str:
    genes = payload.get("genes", []) or []
    if not genes:
        return ""
    return str((genes[0].get("geneName") or {}).get("value", ""))


def enrich_row(
    row: dict[str, str], accession: str, source: str, record: dict[str, Any] | None
) -> dict[str, Any]:
    values = {column: "" for column in ADDED_COLUMNS}
    values.update(
        {
            "uniprot_selected_accession": accession,
            "uniprot_accession_source": source,
            "uniprot_lookup_status": "NO_ACCESSION" if not accession else "",
            "uniprot_validation_status": "NO_ACCESSION" if not accession else "",
        }
    )
    if not accession:
        return values
    if not record or record.get("status") != "FOUND":
        status = str((record or {}).get("status", "ERROR"))
        values["uniprot_lookup_status"] = status
        values["uniprot_validation_status"] = status
        return values
    payload = record["payload"]
    sequence = str((payload.get("sequence") or {}).get("value", ""))
    position = parse_int(row.get("reconciled_protein_position"))
    expected = row.get("reconciled_reference_amino_acid", "")
    if position is None:
        observed = ""
        validation = "POSITION_MISSING"
    elif position < 1 or position > len(sequence):
        observed = ""
        validation = "POSITION_OUT_OF_RANGE"
    else:
        observed = sequence[position - 1]
        validation = "MATCH" if observed == expected else "REFERENCE_AA_MISMATCH"
    feature_types, feature_descriptions = (
        overlap_features(payload, position) if position is not None else ("", "")
    )
    audit = payload.get("entryAudit", {}) or {}
    values.update(
        {
            "uniprot_lookup_status": "FOUND",
            "uniprot_primary_accession": payload.get("primaryAccession", ""),
            "uniprot_reviewed": int(
                "reviewed" in str(payload.get("entryType", "")).lower()
                and "unreviewed" not in str(payload.get("entryType", "")).lower()
            ),
            "uniprot_release": record.get("release", ""),
            "uniprot_entry_version": audit.get("entryVersion", ""),
            "uniprot_sequence_version": audit.get("sequenceVersion", ""),
            "uniprot_protein_name": protein_name(payload),
            "uniprot_gene_name": gene_name(payload),
            "uniprot_protein_length": len(sequence),
            "uniprot_reference_residue": observed,
            "uniprot_validation_status": validation,
            "uniprot_overlapping_feature_types": feature_types,
            "uniprot_overlapping_feature_descriptions": feature_descriptions,
        }
    )
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/06_vep_reconciled.csv"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/07_uniprot_validated.csv"))
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    args = parser.parse_args()
    config = load_config(args.config)
    settings = config["uniprot"]
    input_path = Path(args.input)
    output_path = Path(args.output)
    cache_dir = PROJECT_ROOT / "data/raw/uniprot/cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    rows = read_csv(input_path)
    selected = [
        select_accession(row, bool(settings["allow_trembl_fallback"])) for row in rows
    ]
    unique_accessions = sorted({accession for accession, _ in selected if accession})

    records: dict[str, dict[str, Any]] = {}
    pending: list[str] = []
    for accession in unique_accessions:
        checkpoint = cache_dir / f"{safe_filename(accession)}.json"
        cached: dict[str, Any] | None = None
        if checkpoint.exists():
            try:
                cached = json.loads(checkpoint.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                cached = None
        if is_usable_record(cached):
            records[accession] = cached
        else:
            pending.append(accession)
    with ThreadPoolExecutor(max_workers=int(settings["workers"])) as executor:
        futures = {
            executor.submit(fetch_record, accession, int(settings["max_retries"])): accession
            for accession in pending
        }
        for complete, future in enumerate(as_completed(futures), start=1):
            accession = futures[future]
            record = future.result()
            records[accession] = record
            if record.get("status") in CACHEABLE_STATUSES:
                write_json(cache_dir / f"{safe_filename(accession)}.json", record)
            print(f"UniProt {complete}/{len(pending)} completed", flush=True)

    enriched = []
    for row, (accession, source) in zip(rows, selected):
        combined = dict(row)
        combined.update(enrich_row(row, accession, source, records.get(accession)))
        enriched.append(combined)
    fieldnames = list(rows[0]) + [column for column in ADDED_COLUMNS if column not in rows[0]]

    # If nothing validated, write the table next to the real output under a .FAILED name
    # so it can be inspected without the stage looking DONE in `run_pipeline.py status`.
    has_match = any(row["uniprot_validation_status"] == "MATCH" for row in enriched)
    written_path = output_path if has_match else output_path.with_name(output_path.stem + ".FAILED.csv")
    write_csv(written_path, enriched, fieldnames)
    audit = {
        **audit_basic_rows(enriched),
        "unique_accessions": len(unique_accessions),
        "cache_records_reused": len(unique_accessions) - len(pending),
        "records_requested": len(pending),
        "lookup_status_counts": dict(Counter(row["uniprot_lookup_status"] for row in enriched)),
        "validation_status_counts": dict(Counter(row["uniprot_validation_status"] for row in enriched)),
        "rule": "Prefer VEP Swiss-Prot; validate the VEP reference amino acid against the selected UniProt sequence.",
    }
    audit_path = PROJECT_ROOT / "audits/07_uniprot_audit.json"
    write_json(audit_path, audit)
    print(json.dumps(audit, indent=2))
    if not has_match:
        sys.exit(
            "No MATCH rows after UniProt validation. Inspect "
            f"{audit_path} and {written_path}; stage not marked complete."
        )
    append_stage_manifest(
        "07_enrich_uniprot",
        [input_path],
        [output_path, audit_path],
        {"allow_trembl_fallback": settings["allow_trembl_fallback"]},
    )


if __name__ == "__main__":
    main()