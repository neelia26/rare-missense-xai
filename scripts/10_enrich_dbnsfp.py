#!/usr/bin/env python3
"""Retrieve dbNSFP annotations through MyVariant.info for GRCh38 variants.

The complete dbNSFP object is requested because the MyVariant field filter does
not reliably return keys containing ``++`` (notably ``gerp++``). Raw API
batches and service metadata are cached so enrichment can be audited and
resumed without downloading the full dbNSFP archive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path
from typing import Any, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    load_config,
    normalize_chromosome,
    parse_float,
    parse_int,
    read_csv,
    utc_now,
    write_csv,
    write_json,
)


ACCESS_SERVICE = "MyVariant.info"
DEFAULT_API_ENDPOINT = "https://myvariant.info/v1/variant"
DEFAULT_METADATA_ENDPOINT = "https://myvariant.info/v1/metadata"
REQUEST_FIELDS = "dbnsfp"

CONSERVATION_PATHS = {
    "GERP++_RS": (("gerp++", "rs"), "conservation_gerp_rs"),
    "phyloP100way_vertebrate": (
        ("phylop", "100way_vertebrate", "score"),
        "conservation_phylop100way",
    ),
    "phastCons100way_vertebrate": (
        ("phastcons", "100way_vertebrate", "score"),
        "conservation_phastcons100way",
    ),
}

OPTIONAL_PREDICTOR_PATHS = {
    "SIFT_score": ("sift", "score"),
    "Polyphen2_HDIV_score": ("polyphen2", "hdiv", "score"),
    "REVEL_score": ("revel", "score"),
    "CADD_phred": ("cadd", "phred"),
    "AlphaMissense_score": ("alphamissense", "score"),
}

PROVENANCE_COLUMNS = [
    "dbnsfp_release",
    "dbnsfp_access_service",
    "dbnsfp_source_url",
    "dbnsfp_metadata_url",
    "myvariant_build_version",
    "myvariant_build_date",
    "dbnsfp_query_id",
    "dbnsfp_record_id",
    "dbnsfp_coordinate_verified",
    "dbnsfp_lookup_status",
]


def chunks(values: list[str], size: int) -> Iterable[list[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def nested_value(value: Any, path: tuple[str, ...]) -> Any:
    current = value
    for key in path:
        if not isinstance(current, dict):
            return None
        current = current.get(key)
    return current


def numeric_values(value: Any) -> list[float]:
    if isinstance(value, (list, tuple)):
        output: list[float] = []
        for item in value:
            output.extend(numeric_values(item))
        return output
    if isinstance(value, dict) or isinstance(value, bool):
        return []
    if isinstance(value, str):
        output = []
        for item in re.split(r"[;,&|]", value):
            number = parse_float(item)
            if number is not None:
                output.append(number)
        return output
    number = parse_float(value)
    return [] if number is None else [number]


def numeric_max(values: Iterable[Any]) -> str:
    parsed: list[float] = []
    for value in values:
        parsed.extend(numeric_values(value))
    return "" if not parsed else format(max(parsed), ".15g")


def predictor_column(field: str) -> str:
    return "predictor_" + field.lower().replace("++", "pp").replace("-", "_")


def query_id(row: dict[str, str]) -> str:
    chromosome = normalize_chromosome(row["chromosome"])
    position = int(row["position"])
    reference = row["reference_allele"].strip().upper()
    alternate = row["alternate_allele"].strip().upper()
    if len(reference) != 1 or len(alternate) != 1:
        raise ValueError(
            f"Stage 10 accepts SNVs only; received {row.get('variant_key', '')}: "
            f"{reference}>{alternate}"
        )
    if reference not in "ACGT" or alternate not in "ACGT":
        raise ValueError(
            f"Stage 10 requires A/C/G/T alleles; received {row.get('variant_key', '')}: "
            f"{reference}>{alternate}"
        )
    return f"chr{chromosome}:g.{position}{reference}>{alternate}"


def one_or_many(value: Any) -> list[Any]:
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def contains_text(value: Any, expected: str) -> bool:
    return any(str(item).upper() == expected.upper() for item in one_or_many(value))


def contains_int(value: Any, expected: int) -> bool:
    return any(parse_int(item) == expected for item in one_or_many(value))


def matching_dbnsfp_records(
    row: dict[str, str], response: dict[str, Any]
) -> list[dict[str, Any]]:
    chromosome = normalize_chromosome(row["chromosome"])
    position = int(row["position"])
    reference = row["reference_allele"].strip().upper()
    alternate = row["alternate_allele"].strip().upper()
    matches: list[dict[str, Any]] = []
    for record in one_or_many(response.get("dbnsfp")):
        if not isinstance(record, dict):
            continue
        locations = [
            item for item in one_or_many(record.get("hg38")) if isinstance(item, dict)
        ]
        position_matches = any(
            contains_int(location.get("start"), position)
            and (
                location.get("end") is None
                or contains_int(location.get("end"), position)
            )
            for location in locations
        )
        if (
            contains_text(record.get("chrom"), chromosome)
            and position_matches
            and contains_text(record.get("ref"), reference)
            and contains_text(record.get("alt"), alternate)
        ):
            matches.append(record)
    return matches


def extracted_rows(
    row: dict[str, str],
    response: dict[str, Any],
    source_release: str,
    api_endpoint: str,
    metadata_endpoint: str,
    build_version: str,
    build_date: str,
    conservation_fields: list[str],
    predictor_fields: list[str],
) -> tuple[dict[str, Any], dict[str, Any]]:
    requested_id = query_id(row)
    core: dict[str, Any] = dict(row)
    core.update(
        {
            "dbnsfp_release": source_release,
            "dbnsfp_access_service": ACCESS_SERVICE,
            "dbnsfp_source_url": api_endpoint,
            "dbnsfp_metadata_url": metadata_endpoint,
            "myvariant_build_version": build_version,
            "myvariant_build_date": build_date,
            "dbnsfp_query_id": requested_id,
            "dbnsfp_record_id": response.get("_id", ""),
            "dbnsfp_coordinate_verified": 0,
            "dbnsfp_lookup_status": "",
        }
    )
    optional: dict[str, Any] = {
        "variant_key": row["variant_key"],
        "dbnsfp_release": source_release,
        "dbnsfp_access_service": ACCESS_SERVICE,
        "dbnsfp_query_id": requested_id,
        "dbnsfp_lookup_status": "",
    }
    for field in conservation_fields:
        core[CONSERVATION_PATHS[field][1]] = ""
    for field in predictor_fields:
        optional[predictor_column(field)] = ""

    if response.get("notfound") is True:
        status = "NOT_FOUND"
        matches: list[dict[str, Any]] = []
    elif not response.get("dbnsfp"):
        status = "DBNSFP_NOT_REPORTED"
        matches = []
    else:
        matches = matching_dbnsfp_records(row, response)
        status = "FOUND" if matches else "COORDINATE_MISMATCH"

    core["dbnsfp_lookup_status"] = status
    optional["dbnsfp_lookup_status"] = status
    if status != "FOUND":
        return core, optional

    core["dbnsfp_coordinate_verified"] = 1
    for field in conservation_fields:
        path, destination = CONSERVATION_PATHS[field]
        core[destination] = numeric_max(nested_value(item, path) for item in matches)
    for field in predictor_fields:
        path = OPTIONAL_PREDICTOR_PATHS[field]
        optional[predictor_column(field)] = numeric_max(
            nested_value(item, path) for item in matches
        )
    return core, optional


def request_json(
    url: str,
    retries: int,
    base_delay: float,
    user_agent: str,
    form: dict[str, str] | None = None,
) -> Any:
    encoded = urllib.parse.urlencode(form).encode("utf-8") if form is not None else None
    headers = {"Accept": "application/json", "User-Agent": user_agent}
    if form is not None:
        headers["Content-Type"] = "application/x-www-form-urlencoded"
    last_error = ""
    for attempt in range(retries):
        request = urllib.request.Request(
            url,
            data=encoded,
            headers=headers,
            method="POST" if form is not None else "GET",
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                return json.load(response)
        except urllib.error.HTTPError as error:
            last_error = f"HTTP {error.code}: {error.reason}"
            if error.code not in {429, 500, 502, 503, 504}:
                break
            retry_after = error.headers.get("Retry-After")
            wait = float(retry_after) if retry_after and retry_after.isdigit() else 0.0
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as error:
            last_error = f"{type(error).__name__}: {error}"
            wait = 0.0
        if attempt + 1 < retries:
            time.sleep(max(wait, min(60.0, base_delay * (2**attempt) + random.random())))
    raise RuntimeError(f"MyVariant.info request failed after {retries} attempts: {last_error}")


def batch_cache_path(
    cache_dir: Path,
    ids: list[str],
    assembly: str,
    source_release: str,
    build_version: str,
) -> Path:
    identity = json.dumps(
        {
            "ids": ids,
            "assembly": assembly,
            "fields": REQUEST_FIELDS,
            "source_release": source_release,
            "myvariant_build_version": build_version,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    digest = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:20]
    return cache_dir / "batches" / f"batch_{digest}.json"


def valid_cached_batch(
    cached: Any,
    ids: list[str],
    assembly: str,
    source_release: str,
    build_version: str,
) -> bool:
    return (
        isinstance(cached, dict)
        and cached.get("query_ids") == ids
        and cached.get("assembly") == assembly
        and cached.get("fields") == REQUEST_FIELDS
        and cached.get("dbnsfp_release") == source_release
        and cached.get("myvariant_build_version") == build_version
        and isinstance(cached.get("response"), list)
    )


def obtain_batch(
    ids: list[str],
    cache_dir: Path,
    api_endpoint: str,
    assembly: str,
    source_release: str,
    build_version: str,
    retries: int,
    base_delay: float,
    user_agent: str,
    email: str,
    refresh_cache: bool,
) -> tuple[list[dict[str, Any]], Path, bool]:
    destination = batch_cache_path(
        cache_dir, ids, assembly, source_release, build_version
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not refresh_cache:
        try:
            cached = json.loads(destination.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            cached = None
        if valid_cached_batch(cached, ids, assembly, source_release, build_version):
            return cached["response"], destination, True

    form = {
        "ids": ",".join(ids),
        "fields": REQUEST_FIELDS,
        "assembly": assembly,
    }
    if email:
        form["email"] = email
    response = request_json(api_endpoint, retries, base_delay, user_agent, form=form)
    if not isinstance(response, list):
        raise RuntimeError(
            f"MyVariant.info returned {type(response).__name__}; expected a list for a batch request"
        )
    cache_record = {
        "schema_version": 1,
        "retrieved_at_utc": utc_now(),
        "endpoint": api_endpoint,
        "assembly": assembly,
        "fields": REQUEST_FIELDS,
        "dbnsfp_release": source_release,
        "myvariant_build_version": build_version,
        "query_ids": ids,
        "response": response,
    }
    write_json(destination, cache_record)
    return response, destination, False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", default=str(PROJECT_ROOT / "data/interim/09_constraint.csv")
    )
    parser.add_argument("--output")
    parser.add_argument("--optional-predictor-output")
    parser.add_argument(
        "--config", default=str(PROJECT_ROOT / "config/pipeline_config.json")
    )
    parser.add_argument("--api-endpoint")
    parser.add_argument("--metadata-endpoint")
    parser.add_argument("--expected-source-release")
    parser.add_argument("--batch-size", type=int)
    parser.add_argument("--max-retries", type=int)
    parser.add_argument("--request-pause-seconds", type=float)
    parser.add_argument("--cache-dir")
    parser.add_argument("--email")
    parser.add_argument("--refresh-cache", action="store_true")
    parser.add_argument(
        "--smoke-test",
        type=int,
        metavar="N",
        help="Annotate only the first N rows and write separate smoke-test outputs.",
    )
    args = parser.parse_args()

    config = load_config(args.config)
    settings = config["dbnsfp"]
    assembly = str(settings.get("assembly", "hg38"))
    if assembly.lower() != "hg38":
        raise ValueError("This governed pipeline requires MyVariant assembly=hg38")
    api_endpoint = args.api_endpoint or settings.get(
        "api_endpoint", DEFAULT_API_ENDPOINT
    )
    metadata_endpoint = args.metadata_endpoint or settings.get(
        "metadata_endpoint", DEFAULT_METADATA_ENDPOINT
    )
    batch_size = args.batch_size or int(settings.get("batch_size", 100))
    retries = args.max_retries or int(settings.get("max_retries", 6))
    base_delay = (
        args.request_pause_seconds
        if args.request_pause_seconds is not None
        else float(settings.get("request_pause_seconds", 1.0))
    )
    if not 1 <= batch_size <= 1000:
        raise ValueError("--batch-size must be between 1 and the API maximum of 1000")
    if retries < 1:
        raise ValueError("--max-retries must be at least 1")
    cache_dir = Path(
        args.cache_dir
        or settings.get("cache_dir", PROJECT_ROOT / "data/raw/myvariant/dbnsfp")
    ).expanduser()
    if not cache_dir.is_absolute():
        cache_dir = PROJECT_ROOT / cache_dir
    cache_dir.mkdir(parents=True, exist_ok=True)
    email = args.email if args.email is not None else str(settings.get("email", ""))
    user_agent = str(
        settings.get("user_agent", "rare-missense-evidence-pipeline/1.1")
    )

    input_path = Path(args.input)
    rows = read_csv(input_path)
    if not rows:
        raise ValueError(f"Stage 10 input has no data rows: {input_path}")
    if args.smoke_test is not None:
        if args.smoke_test < 1:
            raise ValueError("--smoke-test must be at least 1")
        rows = rows[: args.smoke_test]
        output_path = Path(
            args.output
            or PROJECT_ROOT / "data/interim/smoke/10_dbnsfp_smoke.csv"
        )
        optional_path = Path(
            args.optional_predictor_output
            or PROJECT_ROOT
            / "data/interim/smoke/10_optional_predictor_scores_smoke.csv"
        )
        audit_path = PROJECT_ROOT / "audits/10_dbnsfp_smoke_audit.json"
        stage_name = "10_enrich_dbnsfp_smoke"
    else:
        output_path = Path(
            args.output or PROJECT_ROOT / "data/interim/10_dbnsfp.csv"
        )
        optional_path = Path(
            args.optional_predictor_output
            or PROJECT_ROOT / "data/interim/10_optional_predictor_scores.csv"
        )
        audit_path = PROJECT_ROOT / "audits/10_dbnsfp_audit.json"
        stage_name = "10_enrich_dbnsfp"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    optional_path.parent.mkdir(parents=True, exist_ok=True)
    audit_path.parent.mkdir(parents=True, exist_ok=True)

    conservation_fields = list(settings["conservation_fields"])
    predictor_fields = list(settings["optional_predictor_fields"])
    unknown_conservation = [
        field for field in conservation_fields if field not in CONSERVATION_PATHS
    ]
    unknown_predictors = [
        field for field in predictor_fields if field not in OPTIONAL_PREDICTOR_PATHS
    ]
    if unknown_conservation or unknown_predictors:
        raise ValueError(
            "Unmapped MyVariant fields in configuration: "
            f"conservation={unknown_conservation}, predictors={unknown_predictors}"
        )

    metadata = request_json(
        metadata_endpoint,
        retries=retries,
        base_delay=base_delay,
        user_agent=user_agent,
    )
    source_metadata = metadata.get("src", {}).get("dbnsfp", {})
    source_release = str(source_metadata.get("version", ""))
    build_version = str(metadata.get("build_version", ""))
    build_date = str(metadata.get("build_date", ""))
    if not source_release or not build_version:
        raise RuntimeError("MyVariant metadata did not contain dbNSFP and build versions")
    expected_release = args.expected_source_release or str(
        settings.get("expected_release", "")
    )
    if expected_release and source_release != expected_release:
        raise RuntimeError(
            "MyVariant currently exposes dbNSFP "
            f"{source_release}, but {expected_release} was expected. Review the new source "
            "before deliberately updating the configuration."
        )
    metadata_path = cache_dir / f"metadata_{build_version}.json"
    write_json(metadata_path, metadata)

    ids = [query_id(row) for row in rows]
    if len(ids) != len(set(ids)):
        duplicates = [key for key, count in Counter(ids).items() if count > 1]
        raise ValueError(
            "Stage 10 input contains duplicate genomic variants: "
            + ", ".join(duplicates[:5])
        )

    responses: dict[str, dict[str, Any]] = {}
    cache_paths: list[Path] = []
    batches_reused = 0
    api_batches_requested = 0
    batch_list = list(chunks(ids, batch_size))
    for number, batch_ids in enumerate(batch_list, start=1):
        batch_response, cache_path, reused = obtain_batch(
            batch_ids,
            cache_dir,
            api_endpoint,
            assembly,
            source_release,
            build_version,
            retries,
            base_delay,
            user_agent,
            email,
            args.refresh_cache,
        )
        cache_paths.append(cache_path)
        batches_reused += int(reused)
        api_batches_requested += int(not reused)
        batch_queries: set[str] = set()
        for item in batch_response:
            if isinstance(item, dict) and item.get("query"):
                identifier = str(item["query"])
                responses[identifier] = item
                batch_queries.add(identifier)
        missing = [item for item in batch_ids if item not in batch_queries]
        if missing:
            raise RuntimeError(
                "MyVariant batch response omitted query IDs: " + ", ".join(missing[:5])
            )
        print(
            f"MyVariant batch {number}/{len(batch_list)} "
            f"({'cache' if reused else 'API'}) completed",
            flush=True,
        )
        if not reused and number < len(batch_list) and base_delay > 0:
            time.sleep(base_delay)

    core_rows: list[dict[str, Any]] = []
    optional_rows: list[dict[str, Any]] = []
    for row, identifier in zip(rows, ids):
        core, optional = extracted_rows(
            row,
            responses[identifier],
            source_release,
            api_endpoint,
            metadata_endpoint,
            build_version,
            build_date,
            conservation_fields,
            predictor_fields,
        )
        core_rows.append(core)
        optional_rows.append(optional)

    conservation_columns = [
        CONSERVATION_PATHS[field][1] for field in conservation_fields
    ]
    core_added = [*PROVENANCE_COLUMNS, *conservation_columns]
    write_csv(
        output_path,
        core_rows,
        list(rows[0]) + [
            column for column in core_added if column not in rows[0]
        ],
    )
    write_csv(optional_path, optional_rows, list(optional_rows[0]))

    status_counts = Counter(row["dbnsfp_lookup_status"] for row in core_rows)
    audit = {
        **audit_basic_rows(core_rows),
        "audit_status": (
            "PASS" if not status_counts.get("COORDINATE_MISMATCH") else "WARNING"
        ),
        "access_service": ACCESS_SERVICE,
        "api_endpoint": api_endpoint,
        "metadata_endpoint": metadata_endpoint,
        "assembly_requested": assembly,
        "coordinate_and_allele_verification": (
            "Required before accepting any dbNSFP record."
        ),
        "source_release": source_release,
        "expected_source_release": expected_release,
        "myvariant_build_version": build_version,
        "myvariant_build_date": build_date,
        "dbnsfp_download_date": source_metadata.get("download_date", ""),
        "dbnsfp_upload_date": source_metadata.get("upload_date", ""),
        "dbnsfp_license_url": source_metadata.get("license_url", ""),
        "request_fields": REQUEST_FIELDS,
        "batch_size": batch_size,
        "total_batches": len(batch_list),
        "api_batches_requested": api_batches_requested,
        "cache_batches_reused": batches_reused,
        "cache_directory": str(cache_dir.resolve()),
        "lookup_status_counts": dict(status_counts),
        "conservation_nonmissing_counts": {
            column: sum(bool(row.get(column)) for row in core_rows)
            for column in conservation_columns
        },
        "conservation_reduction": (
            "Maximum valid numeric value only when MyVariant returns multiple values for "
            "the same coordinate-verified dbNSFP record."
        ),
        "predictor_policy": (
            "Composite predictor scores are written to a separate optional file and are "
            "not included in the core evidence table."
        ),
        "predictor_reduction": (
            "Maximum valid numeric value across values returned by MyVariant without "
            "assuming a selected-transcript alignment; these columns remain optional "
            "pending transcript-alignment and circularity audits."
        ),
        "predictor_nonmissing_counts": {
            predictor_column(field): sum(
                bool(row.get(predictor_column(field))) for row in optional_rows
            )
            for field in predictor_fields
        },
        "smoke_test_rows": len(rows) if args.smoke_test is not None else 0,
    }
    write_json(audit_path, audit)
    append_stage_manifest(
        stage_name,
        [input_path, metadata_path, *cache_paths],
        [output_path, optional_path, audit_path],
        {
            "access_service": ACCESS_SERVICE,
            "api_endpoint": api_endpoint,
            "metadata_endpoint": metadata_endpoint,
            "assembly": assembly,
            "source_release": source_release,
            "myvariant_build_version": build_version,
            "request_fields": REQUEST_FIELDS,
            "batch_size": batch_size,
            "conservation_fields": conservation_fields,
            "optional_predictor_fields": predictor_fields,
            "smoke_test": args.smoke_test,
        },
        notes=[
            "The API response cache is content-addressed by query IDs, assembly, fields, "
            "dbNSFP release and MyVariant build.",
            "MyVariant may expose an older dbNSFP release than the current dbNSFP website; "
            "the API metadata recorded here is authoritative for this stage.",
        ],
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
