#!/usr/bin/env python3
"""Run pinned local VEP, checkpointed REST VEP, or register existing VEP JSONL."""

from __future__ import annotations

import argparse
import hashlib
import json
import random
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    copy_atomic,
    load_config,
    read_csv,
    utc_now,
    write_json,
)


DEFAULT_REST_ENDPOINT = "https://rest.ensembl.org/vep/homo_sapiens/region"
DEFAULT_METADATA_ENDPOINT = "https://rest.ensembl.org/info/software"
DEFAULT_REST_OPTIONS = {
    "mane": 1,
    "canonical": 1,
    "protein": 1,
    "uniprot": 1,
    "hgvs": 1,
    "domains": 1,
    "numbers": 1,
    "variant_class": 1,
    "af": 1,
    "af_gnomade": 1,
    "af_gnomadg": 1,
    "max_af": 1,
}


def rest_variant(row: dict[str, str]) -> str:
    return (
        f"{row['chromosome']} {row['position']} {row['variant_key']} "
        f"{row['reference_allele']} {row['alternate_allele']} . . ."
    )


def canonical_digest(value: Any) -> str:
    rendered = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(rendered.encode("utf-8")).hexdigest()


def request_fingerprint(
    endpoint: str,
    options: dict[str, Any],
    assembly: str,
    api_release: str,
    variants: list[str],
) -> str:
    """Identify one VEP request by its complete scientific and technical inputs."""
    return canonical_digest(
        {
            "endpoint": endpoint,
            "options": options,
            "assembly": assembly,
            "api_release": api_release,
            "variants": variants,
        }
    )


def checkpoint_path(directory: Path, batch_number: int, fingerprint: str) -> Path:
    return directory / f"batch_{batch_number:05d}_{fingerprint[:16]}.json"


def retry_delay(error: urllib.error.HTTPError, attempt: int) -> float:
    retry_after = error.headers.get("Retry-After")
    if retry_after:
        try:
            return min(300.0, float(retry_after))
        except ValueError:
            pass
    return min(60.0, 2**attempt + random.random())


def get_json(url: str, max_retries: int, user_agent: str) -> dict[str, Any]:
    last_error = "unknown"
    for attempt in range(max_retries):
        request = urllib.request.Request(
            url,
            headers={"Accept": "application/json", "User-Agent": user_agent},
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                result = json.load(response)
            if not isinstance(result, dict):
                raise RuntimeError(
                    f"VEP metadata endpoint returned {type(result).__name__}; expected object"
                )
            return result
        except urllib.error.HTTPError as error:
            last_error = f"HTTP {error.code}: {error.reason}"
            if error.code not in {429, 500, 502, 503, 504}:
                break
            delay = retry_delay(error, attempt)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as error:
            last_error = f"{type(error).__name__}: {error}"
            delay = min(60.0, 2**attempt + random.random())
        if attempt + 1 < max_retries:
            time.sleep(delay)
    raise RuntimeError(f"Could not read Ensembl VEP metadata from {url}: {last_error}")


def validate_api_release(metadata: dict[str, Any], expected_release: str) -> str:
    actual = str(metadata.get("release", "")).strip()
    if not actual:
        raise RuntimeError("Ensembl metadata did not report a release")
    if expected_release and actual != str(expected_release):
        raise RuntimeError(
            f"Ensembl REST currently reports release {actual}, but release "
            f"{expected_release} is required. Review the source before changing the config."
        )
    return actual


def validate_response(results: Any, variants: list[str], context: str) -> list[dict[str, Any]]:
    if not isinstance(results, list):
        raise RuntimeError(f"{context}: VEP returned {type(results).__name__}; expected list")
    if len(results) != len(variants):
        raise RuntimeError(
            f"{context}: VEP returned {len(results)} objects for {len(variants)} variants"
        )
    if any(not isinstance(item, dict) for item in results):
        raise RuntimeError(f"{context}: VEP response contains a non-object item")
    return results


def load_checkpoint(
    path: Path,
    fingerprint: str,
    variants: list[str],
    api_release: str,
) -> list[dict[str, Any]] | None:
    if not path.exists():
        return None
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None
    if not isinstance(record, dict):
        return None
    if (
        record.get("schema_version") != 1
        or record.get("request_sha256") != fingerprint
        or str(record.get("api_release", "")) != api_release
        or record.get("variant_count") != len(variants)
    ):
        return None
    try:
        return validate_response(record.get("response"), variants, f"cached {path.name}")
    except RuntimeError:
        return None


def post_batch(
    batch_number: int,
    variants: list[str],
    endpoint: str,
    options: dict[str, Any],
    max_retries: int,
    user_agent: str,
) -> tuple[int, list[dict[str, Any]], int]:
    separator = "&" if "?" in endpoint else "?"
    url = f"{endpoint}{separator}{urllib.parse.urlencode(options)}"
    body = json.dumps({"variants": variants}).encode("utf-8")
    last_error = "unknown"
    for attempt in range(max_retries):
        request = urllib.request.Request(
            url,
            data=body,
            method="POST",
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json",
                "User-Agent": user_agent,
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                results = json.load(response)
            validated = validate_response(results, variants, f"REST batch {batch_number}")
            return batch_number, validated, attempt
        except urllib.error.HTTPError as error:
            last_error = f"HTTP {error.code}: {error.reason}"
            if error.code not in {429, 500, 502, 503, 504}:
                raise RuntimeError(last_error) from error
            delay = retry_delay(error, attempt)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError, OSError) as error:
            last_error = f"{type(error).__name__}: {error}"
            delay = min(60.0, 2**attempt + random.random())
        if attempt + 1 < max_retries:
            time.sleep(delay)
    raise RuntimeError(f"VEP REST batch {batch_number} failed: {last_error}")


def run_rest(
    candidates: Path,
    output: Path,
    settings: dict[str, Any],
    limit: int | None,
) -> tuple[dict[str, Any], Path]:
    endpoint = str(settings.get("rest_endpoint", DEFAULT_REST_ENDPOINT))
    metadata_endpoint = str(
        settings.get("rest_metadata_endpoint", DEFAULT_METADATA_ENDPOINT)
    )
    user_agent = str(
        settings.get("user_agent", "rare-missense-evidence-pipeline/1.2")
    )
    retries = int(settings.get("rest_max_retries", 6))
    batch_size = int(settings.get("rest_batch_size", 200))
    workers = int(settings.get("rest_workers", 2))
    assembly = str(settings.get("assembly", "GRCh38"))
    expected_release = str(settings.get("release", ""))
    options = dict(settings.get("rest_options", DEFAULT_REST_OPTIONS))
    if not 1 <= batch_size <= 200:
        raise ValueError("VEP REST batch size must be between 1 and 200")
    if workers < 1 or retries < 1:
        raise ValueError("VEP REST workers and retries must both be at least 1")

    metadata_started_at = utc_now()
    metadata = get_json(metadata_endpoint, retries, user_agent)
    api_release = validate_api_release(metadata, expected_release)
    checkpoint_setting = settings.get(
        "rest_checkpoint_dir", "data/raw/vep/rest_checkpoints"
    )
    checkpoints = Path(str(checkpoint_setting)).expanduser()
    if not checkpoints.is_absolute():
        checkpoints = PROJECT_ROOT / checkpoints
    checkpoints.mkdir(parents=True, exist_ok=True)
    metadata_path = checkpoints / f"api_metadata_release_{api_release}.json"
    write_json(
        metadata_path,
        {
            "retrieved_at_utc": metadata_started_at,
            "endpoint": metadata_endpoint,
            "expected_release": expected_release,
            "actual_release": api_release,
            "start_response": metadata,
        },
    )

    rows = read_csv(candidates)
    if limit is not None:
        if limit < 1:
            raise ValueError("--limit must be at least 1")
        rows = rows[:limit]
    if not rows:
        raise ValueError(f"No candidate rows were found in {candidates}")
    batches = [rows[offset : offset + batch_size] for offset in range(0, len(rows), batch_size)]
    batch_records: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    reused = 0
    for number, batch_rows in enumerate(batches):
        variants = [rest_variant(row) for row in batch_rows]
        fingerprint = request_fingerprint(
            endpoint, options, assembly, api_release, variants
        )
        path = checkpoint_path(checkpoints, number, fingerprint)
        cached = load_checkpoint(path, fingerprint, variants, api_release)
        record = {
            "number": number,
            "variants": variants,
            "fingerprint": fingerprint,
            "path": path,
            "response": cached,
        }
        batch_records.append(record)
        if cached is None:
            pending.append(record)
        else:
            reused += 1

    retry_counts: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                post_batch,
                record["number"],
                record["variants"],
                endpoint,
                options,
                retries,
                user_agent,
            ): record
            for record in pending
        }
        for complete, future in enumerate(as_completed(futures), start=1):
            record = futures[future]
            number, results, retries_used = future.result()
            cache_record = {
                "schema_version": 1,
                "retrieved_at_utc": utc_now(),
                "request_sha256": record["fingerprint"],
                "batch_number": number,
                "endpoint": endpoint,
                "options": options,
                "assembly": assembly,
                "api_release": api_release,
                "variant_count": len(record["variants"]),
                "variants_sha256": canonical_digest(record["variants"]),
                "response_count": len(results),
                "response": results,
            }
            write_json(record["path"], cache_record)
            record["response"] = results
            retry_counts[str(retries_used)] = retry_counts.get(str(retries_used), 0) + 1
            print(f"REST batch {complete}/{len(pending)} completed", flush=True)

    ending_metadata = get_json(metadata_endpoint, retries, user_agent)
    ending_release = validate_api_release(ending_metadata, expected_release)
    if ending_release != api_release:
        raise RuntimeError(
            f"Ensembl release changed from {api_release} to {ending_release} during Stage 05"
        )
    write_json(
        metadata_path,
        {
            "started_at_utc": metadata_started_at,
            "completed_at_utc": utc_now(),
            "endpoint": metadata_endpoint,
            "expected_release": expected_release,
            "actual_release": api_release,
            "start_response": metadata,
            "end_response": ending_metadata,
        },
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    results_written = 0
    with output.open("w", encoding="utf-8") as out:
        for record in batch_records:
            results = record["response"]
            if results is None:
                raise RuntimeError(f"VEP batch {record['number']} has no validated response")
            for result in results:
                out.write(json.dumps(result, separators=(",", ":")) + "\n")
                results_written += 1
    if results_written != len(rows):
        raise RuntimeError(
            f"VEP output contains {results_written} results for {len(rows)} candidates"
        )
    report = {
        "mode": "rest",
        "candidate_rows": len(rows),
        "result_objects": results_written,
        "batch_size": batch_size,
        "batch_count": len(batches),
        "checkpoint_batches_reused": reused,
        "api_batches_requested": len(pending),
        "retry_counts": retry_counts,
        "endpoint": endpoint,
        "metadata_endpoint": metadata_endpoint,
        "expected_release": expected_release,
        "actual_release": api_release,
        "release_verified_at_start_and_end": True,
        "assembly": assembly,
        "options": options,
        "checkpoint_directory": str(checkpoints.resolve()),
        "checkpoint_policy": (
            "Each checkpoint name is content-addressed by variants, endpoint, options, "
            "assembly and verified Ensembl release."
        ),
        "limited_run_rows": limit or 0,
    }
    return report, metadata_path


def local_command(config: dict[str, Any], input_vcf: Path, output: Path) -> list[str]:
    settings = config["vep"]
    command = [
        settings["executable"],
        "--input_file",
        str(input_vcf),
        "--format",
        "vcf",
        "--output_file",
        str(output),
        "--json",
        "--cache",
        "--offline",
        "--assembly",
        settings["assembly"],
        "--cache_version",
        str(settings["release"]),
        "--symbol",
        "--canonical",
        "--mane",
        "--protein",
        "--uniprot",
        "--domains",
        "--numbers",
        "--hgvs",
        "--variant_class",
        "--af",
        "--af_gnomade",
        "--af_gnomadg",
        "--max_af",
        "--fork",
        str(settings["forks"]),
        "--no_stats",
        "--force_overwrite",
    ]
    if settings.get("cache_dir"):
        command.extend(["--dir_cache", str(Path(settings["cache_dir"]).expanduser())])
    return command


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["local", "rest", "register"])
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    parser.add_argument("--candidates", default=str(PROJECT_ROOT / "data/interim/03_candidate_pool.csv"))
    parser.add_argument("--input-vcf", default=str(PROJECT_ROOT / "data/raw/vep/04_candidates.vcf"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/raw/vep/05_vep.jsonl"))
    parser.add_argument("--existing-jsonl")
    parser.add_argument("--limit", type=int, help="Annotate only the first N candidate rows")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    mode = args.mode or config["vep"]["mode"]
    input_vcf = Path(args.input_vcf)
    candidates = Path(args.candidates)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    extra_outputs: list[Path] = []

    if mode == "register":
        if not args.existing_jsonl:
            raise SystemExit("--existing-jsonl is required in register mode")
        source = Path(args.existing_jsonl).expanduser().resolve()
        copy_atomic(source, output)
        report = {"mode": "register", "registered_source": str(source)}
        inputs = [source]
    elif mode == "rest":
        if args.dry_run:
            settings = config["vep"]
            print(
                json.dumps(
                    {
                        "mode": "rest",
                        "endpoint": settings.get("rest_endpoint", DEFAULT_REST_ENDPOINT),
                        "metadata_endpoint": settings.get(
                            "rest_metadata_endpoint", DEFAULT_METADATA_ENDPOINT
                        ),
                        "expected_release": settings.get("release", ""),
                        "batch_size": settings.get("rest_batch_size", 200),
                        "workers": settings.get("rest_workers", 2),
                        "candidates": str(candidates),
                    },
                    indent=2,
                )
            )
            return
        report, metadata_path = run_rest(candidates, output, config["vep"], args.limit)
        extra_outputs.append(metadata_path)
        inputs = [candidates]
    else:
        command = local_command(config, input_vcf, output)
        if args.dry_run:
            print(json.dumps({"mode": "local", "command": command}, indent=2))
            return
        if not shutil.which(command[0]):
            raise RuntimeError(
                f"Cannot find VEP executable {command[0]!r}. Install VEP/cache or use "
                "--mode rest."
            )
        subprocess.run(command, check=True)
        report = {
            "mode": "local",
            "vep_release": config["vep"]["release"],
            "assembly": config["vep"]["assembly"],
            "command": command,
        }
        inputs = [input_vcf]

    report_path = PROJECT_ROOT / "audits/05_vep_execution_report.json"
    write_json(report_path, report)
    append_stage_manifest(
        "05_run_vep", inputs, [output, report_path, *extra_outputs], report
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
