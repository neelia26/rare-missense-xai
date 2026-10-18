#!/usr/bin/env python3
"""Run pinned local VEP, use checkpointed REST VEP, or register existing VEP JSONL."""

from __future__ import annotations

import argparse
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

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    copy_atomic,
    load_config,
    read_csv,
    write_json,
)


REST_BASE = "https://rest.ensembl.org/vep/homo_sapiens/region"
REST_OPTIONS = {
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


def post_batch(
    batch_number: int, variants: list[str], max_retries: int
) -> tuple[int, list[dict], int]:
    url = f"{REST_BASE}?{urllib.parse.urlencode(REST_OPTIONS)}"
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
                "User-Agent": "rare-missense-rebuild/1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=300) as response:
                return batch_number, json.load(response), attempt
        except urllib.error.HTTPError as error:
            last_error = f"HTTP {error.code}: {error.reason}"
            if error.code not in {429, 500, 502, 503, 504}:
                raise RuntimeError(last_error) from error
            retry_after = error.headers.get("Retry-After")
            delay = float(retry_after) if retry_after else min(60.0, 2**attempt + random.random())
            time.sleep(delay)
        except (urllib.error.URLError, TimeoutError) as error:
            last_error = str(error)
            time.sleep(min(60.0, 2**attempt + random.random()))
    raise RuntimeError(f"VEP REST batch {batch_number} failed: {last_error}")


def run_rest(
    candidates: Path,
    output: Path,
    batch_size: int,
    workers: int,
    retries: int,
    limit: int | None,
) -> dict:
    rows = read_csv(candidates)
    if limit is not None:
        rows = rows[:limit]
    batches = [rows[offset : offset + batch_size] for offset in range(0, len(rows), batch_size)]
    checkpoints = PROJECT_ROOT / "data/raw/vep/rest_checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    pending = [
        (number, rows_in_batch)
        for number, rows_in_batch in enumerate(batches)
        if not (checkpoints / f"batch_{number:05d}.json").exists()
    ]
    retry_counts: dict[str, int] = {}
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                post_batch, number, [rest_variant(row) for row in rows_in_batch], retries
            ): number
            for number, rows_in_batch in pending
        }
        for complete, future in enumerate(as_completed(futures), start=1):
            number, results, retries_used = future.result()
            write_json(checkpoints / f"batch_{number:05d}.json", results)
            retry_counts[str(retries_used)] = retry_counts.get(str(retries_used), 0) + 1
            print(f"REST batch {complete}/{len(pending)} completed", flush=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    results_written = 0
    with output.open("w", encoding="utf-8") as out:
        for number in range(len(batches)):
            checkpoint = checkpoints / f"batch_{number:05d}.json"
            results = json.loads(checkpoint.read_text(encoding="utf-8"))
            for result in results:
                out.write(json.dumps(result, separators=(",", ":")) + "\n")
                results_written += 1
    return {
        "mode": "rest",
        "candidate_rows": len(rows),
        "result_objects": results_written,
        "batch_size": batch_size,
        "batch_count": len(batches),
        "checkpoint_batches_reused": len(batches) - len(pending),
        "retry_counts": retry_counts,
        "endpoint": REST_BASE,
        "options": REST_OPTIONS,
    }


def local_command(config: dict, input_vcf: Path, output: Path) -> list[str]:
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
    parser.add_argument("--limit", type=int, help="REST pilot only")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    mode = args.mode or config["vep"]["mode"]
    input_vcf = Path(args.input_vcf)
    candidates = Path(args.candidates)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    if mode == "register":
        if not args.existing_jsonl:
            raise SystemExit("--existing-jsonl is required in register mode")
        source = Path(args.existing_jsonl).expanduser().resolve()
        copy_atomic(source, output)
        report = {"mode": "register", "registered_source": str(source)}
        inputs = [source]
    elif mode == "rest":
        report = run_rest(
            candidates,
            output,
            int(config["vep"]["rest_batch_size"]),
            int(config["vep"]["rest_workers"]),
            int(config["vep"]["rest_max_retries"]),
            args.limit,
        )
        inputs = [candidates]
    else:
        command = local_command(config, input_vcf, output)
        if args.dry_run:
            print(json.dumps({"mode": "local", "command": command}, indent=2))
            return
        if not shutil.which(command[0]):
            raise RuntimeError(
                f"Cannot find VEP executable {command[0]!r}. Install VEP/cache or use "
                "--mode rest for a small pilot."
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
    append_stage_manifest("05_run_vep", inputs, [output, report_path], report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

