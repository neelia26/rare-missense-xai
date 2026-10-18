#!/usr/bin/env python3
"""Download or register one ClinVar GRCh38 VCF and freeze its provenance."""

from __future__ import annotations

import argparse
import gzip
import json
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    copy_atomic,
    ensure_layout,
    file_record,
    load_config,
    utc_now,
    write_json,
)


def read_vcf_header(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            for key in ("fileDate", "source", "reference"):
                prefix = f"##{key}="
                if line.startswith(prefix):
                    values[key] = line.strip().split("=", 1)[1]
    return values


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "rare-missense-rebuild/1.0"})
    with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as out:
        while chunk := response.read(1024 * 1024):
            out.write(chunk)
    temporary.replace(destination)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Use --local-vcf for an already downloaded release or --url to download one."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--local-vcf", help="Existing ClinVar GRCh38 .vcf.gz file")
    group.add_argument("--url", help="Explicit ClinVar VCF URL")
    parser.add_argument("--expected-sha256", default="")
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    ensure_layout()
    config = load_config(args.config)
    destination = PROJECT_ROOT / "data/raw/clinvar_grch38.vcf.gz"
    source_url = args.url or config["clinvar"]["source_url"]

    if destination.exists() and not args.force:
        raise SystemExit(
            f"{destination} already exists. This protects the frozen release. "
            "Use --force only if you intentionally want to replace it."
        )
    if args.local_vcf:
        source = Path(args.local_vcf).expanduser().resolve()
        if not source.exists():
            raise FileNotFoundError(source)
        copy_atomic(source, destination)
        acquisition = {"method": "local_file", "original_path": str(source)}
    else:
        download(source_url, destination)
        acquisition = {"method": "download", "url": source_url}

    record = file_record(destination)
    if args.expected_sha256 and record["sha256"].lower() != args.expected_sha256.lower():
        destination.unlink(missing_ok=True)
        raise RuntimeError("Downloaded/registered VCF did not match --expected-sha256")

    header = read_vcf_header(destination)
    provenance = {
        "registered_at_utc": utc_now(),
        "source_name": "NCBI ClinVar GRCh38 VCF",
        "acquisition": acquisition,
        "vcf_header": header,
        "file": record,
        "frozen": True,
        "warning": (
            "The default URL points to the live monthly file. Preserve this local file and its "
            "SHA-256; use an archived monthly URL when one is available."
        ),
    }
    provenance_path = PROJECT_ROOT / "data/raw/clinvar_grch38.provenance.json"
    write_json(provenance_path, provenance)
    append_stage_manifest(
        "01_fetch_clinvar",
        [],
        [destination, provenance_path],
        {"source_url": source_url, "expected_sha256_supplied": bool(args.expected_sha256)},
    )
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()

