#!/usr/bin/env python3
"""Check the local environment without changing or downloading research data."""

from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import PROJECT_ROOT, ensure_layout, load_config, utc_now, write_json


def command_version(command: str, arguments: list[str]) -> dict[str, str | bool]:
    executable = shutil.which(command)
    if not executable:
        return {"available": False, "path": "", "version": ""}
    try:
        result = subprocess.run(
            [executable, *arguments],
            check=False,
            capture_output=True,
            text=True,
            timeout=30,
        )
        output = (result.stdout or result.stderr).strip().splitlines()
        version = output[0] if output else "available"
    except Exception as error:  # environment diagnostics must continue
        version = f"could not read version: {error}"
    return {"available": True, "path": executable, "version": version}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    args = parser.parse_args()
    ensure_layout()
    config = load_config(args.config)

    packages = {}
    for package in ("pandas", "numpy", "sklearn", "pysam"):
        spec = importlib.util.find_spec(package)
        packages[package] = {"available": spec is not None}
        if spec is not None:
            module = __import__(package)
            packages[package]["version"] = getattr(module, "__version__", "unknown")

    report = {
        "checked_at_utc": utc_now(),
        "project_root": str(PROJECT_ROOT),
        "platform": platform.platform(),
        "python": sys.version,
        "packages": packages,
        "commands": {
            "vep": command_version(config["vep"]["executable"], ["--help"]),
            "docker": command_version("docker", ["--version"]),
            "git": command_version("git", ["--version"]),
        },
        "python_dependencies_ready": all(packages[name]["available"] for name in packages),
        "local_vep_ready": command_version(config["vep"]["executable"], ["--help"])["available"],
        "notes": [
            "The configured full workflow uses checkpointed Ensembl VEP REST requests.",
            "Stage 05 verifies the live Ensembl release before accepting annotations.",
            "pysam is required for optional direct lookups in indexed gnomAD VCFs.",
        ],
    }
    destination = PROJECT_ROOT / "audits/00_environment_report.json"
    write_json(destination, report)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

