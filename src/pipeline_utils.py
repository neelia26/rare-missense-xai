#!/usr/bin/env python3
"""Small, dependency-light helpers shared by the numbered pipeline scripts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import inspect
import json
import os
import platform
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "config" / "pipeline_config.json"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    config_path = Path(path) if path else DEFAULT_CONFIG
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    config["_config_path"] = str(config_path.resolve())
    return config


def project_path(value: str | Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else PROJECT_ROOT / path


def ensure_layout() -> None:
    for relative in (
        "data/raw",
        "data/interim",
        "data/processed",
        "data/processed/splits",
        "audits",
        "logs",
    ):
        (PROJECT_ROOT / relative).mkdir(parents=True, exist_ok=True)


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    return {
        "path": str(resolved),
        "size_bytes": resolved.stat().st_size,
        "sha256": sha256_file(resolved),
    }


def write_json(path: str | Path, value: Any) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=destination.parent, delete=False
    ) as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        temporary = Path(handle.name)
    os.replace(temporary, destination)


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def iter_csv(path: str | Path) -> Iterator[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as handle:
        yield from csv.DictReader(handle)


def write_csv(
    path: str | Path,
    rows: Iterable[dict[str, Any]],
    fieldnames: list[str],
) -> int:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with tempfile.NamedTemporaryFile(
        "w", newline="", encoding="utf-8", dir=destination.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: empty_if_none(row.get(key)) for key in fieldnames})
            count += 1
        temporary = Path(handle.name)
    os.replace(temporary, destination)
    return count


def copy_atomic(source: str | Path, destination: str | Path) -> None:
    source_path = Path(source)
    destination_path = Path(destination)
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("wb", dir=destination_path.parent, delete=False) as out:
        temporary = Path(out.name)
        with source_path.open("rb") as incoming:
            shutil.copyfileobj(incoming, out)
    os.replace(temporary, destination_path)


def open_text_auto(path: str | Path):
    source = Path(path)
    if source.suffix == ".gz":
        return gzip.open(source, "rt", encoding="utf-8")
    return source.open("r", encoding="utf-8")


def empty_if_none(value: Any) -> Any:
    return "" if value is None else value


def is_missing(value: Any) -> bool:
    if value is None:
        return True
    return str(value).strip().lower() in {"", ".", "na", "nan", "none", "null"}


def parse_float(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        return float(str(value).strip())
    except ValueError:
        return None


def parse_int(value: Any) -> int | None:
    if is_missing(value):
        return None
    try:
        return int(float(str(value).strip()))
    except ValueError:
        return None


def normalize_chromosome(value: Any) -> str:
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    if text == "M":
        return "MT"
    return text


def variant_key(chromosome: Any, position: Any, reference: Any, alternate: Any) -> str:
    return (
        f"GRCh38:{normalize_chromosome(chromosome)}:{int(position)}:"
        f"{str(reference).upper()}:{str(alternate).upper()}"
    )


def split_multi_value(value: Any) -> list[str]:
    if is_missing(value):
        return []
    text = str(value).replace(",", ";")
    return [item.strip() for item in text.split(";") if item.strip()]


def unique_join(values: Iterable[Any], separator: str = ";") -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for raw in values:
        for value in split_multi_value(raw):
            if value not in seen:
                seen.add(value)
                ordered.append(value)
    return separator.join(ordered)


def append_stage_manifest(
    stage: str,
    inputs: list[str | Path],
    outputs: list[str | Path],
    parameters: dict[str, Any],
    notes: list[str] | None = None,
) -> Path:
    ensure_layout()
    caller = Path(inspect.stack()[1].filename).resolve()
    version_path = PROJECT_ROOT / "VERSION"
    record = {
        "stage": stage,
        "completed_at_utc": utc_now(),
        "pipeline_version": version_path.read_text(encoding="utf-8").strip() if version_path.exists() else "unknown",
        "python_version": sys.version.split()[0],
        "platform": platform.platform(),
        "executing_script": file_record(caller) if caller.exists() else {"path": str(caller)},
        "inputs": [file_record(path) for path in inputs if Path(path).exists()],
        "outputs": [file_record(path) for path in outputs if Path(path).exists()],
        "parameters": parameters,
        "notes": notes or [],
    }
    destination = PROJECT_ROOT / "audits" / f"{stage}_manifest.json"
    write_json(destination, record)
    return destination


def require_columns(rows: list[dict[str, Any]], columns: Iterable[str], context: str) -> None:
    if not rows:
        raise ValueError(f"{context}: input has no data rows")
    missing = [column for column in columns if column not in rows[0]]
    if missing:
        raise ValueError(f"{context}: missing required columns: {', '.join(missing)}")


def audit_basic_rows(rows: list[dict[str, Any]], key: str = "variant_key") -> dict[str, Any]:
    keys = [str(row.get(key, "")) for row in rows]
    labels: dict[str, int] = {}
    for row in rows:
        label = str(row.get("binary_label", ""))
        labels[label] = labels.get(label, 0) + 1
    return {
        "rows": len(rows),
        "columns": len(rows[0]) if rows else 0,
        "unique_keys": len(set(keys)),
        "blank_keys": sum(not key_value for key_value in keys),
        "duplicate_key_rows": len(keys) - len(set(keys)),
        "class_counts": labels,
    }
