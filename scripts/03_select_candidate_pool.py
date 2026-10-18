#!/usr/bin/env python3
"""Select a reproducible, label-stratified pool before costly enrichment."""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    audit_basic_rows,
    load_config,
    read_csv,
    write_csv,
    write_json,
)


def sample(rows: list[dict[str, str]], target: int, seed: int) -> list[dict[str, str]]:
    if len(rows) < target:
        raise RuntimeError(f"Requested {target:,} rows, but only {len(rows):,} are available")
    rng = random.Random(seed)
    selected = rng.sample(rows, target)
    return sorted(selected, key=lambda row: row["variant_key"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=str(PROJECT_ROOT / "data/interim/02_clinvar_candidates.csv"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "data/interim/03_candidate_pool.csv"))
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    parser.add_argument("--benign", type=int, help="Override candidate benign count")
    parser.add_argument("--pathogenic", type=int, help="Override candidate pathogenic count")
    args = parser.parse_args()
    config = load_config(args.config)
    seed = int(config["project"]["random_seed"])
    benign_target = args.benign or int(config["sampling"]["candidate_benign"])
    pathogenic_target = args.pathogenic or int(config["sampling"]["candidate_pathogenic"])
    input_path = Path(args.input)
    output_path = Path(args.output)
    rows = read_csv(input_path)
    by_label = {
        0: [row for row in rows if row["binary_label"] == "0"],
        1: [row for row in rows if row["binary_label"] == "1"],
    }
    selected = sample(by_label[0], benign_target, seed) + sample(
        by_label[1], pathogenic_target, seed + 1
    )
    random.Random(seed + 2).shuffle(selected)
    fieldnames = list(rows[0])
    write_csv(output_path, selected, fieldnames)
    audit = {
        **audit_basic_rows(selected),
        "eligible_class_counts": dict(Counter(row["binary_label"] for row in rows)),
        "candidate_targets": {"benign": benign_target, "pathogenic": pathogenic_target},
        "random_seed": seed,
        "purpose": (
            "Oversample before rarity, transcript and protein validation so that the final "
            "balanced cohort can be selected only after biological checks."
        ),
    }
    audit_path = PROJECT_ROOT / "audits/03_candidate_pool_audit.json"
    write_json(audit_path, audit)
    append_stage_manifest(
        "03_select_candidate_pool",
        [input_path],
        [output_path, audit_path],
        {"benign_target": benign_target, "pathogenic_target": pathogenic_target, "seed": seed},
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
