#!/usr/bin/env python3
"""Create deterministic train/validation/test partitions with no gene overlap."""

from __future__ import annotations

import argparse
import json
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.pipeline_utils import (
    PROJECT_ROOT,
    append_stage_manifest,
    load_config,
    read_csv,
    write_csv,
    write_json,
)


def group_statistics(rows: list[dict[str, str]]) -> dict[str, dict[str, Any]]:
    groups: dict[str, dict[str, Any]] = {}
    for row in rows:
        gene = row.get("reconciled_gene_symbol", "").strip()
        if not gene:
            raise ValueError(f"Variant {row.get('variant_key')} has no grouping gene")
        stats = groups.setdefault(gene, {"rows": [], "0": 0, "1": 0, "total": 0})
        stats["rows"].append(row)
        stats[row["binary_label"]] += 1
        stats["total"] += 1
    return groups


def assign_groups(
    groups: dict[str, dict[str, Any]], fractions: dict[str, float], seed: int
) -> dict[str, str]:
    totals = {"0": sum(group["0"] for group in groups.values()), "1": sum(group["1"] for group in groups.values())}
    targets = {
        split: {
            "0": totals["0"] * fraction,
            "1": totals["1"] * fraction,
            "total": (totals["0"] + totals["1"]) * fraction,
        }
        for split, fraction in fractions.items()
    }
    current = {split: {"0": 0, "1": 0, "total": 0} for split in fractions}
    rng = random.Random(seed)
    tie_breakers = {gene: rng.random() for gene in groups}
    ordered = sorted(
        groups,
        key=lambda gene: (
            -groups[gene]["total"],
            -abs(groups[gene]["0"] - groups[gene]["1"]),
            tie_breakers[gene],
        ),
    )
    assignments: dict[str, str] = {}
    for gene in ordered:
        group = groups[gene]
        candidates: list[tuple[float, float, str]] = []
        for split in fractions:
            score = 0.0
            overfill = 0.0
            for label in ("0", "1", "total"):
                projected = current[split][label] + group[label]
                target = max(targets[split][label], 1.0)
                score += ((projected - target) / target) ** 2
                overfill += max(0.0, projected - target) / target
            deficit = targets[split]["total"] - current[split]["total"]
            candidates.append((score + 10.0 * overfill, -deficit, split))
        _, _, selected = min(candidates)
        assignments[gene] = selected
        for label in ("0", "1", "total"):
            current[selected][label] += group[label]
    return assignments


def split_table(
    rows: list[dict[str, str]], assignments: dict[str, str], output_stem: str
) -> list[Path]:
    outputs: list[Path] = []
    fieldnames = list(rows[0])
    for split in ("train", "validation", "test"):
        selected = [
            row for row in rows if assignments[row["reconciled_gene_symbol"]] == split
        ]
        path = PROJECT_ROOT / f"data/processed/splits/{split}_{output_stem}.csv"
        write_csv(path, selected, fieldnames)
        outputs.append(path)
    return outputs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", default=str(PROJECT_ROOT / "data/processed/final_dataset_full.csv"))
    parser.add_argument("--core", default=str(PROJECT_ROOT / "data/processed/model_table_core.csv"))
    parser.add_argument(
        "--optional", default=str(PROJECT_ROOT / "data/processed/model_table_with_optional_predictors.csv")
    )
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/pipeline_config.json"))
    args = parser.parse_args()
    config = load_config(args.config)
    fractions = {
        "train": float(config["splits"]["train_fraction"]),
        "validation": float(config["splits"]["validation_fraction"]),
        "test": float(config["splits"]["test_fraction"]),
    }
    if abs(sum(fractions.values()) - 1.0) > 1e-9:
        raise ValueError("Split fractions must sum to 1")
    seed = int(config["project"]["random_seed"])
    full_path = Path(args.full)
    full_rows = read_csv(full_path)
    groups = group_statistics(full_rows)
    assignments = assign_groups(groups, fractions, seed)

    assignment_rows = [
        {
            "gene": gene,
            "split": split,
            "variant_count": groups[gene]["total"],
            "benign_count": groups[gene]["0"],
            "pathogenic_count": groups[gene]["1"],
        }
        for gene, split in sorted(assignments.items())
    ]
    assignment_path = PROJECT_ROOT / "data/processed/split_gene_assignments.csv"
    write_csv(assignment_path, assignment_rows, list(assignment_rows[0]))
    outputs = [assignment_path, *split_table(full_rows, assignments, "full")]

    for path_value, stem in ((args.core, "core"), (args.optional, "optional_predictors")):
        path = Path(path_value)
        if path.exists():
            table_rows = read_csv(path)
            outputs.extend(split_table(table_rows, assignments, stem))

    split_counts: dict[str, Counter[str]] = defaultdict(Counter)
    split_genes: dict[str, set[str]] = defaultdict(set)
    for row in full_rows:
        gene = row["reconciled_gene_symbol"]
        split = assignments[gene]
        split_counts[split][row["binary_label"]] += 1
        split_genes[split].add(gene)
    overlap = {
        "train_validation": sorted(split_genes["train"] & split_genes["validation"]),
        "train_test": sorted(split_genes["train"] & split_genes["test"]),
        "validation_test": sorted(split_genes["validation"] & split_genes["test"]),
    }
    audit = {
        "random_seed": seed,
        "target_fractions": fractions,
        "variant_counts": {split: dict(counts) for split, counts in split_counts.items()},
        "gene_counts": {split: len(values) for split, values in split_genes.items()},
        "gene_overlap": overlap,
        "gene_overlap_pass": not any(overlap.values()),
        "selection_rule": (
            "The deterministic greedy assignment balances class and total counts while keeping "
            "every reconciled gene in exactly one partition."
        ),
        "usage_rule": (
            "Fit preprocessing and models on train genes, choose/tune on validation genes, and "
            "use test genes once for the final locked evaluation."
        ),
    }
    audit_path = PROJECT_ROOT / "audits/13_gene_split_audit.json"
    write_json(audit_path, audit)
    outputs.append(audit_path)
    append_stage_manifest("13_split_by_gene", [full_path, Path(args.core)], outputs, {"fractions": fractions, "seed": seed})
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()

