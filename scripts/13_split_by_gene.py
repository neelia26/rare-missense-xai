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
    groups: dict[str, dict[str, Any]],
    fractions: dict[str, float],
    seed: int,
    search_attempts: int = 2000,
) -> dict[str, str]:
    if search_attempts < 1:
        raise ValueError("search_attempts must be positive")
    split_order = list(fractions)
    genes = sorted(groups)
    total_variants = sum(group["total"] for group in groups.values())
    target_pathogenic = sum(group["1"] for group in groups.values()) / max(total_variants, 1)
    gene_targets = {
        split: int(round(len(genes) * fractions[split])) for split in split_order
    }
    gene_targets[split_order[0]] += len(genes) - sum(gene_targets.values())
    best: tuple[float, dict[str, str]] | None = None
    for attempt in range(search_attempts):
        shuffled = genes.copy()
        random.Random(seed + attempt).shuffle(shuffled)
        assignments: dict[str, str] = {}
        offset = 0
        for split in split_order:
            count = gene_targets[split]
            for gene in shuffled[offset : offset + count]:
                assignments[gene] = split
            offset += count
        score = 0.0
        for split, fraction in fractions.items():
            assigned = [groups[gene] for gene in genes if assignments[gene] == split]
            variants = sum(group["total"] for group in assigned)
            pathogenic = sum(group["1"] for group in assigned)
            observed_fraction = variants / max(total_variants, 1)
            observed_pathogenic = pathogenic / max(variants, 1)
            score += ((observed_fraction - fraction) / 0.01) ** 2
            score += ((observed_pathogenic - target_pathogenic) / 0.01) ** 2
        candidate = (score, assignments)
        if best is None or candidate[0] < best[0]:
            best = candidate
    assert best is not None
    return best[1]


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
    parser.add_argument("--no-population")
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "data/processed/splits"))
    parser.add_argument(
        "--assignment-output",
        default=str(PROJECT_ROOT / "data/processed/split_gene_assignments.csv"),
    )
    parser.add_argument(
        "--audit-output", default=str(PROJECT_ROOT / "audits/13_gene_split_audit.json")
    )
    parser.add_argument("--search-attempts", type=int, default=2000)
    parser.add_argument("--minimum-evaluation-genes", type=int, default=0)
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
    assignments = assign_groups(groups, fractions, seed, args.search_attempts)

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
    assignment_path = Path(args.assignment_output)
    write_csv(assignment_path, assignment_rows, list(assignment_rows[0]))
    output_dir = Path(args.output_dir)
    allowed_variant_keys = {row["variant_key"] for row in full_rows}

    def split_to_output(rows: list[dict[str, str]], stem: str) -> list[Path]:
        outputs_local: list[Path] = []
        fieldnames = list(rows[0])
        for split in ("train", "validation", "test"):
            selected = [
                row
                for row in rows
                if row.get("variant_key") in allowed_variant_keys
                if assignments[row["reconciled_gene_symbol"]] == split
            ]
            path = output_dir / f"{split}_{stem}.csv"
            write_csv(path, selected, fieldnames)
            outputs_local.append(path)
        return outputs_local

    outputs = [assignment_path, *split_to_output(full_rows, "full")]

    for path_value, stem in ((args.core, "core"), (args.optional, "optional_predictors")):
        path = Path(path_value)
        if path.exists():
            table_rows = read_csv(path)
            outputs.extend(split_to_output(table_rows, stem))

    if args.no_population:
        no_population_path = Path(args.no_population)
        if no_population_path.exists():
            outputs.extend(split_to_output(read_csv(no_population_path), "core_no_population"))

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
    total_rows = len(full_rows)
    overall_pathogenic = sum(row["binary_label"] == "1" for row in full_rows) / total_rows
    split_metrics = {}
    for split in fractions:
        total = sum(split_counts[split].values())
        split_metrics[split] = {
            "variant_count": total,
            "variant_fraction": round(total / total_rows, 6),
            "pathogenic_proportion": round(split_counts[split]["1"] / max(total, 1), 6),
            "pathogenic_proportion_delta_from_overall": round(
                abs(split_counts[split]["1"] / max(total, 1) - overall_pathogenic), 6
            ),
            "gene_count": len(split_genes[split]),
        }
    minimum_genes_pass = all(
        len(split_genes[split]) >= args.minimum_evaluation_genes
        for split in ("validation", "test")
    )
    audit = {
        "random_seed": seed,
        "target_fractions": fractions,
        "variant_counts": {split: dict(counts) for split, counts in split_counts.items()},
        "gene_counts": {split: len(values) for split, values in split_genes.items()},
        "gene_overlap": overlap,
        "gene_overlap_pass": not any(overlap.values()),
        "split_metrics": split_metrics,
        "minimum_evaluation_genes": args.minimum_evaluation_genes,
        "minimum_evaluation_genes_pass": minimum_genes_pass,
        "search_attempts": args.search_attempts,
        "selection_rule": (
            "Deterministic multi-start allocation fixes the requested gene proportions, then "
            "selects the allocation that best matches variant and class targets while keeping "
            "every reconciled gene in exactly one partition."
        ),
        "usage_rule": (
            "Fit preprocessing and models on train genes, choose/tune on validation genes, and "
            "use test genes once for the final locked evaluation."
        ),
    }
    if not minimum_genes_pass:
        raise RuntimeError("Evaluation partitions do not meet the minimum gene-count requirement")
    audit_path = Path(args.audit_output)
    write_json(audit_path, audit)
    outputs.append(audit_path)
    append_stage_manifest(
        "13_split_by_gene",
        [full_path, Path(args.core)],
        outputs,
        {
            "fractions": fractions,
            "seed": seed,
            "search_attempts": args.search_attempts,
            "minimum_evaluation_genes": args.minimum_evaluation_genes,
        },
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
