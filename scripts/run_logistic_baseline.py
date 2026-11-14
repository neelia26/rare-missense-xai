#!/usr/bin/env python3
"""Run the governed training-only logistic-regression baseline."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--data-root",
        type=Path,
        default=None,
        help="Optional root containing the frozen data paths, useful for Google Drive in Colab.",
    )
    parser.add_argument(
        "--output-dir",
        default="outputs/grouped_additive_model/01_logistic_baseline",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    repo_root = args.repo_root.resolve()
    sys.path.insert(0, str(repo_root))
    from src.modeling.logistic_baseline import run_logistic_baseline

    summary = run_logistic_baseline(
        repo_root,
        data_root=args.data_root,
        output_dir=args.output_dir,
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
