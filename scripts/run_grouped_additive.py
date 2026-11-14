#!/usr/bin/env python3
"""Run the frozen training-only grouped additive experiment."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.modeling.grouped_additive import run_grouped_additive


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument(
        "--output-dir", default="outputs/grouped_additive_model/02_grouped_additive"
    )
    args = parser.parse_args()
    summary = run_grouped_additive(
        args.repo_root, data_root=args.data_root, output_dir=args.output_dir
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
