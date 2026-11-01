#!/usr/bin/env python3
"""Dedicated entry point for Phase A cohort and split validation."""

from src.pre_modeling_audit.run_pre_modeling_audit import inventory_and_leakage, validate_cohort, validate_split


def main() -> None:
    cohort, _ = validate_cohort()
    validate_split(cohort)
    inventory_and_leakage(cohort)
    print("Final cohort and locked split: PASS")


if __name__ == "__main__":
    main()
