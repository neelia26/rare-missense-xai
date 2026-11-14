"""Resolve and validate the frozen grouped-additive modelling inputs."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


class FrozenInputError(ValueError):
    """Raised when a frozen modelling input or governance rule is violated."""


@dataclass(frozen=True)
class FrozenTrainingData:
    frame: pd.DataFrame
    features: tuple[str, ...]
    evidence_groups: dict[str, tuple[str, ...]]
    continuous_features: tuple[str, ...]
    binary_features: tuple[str, ...]
    variant_key: str
    gene: str
    label: str
    input_config: dict[str, Any]
    protocol: dict[str, Any]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def _require_hash(path: Path, expected: str) -> None:
    if not path.is_file():
        raise FrozenInputError(f"Frozen input is missing: {path}")
    actual = sha256_file(path)
    if actual != expected:
        raise FrozenInputError(
            f"Frozen input hash mismatch for {path}: expected {expected}, got {actual}"
        )


def _flatten_groups(groups: dict[str, list[str]]) -> list[str]:
    return [feature for features in groups.values() for feature in features]


def _feature_types(features: list[str]) -> tuple[list[str], list[str]]:
    explicit_binary = {"aa_charge_changed", "aa_polarity_class_changed"}
    binary = [
        feature
        for feature in features
        if feature in explicit_binary
        or feature.startswith("interpro_type_")
        or feature.startswith("uniprot_feature_")
    ]
    continuous = [feature for feature in features if feature not in set(binary)]
    return continuous, binary


def _validate_feature_governance(config: dict[str, Any]) -> tuple[list[str], list[str], list[str]]:
    groups = config["evidence_groups"]
    features = _flatten_groups(groups)
    expected = int(config["feature_count"])
    if len(features) != expected:
        raise FrozenInputError(
            f"Expected {expected} primary predictors, but configuration lists {len(features)}"
        )
    if len(features) != len(set(features)):
        raise FrozenInputError("A primary predictor belongs to more than one evidence group")

    forbidden = {
        "conservation_gerp_rs",
        "aa_polarity_delta",
        "gnomad_global_af",
        "gnomad_max_population_af",
        "binary_label",
        "variant_key",
        "reconciled_gene_symbol",
    }
    forbidden_present = sorted(forbidden.intersection(features))
    if forbidden_present:
        raise FrozenInputError(
            f"Forbidden primary predictors are present: {', '.join(forbidden_present)}"
        )
    if any(feature.startswith("predictor_") for feature in features):
        raise FrozenInputError("Composite computational predictors are forbidden in primary input")

    continuous, binary = _feature_types(features)
    return features, continuous, binary


def _validate_training_rows(
    frame: pd.DataFrame,
    *,
    variant_key: str,
    gene: str,
    label: str,
    features: list[str],
    expected_rows: int,
) -> None:
    required = [variant_key, gene, label, *features]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise FrozenInputError(f"Training data is missing columns: {missing}")
    if len(frame) != expected_rows:
        raise FrozenInputError(
            f"Expected {expected_rows} training rows, found {len(frame)}"
        )
    if frame[variant_key].isna().any() or frame[variant_key].astype(str).str.strip().eq("").any():
        raise FrozenInputError("Blank variant keys are not permitted")
    if frame[gene].isna().any() or frame[gene].astype(str).str.strip().eq("").any():
        raise FrozenInputError("Blank reconciled genes are not permitted")
    if frame[variant_key].duplicated().any():
        raise FrozenInputError("Duplicate variant keys are not permitted")
    labels = set(frame[label].dropna().astype(int).unique())
    if labels != {0, 1}:
        raise FrozenInputError(f"Expected binary labels {{0, 1}}, found {sorted(labels)}")


def _attach_fixed_folds(
    frame: pd.DataFrame,
    fold_path: Path,
    *,
    variant_key: str,
    gene: str,
) -> pd.DataFrame:
    folds = pd.read_csv(fold_path)
    required = {variant_key, gene, "cv_fold"}
    if not required.issubset(folds.columns):
        raise FrozenInputError(f"Fold assignment is missing columns: {sorted(required - set(folds.columns))}")
    if folds[variant_key].duplicated().any():
        raise FrozenInputError("Fold assignment contains duplicate variant keys")
    if folds[variant_key].isna().any() or folds[gene].isna().any():
        raise FrozenInputError("Fold assignment contains blank keys or genes")
    if folds.groupby(gene)["cv_fold"].nunique().max() != 1:
        raise FrozenInputError("At least one gene appears in more than one CV fold")

    merged = frame.merge(
        folds[[variant_key, gene, "cv_fold"]],
        on=[variant_key, gene],
        how="left",
        validate="one_to_one",
        sort=False,
    )
    if merged["cv_fold"].isna().any():
        raise FrozenInputError("Some training variants have no fixed CV-fold assignment")
    observed = set(merged["cv_fold"].astype(int).unique())
    if observed != {0, 1, 2, 3, 4}:
        raise FrozenInputError(f"Expected folds 0-4, found {sorted(observed)}")
    merged["cv_fold"] = merged["cv_fold"].astype(int)
    return merged


def load_frozen_training_data(
    repo_root: str | Path,
    *,
    input_config_path: str | Path = "config/grouped_additive_primary_input.json",
    protocol_path: str | Path = "config/grouped_additive_experiment_protocol.json",
    data_root: str | Path | None = None,
    verify_reserved_hashes: bool = True,
) -> FrozenTrainingData:
    """Load training data only; validation and test contents are never read."""

    repo_root = Path(repo_root).resolve()
    data_root = Path(data_root).resolve() if data_root else repo_root
    input_config = load_json(repo_root / input_config_path)
    protocol = load_json(repo_root / protocol_path)

    features, continuous, binary = _validate_feature_governance(input_config)
    dataset = input_config["dataset"]

    paths: dict[str, Path] = {}
    for key in ["cohort", "training", "validation", "test", "split_assignments"]:
        paths[key] = data_root / dataset[f"{key}_path"]
    _require_hash(paths["training"], dataset["training_sha256"])
    _require_hash(paths["cohort"], dataset["cohort_sha256"])
    _require_hash(paths["split_assignments"], dataset["split_assignments_sha256"])
    if verify_reserved_hashes:
        # Byte-level verification is allowed; no reserved labels or rows are parsed.
        _require_hash(paths["validation"], dataset["validation_sha256"])
        _require_hash(paths["test"], dataset["test_sha256"])

    columns = input_config["columns"]
    variant_key = columns["variant_key"]
    gene = columns["gene"]
    label = columns["label"]
    frame = pd.read_csv(paths["training"], usecols=[variant_key, gene, label, *features])
    _validate_training_rows(
        frame,
        variant_key=variant_key,
        gene=gene,
        label=label,
        features=features,
        expected_rows=int(dataset["training_rows"]),
    )

    fold_path = repo_root / protocol["cross_validation"]["fold_assignment_path"]
    _require_hash(fold_path, protocol["cross_validation"]["fold_assignment_sha256"])
    frame = _attach_fixed_folds(frame, fold_path, variant_key=variant_key, gene=gene)

    return FrozenTrainingData(
        frame=frame,
        features=tuple(features),
        evidence_groups={name: tuple(values) for name, values in input_config["evidence_groups"].items()},
        continuous_features=tuple(continuous),
        binary_features=tuple(binary),
        variant_key=variant_key,
        gene=gene,
        label=label,
        input_config=input_config,
        protocol=protocol,
    )
