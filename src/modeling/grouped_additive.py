"""Training-only neural grouped additive model with exact group explanations."""

from __future__ import annotations

import copy
import hashlib
import json
import platform
import random
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
import torch
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from torch import Tensor, nn
from torch.utils.data import DataLoader, TensorDataset

from .common import FrozenTrainingData, load_frozen_training_data
from .logistic_baseline import (
    calibration_intercept_slope,
    classification_metrics,
    select_mcc_threshold,
)


@dataclass(frozen=True)
class GroupedAdditiveSettings:
    hidden_width: int = 8
    hidden_layers: int = 1
    activation: str = "relu"
    learning_rate: float = 0.001
    weight_decay: float = 0.0001
    dropout: float = 0.0
    batch_size: int = 512
    maximum_epochs: int = 200
    patience: int = 20

    @property
    def name(self) -> str:
        return (
            f"w{self.hidden_width}_l{self.hidden_layers}_{self.activation}_"
            f"lr{self.learning_rate:g}_wd{self.weight_decay:g}_d{self.dropout:g}"
        )


COMPACT_SEARCH = (
    GroupedAdditiveSettings(hidden_width=4, activation="relu"),
    GroupedAdditiveSettings(hidden_width=8, activation="relu"),
    GroupedAdditiveSettings(hidden_width=8, activation="tanh"),
    GroupedAdditiveSettings(hidden_width=16, activation="relu", dropout=0.1),
)

RECONSTRUCTION_TOLERANCE = 1e-5


def set_deterministic_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)


def build_preprocessor(data: FrozenTrainingData) -> ColumnTransformer:
    from sklearn.pipeline import Pipeline

    continuous = Pipeline(
        [("imputer", SimpleImputer(strategy="median")), ("scaler", StandardScaler())]
    )
    binary = Pipeline([("imputer", SimpleImputer(strategy="most_frequent"))])
    return ColumnTransformer(
        [
            ("continuous", continuous, list(data.continuous_features)),
            ("binary", binary, list(data.binary_features)),
        ],
        remainder="drop",
        verbose_feature_names_out=False,
    )


def group_indices(
    transformed_feature_names: list[str], evidence_groups: dict[str, tuple[str, ...]]
) -> dict[str, tuple[int, ...]]:
    positions = {name: index for index, name in enumerate(transformed_feature_names)}
    resolved: dict[str, tuple[int, ...]] = {}
    for group, features in evidence_groups.items():
        missing = sorted(set(features) - set(positions))
        if missing:
            raise ValueError(f"Transformed data is missing {group} features: {missing}")
        resolved[group] = tuple(positions[feature] for feature in features)
    flattened = [index for indexes in resolved.values() for index in indexes]
    if len(flattened) != len(set(flattened)):
        raise ValueError("A transformed predictor belongs to multiple evidence groups")
    if set(flattened) != set(range(len(transformed_feature_names))):
        raise ValueError("Not every transformed predictor belongs to an evidence group")
    return resolved


class GroupSubnetwork(nn.Module):
    def __init__(self, input_width: int, settings: GroupedAdditiveSettings) -> None:
        super().__init__()
        activation_class = nn.ReLU if settings.activation == "relu" else nn.Tanh
        layers: list[nn.Module] = []
        width = input_width
        for _ in range(settings.hidden_layers):
            layers.extend([nn.Linear(width, settings.hidden_width), activation_class()])
            if settings.dropout:
                layers.append(nn.Dropout(settings.dropout))
            width = settings.hidden_width
        layers.append(nn.Linear(width, 1))
        self.network = nn.Sequential(*layers)

    def forward(self, values: Tensor) -> Tensor:
        return self.network(values).squeeze(1)


class GroupedAdditiveNetwork(nn.Module):
    """No cross-group pathway: logit equals intercept plus four group outputs."""

    def __init__(
        self, group_indexes: dict[str, tuple[int, ...]], settings: GroupedAdditiveSettings
    ) -> None:
        super().__init__()
        self.group_names = tuple(group_indexes)
        self.group_indexes = {name: tuple(values) for name, values in group_indexes.items()}
        self.subnetworks = nn.ModuleDict(
            {
                name: GroupSubnetwork(len(self.group_indexes[name]), settings)
                for name in self.group_names
            }
        )
        self.intercept = nn.Parameter(torch.zeros(()))

    def forward(self, values: Tensor) -> tuple[Tensor, Tensor]:
        contributions = torch.stack(
            [
                self.subnetworks[name](values[:, self.group_indexes[name]])
                for name in self.group_names
            ],
            dim=1,
        )
        logit = self.intercept + contributions.sum(dim=1)
        return logit, contributions


def _internal_gene_split(genes: pd.Series, seed: int) -> tuple[np.ndarray, np.ndarray]:
    unique = np.array(sorted(genes.astype(str).unique()))
    rng = np.random.default_rng(seed)
    rng.shuffle(unique)
    validation_genes = set(unique[: max(1, int(round(len(unique) * 0.1)))])
    validation = genes.astype(str).isin(validation_genes).to_numpy()
    fitting = ~validation
    if set(genes[fitting]).intersection(set(genes[validation])):
        raise RuntimeError("Internal early-stopping gene leakage detected")
    return fitting, validation


def _tensor(values: np.ndarray) -> Tensor:
    return torch.as_tensor(np.asarray(values, dtype=np.float32))


def _train_epochs(
    model: GroupedAdditiveNetwork,
    x: np.ndarray,
    y: np.ndarray,
    settings: GroupedAdditiveSettings,
    seed: int,
    epochs: int,
) -> list[float]:
    set_deterministic_seed(seed)
    generator = torch.Generator().manual_seed(seed)
    loader = DataLoader(
        TensorDataset(_tensor(x), _tensor(y)),
        batch_size=settings.batch_size,
        shuffle=True,
        generator=generator,
    )
    optimiser = torch.optim.Adam(
        model.parameters(), lr=settings.learning_rate, weight_decay=settings.weight_decay
    )
    loss_function = nn.BCEWithLogitsLoss()
    history: list[float] = []
    for _ in range(epochs):
        model.train()
        losses: list[float] = []
        for batch_x, batch_y in loader:
            optimiser.zero_grad(set_to_none=True)
            logits, _ = model(batch_x)
            loss = loss_function(logits, batch_y)
            loss.backward()
            optimiser.step()
            losses.append(float(loss.detach()))
        history.append(float(np.mean(losses)))
    return history


def _find_best_epoch(
    model: GroupedAdditiveNetwork,
    x_fit: np.ndarray,
    y_fit: np.ndarray,
    x_stop: np.ndarray,
    y_stop: np.ndarray,
    settings: GroupedAdditiveSettings,
    seed: int,
) -> tuple[int, list[dict[str, float]]]:
    set_deterministic_seed(seed)
    generator = torch.Generator().manual_seed(seed)
    loader = DataLoader(
        TensorDataset(_tensor(x_fit), _tensor(y_fit)),
        batch_size=settings.batch_size,
        shuffle=True,
        generator=generator,
    )
    optimiser = torch.optim.Adam(
        model.parameters(), lr=settings.learning_rate, weight_decay=settings.weight_decay
    )
    loss_function = nn.BCEWithLogitsLoss()
    stop_x, stop_y = _tensor(x_stop), _tensor(y_stop)
    best_loss = float("inf")
    best_epoch = 1
    best_state = copy.deepcopy(model.state_dict())
    without_improvement = 0
    history: list[dict[str, float]] = []
    for epoch in range(1, settings.maximum_epochs + 1):
        model.train()
        train_losses: list[float] = []
        for batch_x, batch_y in loader:
            optimiser.zero_grad(set_to_none=True)
            logits, _ = model(batch_x)
            loss = loss_function(logits, batch_y)
            loss.backward()
            optimiser.step()
            train_losses.append(float(loss.detach()))
        model.eval()
        with torch.no_grad():
            stopping_logits, _ = model(stop_x)
            stopping_loss = float(loss_function(stopping_logits, stop_y))
        history.append(
            {"epoch": epoch, "train_loss": float(np.mean(train_losses)), "stopping_loss": stopping_loss}
        )
        if stopping_loss < best_loss - 1e-6:
            best_loss = stopping_loss
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            without_improvement = 0
        else:
            without_improvement += 1
            if without_improvement >= settings.patience:
                break
    model.load_state_dict(best_state)
    return best_epoch, history


def _predict(model: GroupedAdditiveNetwork, x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    with torch.no_grad():
        logits, contributions = model(_tensor(x))
        probability = torch.sigmoid(logits)
    return (
        logits.numpy().astype(float),
        probability.numpy().astype(float),
        contributions.numpy().astype(float),
    )


def _fit_outer_fold(
    data: FrozenTrainingData,
    fold: int,
    settings: GroupedAdditiveSettings,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any], pd.DataFrame]:
    frame = data.frame
    outer_train = frame["cv_fold"].ne(fold).to_numpy()
    outer_holdout = ~outer_train
    train_frame = frame.loc[outer_train].reset_index(drop=True)
    fit_mask, stopping_mask = _internal_gene_split(train_frame[data.gene], seed + fold)

    inner_preprocessor = build_preprocessor(data)
    inner_preprocessor.fit(train_frame.loc[fit_mask, list(data.features)])
    x_fit = inner_preprocessor.transform(train_frame.loc[fit_mask, list(data.features)])
    x_stop = inner_preprocessor.transform(train_frame.loc[stopping_mask, list(data.features)])
    names = list(inner_preprocessor.get_feature_names_out())
    indexes = group_indices(names, data.evidence_groups)
    set_deterministic_seed(seed + fold)
    model = GroupedAdditiveNetwork(indexes, settings)
    best_epoch, history = _find_best_epoch(
        model,
        x_fit,
        train_frame.loc[fit_mask, data.label].to_numpy(dtype=float),
        x_stop,
        train_frame.loc[stopping_mask, data.label].to_numpy(dtype=float),
        settings,
        seed + fold,
    )

    preprocessor = build_preprocessor(data)
    x_train = preprocessor.fit_transform(frame.loc[outer_train, list(data.features)])
    x_holdout = preprocessor.transform(frame.loc[outer_holdout, list(data.features)])
    names = list(preprocessor.get_feature_names_out())
    indexes = group_indices(names, data.evidence_groups)
    set_deterministic_seed(seed + 1000 + fold)
    model = GroupedAdditiveNetwork(indexes, settings)
    _train_epochs(
        model,
        x_train,
        frame.loc[outer_train, data.label].to_numpy(dtype=float),
        settings,
        seed + 1000 + fold,
        best_epoch,
    )
    logits, probability, contributions = _predict(model, x_holdout)
    rows = frame.loc[outer_holdout, [data.variant_key, data.gene, data.label, "cv_fold"]].copy()
    rows["logit"] = logits
    rows["probability"] = probability
    rows["intercept"] = float(model.intercept.detach())
    for column, group in enumerate(model.group_names):
        rows[f"contribution_{group}"] = contributions[:, column]
    reconstructed = rows["intercept"].to_numpy() + contributions.sum(axis=1)
    rows["reconstructed_logit"] = reconstructed
    rows["reconstruction_error"] = np.abs(logits - reconstructed)
    if rows["reconstruction_error"].max() > RECONSTRUCTION_TOLERANCE:
        raise RuntimeError("Saved group contributions do not reconstruct the model logit")
    metrics = classification_metrics(
        rows[data.label].to_numpy(dtype=int), probability, threshold=0.5
    )
    metrics.update(
        {
            "configuration": settings.name,
            "fold": fold,
            "best_epoch": best_epoch,
            "train_rows": int(outer_train.sum()),
            "holdout_rows": int(outer_holdout.sum()),
            "max_reconstruction_error": float(rows["reconstruction_error"].max()),
        }
    )
    history_frame = pd.DataFrame(history)
    history_frame.insert(0, "fold", fold)
    history_frame.insert(0, "configuration", settings.name)
    return rows, metrics, history_frame


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_grouped_additive(
    repo_root: str | Path,
    *,
    data_root: str | Path | None = None,
    output_dir: str | Path = "outputs/grouped_additive_model/02_grouped_additive",
    search: tuple[GroupedAdditiveSettings, ...] = COMPACT_SEARCH,
) -> dict[str, Any]:
    started = time.perf_counter()
    repo_root = Path(repo_root).resolve()
    output_dir = repo_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    data = load_frozen_training_data(repo_root, data_root=data_root)
    seed = int(data.protocol["random_seed"])

    candidate_predictions: dict[str, pd.DataFrame] = {}
    all_fold_metrics: list[dict[str, Any]] = []
    all_histories: list[pd.DataFrame] = []
    settings_by_name = {settings.name: settings for settings in search}
    for settings in search:
        fold_predictions: list[pd.DataFrame] = []
        for fold in range(5):
            predictions, metrics, history = _fit_outer_fold(data, fold, settings, seed)
            fold_predictions.append(predictions)
            all_fold_metrics.append(metrics)
            all_histories.append(history)
        combined = pd.concat(fold_predictions).sort_index(kind="stable")
        if len(combined) != len(data.frame) or combined[data.variant_key].nunique() != len(data.frame):
            raise RuntimeError("Out-of-fold predictions are not complete and unique")
        candidate_predictions[settings.name] = combined

    fold_metrics = pd.DataFrame(all_fold_metrics)
    search_summary = (
        fold_metrics.groupby("configuration", as_index=False)
        .agg(mean_fold_mcc=("mcc", "mean"), sd_fold_mcc=("mcc", "std"), median_best_epoch=("best_epoch", "median"))
        .sort_values(["mean_fold_mcc", "configuration"], ascending=[False, True], kind="stable")
    )
    selected_name = str(search_summary.iloc[0]["configuration"])
    selected_settings = settings_by_name[selected_name]
    predictions = candidate_predictions[selected_name].copy()
    threshold, _ = select_mcc_threshold(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy()
    )
    predictions["prediction_at_0_5"] = (predictions["probability"] >= 0.5).astype(int)
    predictions["selected_threshold"] = threshold
    predictions["prediction_at_selected_threshold"] = (
        predictions["probability"] >= threshold
    ).astype(int)
    calibration_intercept, calibration_slope = calibration_intercept_slope(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy()
    )
    pooled_half = classification_metrics(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy(), 0.5
    )
    pooled_selected = classification_metrics(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy(), threshold
    )

    final_preprocessor = build_preprocessor(data)
    x_full = final_preprocessor.fit_transform(data.frame[list(data.features)])
    names = list(final_preprocessor.get_feature_names_out())
    indexes = group_indices(names, data.evidence_groups)
    final_epoch = int(round(float(search_summary.iloc[0]["median_best_epoch"])))
    set_deterministic_seed(seed)
    final_model = GroupedAdditiveNetwork(indexes, selected_settings)
    final_history = _train_epochs(
        final_model,
        x_full,
        data.frame[data.label].to_numpy(dtype=float),
        selected_settings,
        seed,
        final_epoch,
    )
    joblib.dump(final_preprocessor, output_dir / "fitted_training_preprocessor.joblib")
    torch.save(
        {
            "state_dict": final_model.state_dict(),
            "group_indexes": indexes,
            "group_names": final_model.group_names,
            "settings": asdict(selected_settings),
            "feature_names": names,
        },
        output_dir / "fitted_training_model.pt",
    )

    explanation_columns = [f"contribution_{name}" for name in data.evidence_groups]
    summary = {
        "status": "training_only_complete",
        "model": "neural_grouped_additive",
        "equation": "logit = intercept + conservation + gene_constraint + amino_acid_properties + protein_context",
        "training_rows": len(data.frame),
        "feature_count": len(data.features),
        "evidence_groups": list(data.evidence_groups),
        "selected_configuration": asdict(selected_settings),
        "selected_configuration_name": selected_name,
        "selected_final_epoch": final_epoch,
        "selected_threshold": threshold,
        "pooled_oof_metrics_at_0_5": pooled_half,
        "pooled_oof_metrics_at_selected_threshold": pooled_selected,
        "calibration_intercept": calibration_intercept,
        "calibration_slope": calibration_slope,
        "maximum_oof_reconstruction_error": float(predictions["reconstruction_error"].max()),
        "reconstruction_tolerance": RECONSTRUCTION_TOLERANCE,
        "explanation_columns": explanation_columns,
        "validation_loaded": False,
        "test_loaded": False,
        "random_seed": seed,
        "python": sys.version,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "torch": torch.__version__,
        "platform": platform.platform(),
        "runtime_seconds": float(time.perf_counter() - started),
    }
    predictions.to_csv(output_dir / "oof_predictions_and_group_contributions.csv", index=False)
    fold_metrics.to_csv(output_dir / "fold_metrics_all_candidates.csv", index=False)
    search_summary.to_csv(output_dir / "hyperparameter_search.csv", index=False)
    pd.concat(all_histories, ignore_index=True).to_csv(output_dir / "training_histories.csv", index=False)
    pd.DataFrame({"epoch": np.arange(1, len(final_history) + 1), "training_loss": final_history}).to_csv(
        output_dir / "final_training_history.csv", index=False
    )
    with (output_dir / "run_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    manifest = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "output_sha256.csv":
            manifest.append({"file": path.name, "bytes": path.stat().st_size, "sha256": _sha256(path)})
    pd.DataFrame(manifest).to_csv(output_dir / "output_sha256.csv", index=False)
    return summary
