"""Training-only gene-grouped logistic-regression baseline."""

from __future__ import annotations

import json
import platform
import sys
import time
import hashlib
from pathlib import Path
from typing import Any

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    matthews_corrcoef,
    precision_score,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from .common import FrozenTrainingData, load_frozen_training_data


DEFAULT_C_GRID = (0.01, 0.1, 1.0, 10.0)


def build_pipeline(data: FrozenTrainingData, *, c_value: float, seed: int) -> Pipeline:
    continuous = Pipeline(
        [
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )
    binary = Pipeline([("imputer", SimpleImputer(strategy="most_frequent"))])
    preprocessing = ColumnTransformer(
        [
            ("continuous", continuous, list(data.continuous_features)),
            ("binary", binary, list(data.binary_features)),
        ],
        remainder="drop",
        verbose_feature_names_out=False,
    )
    model = LogisticRegression(
        C=float(c_value),
        solver="liblinear",
        max_iter=3000,
        random_state=int(seed),
    )
    return Pipeline([("preprocessing", preprocessing), ("model", model)])


def classification_metrics(y_true: np.ndarray, probability: np.ndarray, threshold: float) -> dict[str, float]:
    prediction = (probability >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, prediction, labels=[0, 1]).ravel()
    return {
        "threshold": float(threshold),
        "mcc": float(matthews_corrcoef(y_true, prediction)),
        "roc_auc": float(roc_auc_score(y_true, probability)),
        "pr_auc": float(average_precision_score(y_true, probability)),
        "brier": float(brier_score_loss(y_true, probability)),
        "accuracy": float(accuracy_score(y_true, prediction)),
        "sensitivity": float(tp / (tp + fn)) if tp + fn else float("nan"),
        "specificity": float(tn / (tn + fp)) if tn + fp else float("nan"),
        "precision": float(precision_score(y_true, prediction, zero_division=0)),
        "f1": float(f1_score(y_true, prediction, zero_division=0)),
        "tn": int(tn),
        "fp": int(fp),
        "fn": int(fn),
        "tp": int(tp),
    }


def select_mcc_threshold(y_true: np.ndarray, probability: np.ndarray) -> tuple[float, float]:
    thresholds = np.round(np.arange(0.05, 0.951, 0.005), 3)
    scored = [
        (float(matthews_corrcoef(y_true, probability >= threshold)), float(threshold))
        for threshold in thresholds
    ]
    best_mcc = max(score for score, _ in scored)
    tied = [threshold for score, threshold in scored if np.isclose(score, best_mcc)]
    return min(tied, key=lambda threshold: (abs(threshold - 0.5), threshold)), best_mcc


def calibration_intercept_slope(
    y_true: np.ndarray, probability: np.ndarray
) -> tuple[float, float]:
    clipped = np.clip(probability, 1e-8, 1 - 1e-8)
    logit = np.log(clipped / (1 - clipped)).reshape(-1, 1)
    calibrator = LogisticRegression(C=1e6, solver="lbfgs", max_iter=2000)
    calibrator.fit(logit, y_true)
    return float(calibrator.intercept_[0]), float(calibrator.coef_[0, 0])


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _fit_candidate(
    data: FrozenTrainingData,
    *,
    c_value: float,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frame = data.frame
    oof = np.full(len(frame), np.nan, dtype=float)
    fold_rows: list[dict[str, Any]] = []
    coefficient_rows: list[dict[str, Any]] = []

    for fold in range(5):
        train_mask = frame["cv_fold"].ne(fold).to_numpy()
        holdout_mask = ~train_mask
        train_genes = set(frame.loc[train_mask, data.gene])
        holdout_genes = set(frame.loc[holdout_mask, data.gene])
        if train_genes.intersection(holdout_genes):
            raise RuntimeError(f"Gene leakage detected in fold {fold}")

        pipeline = build_pipeline(data, c_value=c_value, seed=seed + fold)
        started = time.perf_counter()
        pipeline.fit(frame.loc[train_mask, list(data.features)], frame.loc[train_mask, data.label])
        elapsed = time.perf_counter() - started
        probability = pipeline.predict_proba(frame.loc[holdout_mask, list(data.features)])[:, 1]
        oof[holdout_mask] = probability
        metrics = classification_metrics(
            frame.loc[holdout_mask, data.label].to_numpy(dtype=int), probability, 0.5
        )
        metrics.update(
            {
                "C": float(c_value),
                "fold": int(fold),
                "train_rows": int(train_mask.sum()),
                "holdout_rows": int(holdout_mask.sum()),
                "train_genes": len(train_genes),
                "holdout_genes": len(holdout_genes),
                "fit_seconds": float(elapsed),
            }
        )
        fold_rows.append(metrics)

        transformed_names = pipeline.named_steps["preprocessing"].get_feature_names_out()
        coefficients = pipeline.named_steps["model"].coef_[0]
        for feature, coefficient in zip(transformed_names, coefficients, strict=True):
            coefficient_rows.append(
                {
                    "C": float(c_value),
                    "fold": int(fold),
                    "feature": str(feature),
                    "coefficient": float(coefficient),
                }
            )

    if np.isnan(oof).any():
        raise RuntimeError("Out-of-fold probabilities are incomplete")
    predictions = frame[[data.variant_key, data.gene, data.label, "cv_fold"]].copy()
    predictions["probability"] = oof
    predictions["prediction_at_0_5"] = (oof >= 0.5).astype(int)
    return predictions, pd.DataFrame(fold_rows), pd.DataFrame(coefficient_rows)


def _write_diagnostic_plot(predictions: pd.DataFrame, output_path: Path, label: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    for target, colour, name in [(0, "#5B6C8F", "Benign"), (1, "#7A4EAB", "Pathogenic")]:
        values = predictions.loc[predictions[label].eq(target), "probability"]
        ax.hist(values, bins=35, alpha=0.55, density=True, color=colour, label=name)
    ax.set_xlabel("Out-of-fold pathogenicity probability")
    ax.set_ylabel("Density")
    ax.set_title("Logistic-regression baseline: training-only OOF probabilities")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def run_logistic_baseline(
    repo_root: str | Path,
    *,
    data_root: str | Path | None = None,
    output_dir: str | Path = "outputs/grouped_additive_model/01_logistic_baseline",
    c_grid: tuple[float, ...] = DEFAULT_C_GRID,
) -> dict[str, Any]:
    run_started = time.perf_counter()
    repo_root = Path(repo_root).resolve()
    output_dir = repo_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    data = load_frozen_training_data(repo_root, data_root=data_root)
    seed = int(data.protocol["random_seed"])

    candidate_results: dict[float, tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = {}
    search_rows: list[dict[str, float]] = []
    for c_value in c_grid:
        predictions, fold_metrics, coefficients = _fit_candidate(
            data, c_value=float(c_value), seed=seed
        )
        candidate_results[float(c_value)] = (predictions, fold_metrics, coefficients)
        search_rows.append(
            {
                "C": float(c_value),
                "mean_fold_mcc": float(fold_metrics["mcc"].mean()),
                "sd_fold_mcc": float(fold_metrics["mcc"].std(ddof=1)),
                "pooled_oof_mcc_at_0_5": float(
                    matthews_corrcoef(
                        predictions[data.label], predictions["prediction_at_0_5"]
                    )
                ),
            }
        )

    search = pd.DataFrame(search_rows).sort_values(
        ["mean_fold_mcc", "C"], ascending=[False, True], kind="stable"
    )
    selected_c = float(search.iloc[0]["C"])
    predictions, fold_metrics, coefficients = candidate_results[selected_c]
    selected_threshold, _ = select_mcc_threshold(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy()
    )
    predictions["selected_threshold"] = selected_threshold
    predictions["prediction_at_selected_threshold"] = (
        predictions["probability"] >= selected_threshold
    ).astype(int)

    pooled_at_half = classification_metrics(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy(), 0.5
    )
    pooled_selected = classification_metrics(
        predictions[data.label].to_numpy(dtype=int),
        predictions["probability"].to_numpy(),
        selected_threshold,
    )
    calibration_intercept, calibration_slope = calibration_intercept_slope(
        predictions[data.label].to_numpy(dtype=int), predictions["probability"].to_numpy()
    )

    final_pipeline = build_pipeline(data, c_value=selected_c, seed=seed)
    final_pipeline.fit(data.frame[list(data.features)], data.frame[data.label])
    joblib.dump(final_pipeline, output_dir / "fitted_training_pipeline.joblib")
    final_feature_names = final_pipeline.named_steps["preprocessing"].get_feature_names_out()
    final_coefficients = pd.DataFrame(
        {
            "feature": final_feature_names,
            "coefficient": final_pipeline.named_steps["model"].coef_[0],
        }
    )
    final_coefficients.to_csv(output_dir / "fitted_training_coefficients.csv", index=False)
    summary = {
        "status": "training_only_complete",
        "model": "logistic_regression_baseline",
        "feature_count": len(data.features),
        "training_rows": len(data.frame),
        "selected_C": selected_c,
        "selected_threshold": selected_threshold,
        "calibration_intercept": calibration_intercept,
        "calibration_slope": calibration_slope,
        "pooled_oof_metrics_at_0_5": pooled_at_half,
        "pooled_oof_metrics_at_selected_threshold": pooled_selected,
        "validation_loaded": False,
        "test_loaded": False,
        "random_seed": seed,
        "python": sys.version,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scikit_learn": sklearn.__version__,
        "joblib": joblib.__version__,
        "platform": platform.platform(),
        "total_runtime_seconds": float(time.perf_counter() - run_started),
    }

    predictions.to_csv(output_dir / "oof_predictions.csv", index=False)
    fold_metrics.to_csv(output_dir / "fold_metrics.csv", index=False)
    coefficients.to_csv(output_dir / "fold_coefficients.csv", index=False)
    search.to_csv(output_dir / "hyperparameter_search.csv", index=False)
    with (output_dir / "run_summary.json").open("w", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True)
        stream.write("\n")
    _write_diagnostic_plot(
        predictions, output_dir / "oof_probability_distribution.png", data.label
    )
    manifest_rows = []
    for path in sorted(output_dir.iterdir()):
        if path.is_file() and path.name != "output_sha256.csv":
            manifest_rows.append(
                {"file": path.name, "bytes": path.stat().st_size, "sha256": _sha256(path)}
            )
    pd.DataFrame(manifest_rows).to_csv(output_dir / "output_sha256.csv", index=False)
    return summary
