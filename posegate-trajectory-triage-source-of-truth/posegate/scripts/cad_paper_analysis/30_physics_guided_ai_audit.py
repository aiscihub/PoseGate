#!/usr/bin/env python3
"""Nested protein-held-out audit of physics-guided early-trajectory models.

All predictors end at 20 ns. The outcome is median corrected ligand RMSD over
(70,100] ns. Hyperparameters are selected using only training proteins inside
each outer leave-one-protein-out fold.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import SplineTransformer, StandardScaler


CHECKPOINTS = (5, 10, 15, 20)
ENDPOINT_FEATURES = (
    "pose_rmsd_mean_18_20_A",
    "centroid_displacement_mean_18_20_A",
)
TEMPORAL_FEATURES = tuple(
    feature
    for checkpoint in CHECKPOINTS
    for feature in (
        f"pose_rmsd_mean_{checkpoint - 2}_{checkpoint}_A",
        f"centroid_displacement_mean_{checkpoint - 2}_{checkpoint}_A",
    )
)
DERIVED_FEATURES = (
    "rmsd_change_5_20_A",
    "rmsd_change_15_20_A",
    "rmsd_acceleration_A",
    "centroid_change_5_20_A",
    "centroid_change_15_20_A",
    "centroid_acceleration_A",
)
ALL_FEATURES = TEMPORAL_FEATURES + DERIVED_FEATURES


def normalize_input(table: pd.DataFrame) -> pd.DataFrame:
    """Accept the canonical master schema while preserving legacy compatibility."""
    if "source_of_truth_version" not in table.columns:
        return table
    rename = {
        "measurement_status": "status",
        "complete_100ns_coverage": "complete_70_100_coverage",
    }
    for checkpoint in (5, 10, 15, 20, 30):
        rename[f"corrected_pose_rmsd_mean_{checkpoint - 2}_{checkpoint}_A"] = (
            f"pose_rmsd_mean_{checkpoint - 2}_{checkpoint}_A"
        )
        rename[f"corrected_centroid_displacement_mean_{checkpoint - 2}_{checkpoint}_A"] = (
            f"centroid_displacement_mean_{checkpoint - 2}_{checkpoint}_A"
        )
    return table.rename(columns=rename)


def add_temporal_features(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    r5, r10, r15, r20 = [f"pose_rmsd_mean_{c - 2}_{c}_A" for c in CHECKPOINTS]
    d5, d10, d15, d20 = [
        f"centroid_displacement_mean_{c - 2}_{c}_A" for c in CHECKPOINTS
    ]
    table["rmsd_change_5_20_A"] = table[r20] - table[r5]
    table["rmsd_change_15_20_A"] = table[r20] - table[r15]
    table["rmsd_acceleration_A"] = (table[r20] - table[r15]) - (table[r15] - table[r10])
    table["centroid_change_5_20_A"] = table[d20] - table[d5]
    table["centroid_change_15_20_A"] = table[d20] - table[d15]
    table["centroid_acceleration_A"] = (table[d20] - table[d15]) - (table[d15] - table[d10])
    return table


def model_candidates() -> dict[str, list[tuple[str, object]]]:
    endpoint_logistic = Pipeline(
        [("scale", StandardScaler()), ("model", LogisticRegression(class_weight="balanced", C=1.0, max_iter=2000))]
    )
    temporal_logistic = [
        (
            f"C={c}",
            Pipeline(
                [
                    ("impute", SimpleImputer(strategy="median")),
                    ("scale", StandardScaler()),
                    ("model", LogisticRegression(class_weight="balanced", C=c, max_iter=2000)),
                ]
            ),
        )
        for c in (0.03, 0.1, 0.3, 1.0, 3.0)
    ]
    spline = [
        (
            f"knots={knots},C={c}",
            Pipeline(
                [
                    ("impute", SimpleImputer(strategy="median")),
                    ("spline", SplineTransformer(n_knots=knots, degree=2, include_bias=False)),
                    ("scale", StandardScaler()),
                    ("model", LogisticRegression(class_weight="balanced", C=c, max_iter=3000)),
                ]
            ),
        )
        for knots in (3, 4)
        for c in (0.03, 0.1, 0.3, 1.0)
    ]
    hist = [
        (
            f"leaf={leaf},l2={l2}",
            Pipeline(
                [
                    ("impute", SimpleImputer(strategy="median")),
                    (
                        "model",
                        HistGradientBoostingClassifier(
                            max_iter=100,
                            learning_rate=0.05,
                            max_leaf_nodes=leaf,
                            min_samples_leaf=8,
                            l2_regularization=l2,
                            random_state=7,
                        ),
                    ),
                ]
            ),
        )
        for leaf in (2, 3, 5)
        for l2 in (0.0, 1.0, 5.0)
    ]
    forest = [
        (
            f"depth={depth},leaf={leaf}",
            Pipeline(
                [
                    ("impute", SimpleImputer(strategy="median")),
                    (
                        "model",
                        RandomForestClassifier(
                            n_estimators=200,
                            max_depth=depth,
                            min_samples_leaf=leaf,
                            max_features="sqrt",
                            class_weight="balanced_subsample",
                            random_state=7,
                            n_jobs=1,
                        ),
                    ),
                ]
            ),
        )
        for depth in (2, 3, None)
        for leaf in (3, 5, 8)
    ]
    return {
        "endpoint_logistic": [("fixed", endpoint_logistic)],
        "temporal_logistic": temporal_logistic,
        "endpoint_spline": spline,
        "temporal_hist_gradient_boosting": hist,
        "temporal_random_forest": forest,
    }


def feature_set(model_name: str) -> tuple[str, ...]:
    if model_name in {"endpoint_logistic", "endpoint_spline"}:
        return ENDPOINT_FEATURES
    return ALL_FEATURES


def inner_group_auc(
    table: pd.DataFrame, features: tuple[str, ...], estimator: object
) -> float:
    predictions = np.full(len(table), np.nan)
    for protein in sorted(table["protein"].unique()):
        test = table["protein"].eq(protein).to_numpy()
        train = ~test
        if table.loc[train, "late_pose_retained"].nunique() < 2:
            continue
        fitted = clone(estimator).fit(
            table.loc[train, list(features)], table.loc[train, "late_pose_retained"]
        )
        predictions[test] = fitted.predict_proba(table.loc[test, list(features)])[:, 1]
    valid = np.isfinite(predictions)
    if table.loc[valid, "late_pose_retained"].nunique() < 2:
        return -np.inf
    return roc_auc_score(table.loc[valid, "late_pose_retained"], predictions[valid])


def clustered_auc_interval(
    table: pd.DataFrame, score: np.ndarray, draws: int = 5000, seed: int = 7
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    proteins = table["protein"].unique()
    values = []
    scored = table.assign(_score=score)
    for _ in range(draws):
        sampled = rng.choice(proteins, size=len(proteins), replace=True)
        boot = pd.concat([scored[scored["protein"].eq(p)] for p in sampled])
        if boot["late_pose_retained"].nunique() < 2:
            continue
        values.append(roc_auc_score(boot["late_pose_retained"], boot["_score"]))
    return tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    full = add_temporal_features(normalize_input(pd.read_csv(args.input)))
    canonical_eligibility = (
        full["eligible_primary_cohort"].fillna(False)
        if "eligible_primary_cohort" in full.columns
        else pd.Series(True, index=full.index)
    )
    required = list(set(ALL_FEATURES + ENDPOINT_FEATURES))
    eligible = (
        canonical_eligibility
        &
        full["status"].eq("ok")
        & full["complete_70_100_coverage"].fillna(False)
        & full["late_pose_rmsd_median_70_100_A"].notna()
        & full[required].notna().all(axis=1)
    )
    table = full.loc[eligible].reset_index(drop=True).copy()
    table["late_pose_retained"] = (table["late_pose_rmsd_median_70_100_A"] < 3.0).astype(int)
    candidates = model_candidates()
    prediction_rows = table[["complex_id", "protein", "late_pose_retained"]].copy()
    selection_rows = []

    for model_name, configurations in candidates.items():
        features = feature_set(model_name)
        outer_scores = np.full(len(table), np.nan)
        for held_out in sorted(table["protein"].unique()):
            outer_test = table["protein"].eq(held_out).to_numpy()
            outer_train = table.loc[~outer_test].copy()
            ranked = []
            for config_name, estimator in configurations:
                ranked.append((inner_group_auc(outer_train, features, estimator), config_name, estimator))
            best_auc, best_name, best_estimator = max(ranked, key=lambda item: (item[0], item[1]))
            fitted = clone(best_estimator).fit(
                outer_train[list(features)], outer_train["late_pose_retained"]
            )
            outer_scores[outer_test] = fitted.predict_proba(table.loc[outer_test, list(features)])[:, 1]
            selection_rows.append(
                {
                    "model": model_name,
                    "held_out_protein": held_out,
                    "selected_configuration": best_name,
                    "inner_grouped_auroc": best_auc,
                    "outer_train_n": len(outer_train),
                    "feature_n": len(features),
                }
            )
        prediction_rows[f"score_{model_name}"] = outer_scores

    prediction_rows["score_raw_corrected_rmsd"] = -table[ENDPOINT_FEATURES[0]].to_numpy()
    metrics = []
    for column in [c for c in prediction_rows if c.startswith("score_")]:
        score = prediction_rows[column].to_numpy()
        auc = roc_auc_score(table["late_pose_retained"], score)
        ci = clustered_auc_interval(table, score)
        decision = (score >= 0.5).astype(int) if column != "score_raw_corrected_rmsd" else np.full(len(table), -1)
        metrics.append(
            {
                "model": column.removeprefix("score_"),
                "n": len(table),
                "proteins": table["protein"].nunique(),
                "auroc": auc,
                "clustered_ci_low": ci[0],
                "clustered_ci_high": ci[1],
                "balanced_accuracy_at_0.5": (
                    balanced_accuracy_score(table["late_pose_retained"], decision)
                    if column != "score_raw_corrected_rmsd"
                    else np.nan
                ),
            }
        )

    metrics_table = pd.DataFrame(metrics).sort_values("auroc", ascending=False)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    prediction_rows.to_csv(args.output_dir / "nested_lopo_predictions.csv", index=False)
    pd.DataFrame(selection_rows).to_csv(args.output_dir / "nested_model_selections.csv", index=False)
    metrics_table.to_csv(args.output_dir / "model_comparison_metrics.csv", index=False)
    print(metrics_table.to_string(index=False))


if __name__ == "__main__":
    main()
