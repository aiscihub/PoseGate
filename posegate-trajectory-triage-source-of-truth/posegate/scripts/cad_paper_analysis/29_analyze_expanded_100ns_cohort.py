#!/usr/bin/env python3
"""Audit primary same-trajectory results on the expanded raw-data cohort."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, recall_score, roc_auc_score
from sklearn.preprocessing import StandardScaler


FEATURES = (
    "pose_rmsd_mean_18_20_A",
    "centroid_displacement_mean_18_20_A",
)


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


def cluster_bootstrap_auc(
    table: pd.DataFrame, score_col: str, draws: int = 5000, seed: int = 7
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    proteins = table["protein"].unique()
    values = []
    for _ in range(draws):
        sampled = rng.choice(proteins, size=len(proteins), replace=True)
        boot = pd.concat([table[table["protein"].eq(p)] for p in sampled])
        if boot["late_pose_retained"].nunique() < 2:
            continue
        values.append(roc_auc_score(boot["late_pose_retained"], boot[score_col]))
    return tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    full = normalize_input(pd.read_csv(args.input))
    canonical_eligibility = (
        full["eligible_primary_cohort"].fillna(False)
        if "eligible_primary_cohort" in full.columns
        else pd.Series(True, index=full.index)
    )
    eligible = (
        canonical_eligibility
        &
        full["status"].eq("ok")
        & full["complete_70_100_coverage"].fillna(False)
        & full[list(FEATURES)].notna().all(axis=1)
        & full["late_pose_rmsd_median_70_100_A"].notna()
    )
    table = full.loc[eligible].copy()
    table["late_pose_retained"] = (
        table["late_pose_rmsd_median_70_100_A"] < 3.0
    ).astype(int)
    table["corrected_rmsd_score"] = -table[FEATURES[0]]

    table["lopo_score_retained"] = np.nan
    for held_out in sorted(table["protein"].unique()):
        test = table["protein"].eq(held_out)
        train = ~test
        scaler = StandardScaler()
        x_train = scaler.fit_transform(table.loc[train, list(FEATURES)])
        model = LogisticRegression(
            penalty="l2",
            C=1.0,
            fit_intercept=True,
            class_weight="balanced",
            solver="lbfgs",
            max_iter=100,
            tol=1e-4,
            random_state=7,
        )
        model.fit(x_train, table.loc[train, "late_pose_retained"])
        table.loc[test, "lopo_score_retained"] = model.predict_proba(
            scaler.transform(table.loc[test, list(FEATURES)])
        )[:, 1]

    table["lopo_continue_p05"] = (table["lopo_score_retained"] >= 0.5).astype(int)
    raw_auc = roc_auc_score(table["late_pose_retained"], table["corrected_rmsd_score"])
    raw_ci = cluster_bootstrap_auc(table, "corrected_rmsd_score")
    model_auc = roc_auc_score(table["late_pose_retained"], table["lopo_score_retained"])
    model_ci = cluster_bootstrap_auc(table, "lopo_score_retained")
    y = table["late_pose_retained"]
    decision = table["lopo_continue_p05"]
    stopped = int((decision == 0).sum())
    false_stops = int(((decision == 0) & (y == 1)).sum())
    nominal_saving = stopped * 80.0 / (len(table) * 100.0)

    metrics = pd.DataFrame(
        [
            {
                "cohort_n": len(table),
                "protein_n": table["protein"].nunique(),
                "retained_n": int(y.sum()),
                "non_retained_n": int((y == 0).sum()),
                "raw_corrected_rmsd_auroc": raw_auc,
                "raw_corrected_rmsd_ci_low": raw_ci[0],
                "raw_corrected_rmsd_ci_high": raw_ci[1],
                "two_feature_lopo_auroc": model_auc,
                "two_feature_lopo_ci_low": model_ci[0],
                "two_feature_lopo_ci_high": model_ci[1],
                "balanced_accuracy_p05": balanced_accuracy_score(y, decision),
                "retained_sensitivity_p05": recall_score(y, decision, pos_label=1),
                "non_retained_specificity_p05": recall_score(y, decision, pos_label=0),
                "stopped_n_p05": stopped,
                "false_stops_n_p05": false_stops,
                "nominal_compute_saving_p05": nominal_saving,
            }
        ]
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "expanded_cohort_predictions.csv", index=False)
    metrics.to_csv(args.output_dir / "expanded_cohort_metrics.csv", index=False)
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
