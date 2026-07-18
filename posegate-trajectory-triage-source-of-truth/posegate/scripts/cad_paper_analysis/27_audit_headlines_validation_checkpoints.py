#!/usr/bin/env python3
"""Resolve manuscript TODO items A3, A4, and A5 from the locked audit table.

This script has three jobs:

1. A3 -- build a machine-readable headline registry and explain why the
   legacy headline numbers disagreed;
2. A4 -- audit the outer-LOPO model, preprocessing, regularization, decision
   threshold, policy use, and deployment status; and
3. A5 -- rerun every checkpoint model on the identical locked 52-trajectory
   cohort and freeze 20 ns for a future prospective shadow test.

Locked results are recomputed only from
``same_trajectory_source_of_truth.csv``. Legacy aggregate files are read only
to document provenance for superseded numbers; they are never used to compute
the locked results.

Outputs
-------
pipeline/training/outputs/yau_recapture_screen/cad_states/
same_trajectory_source_of_truth/
    headline_result_registry.csv
    legacy_headline_reconciliation.csv
    validation_procedure_audit.csv
    checkpoint_model_results.csv
    checkpoint_oof_predictions.csv
    analysis_lock_a3_a5.csv
    RESULTS_LOCK_A3_A5.md

Usage:
    python 27_audit_headlines_validation_checkpoints.py
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, recall_score
from sklearn.preprocessing import StandardScaler

from common import OUT_DIR


RESULT_DIR = OUT_DIR / "same_trajectory_source_of_truth"
SOURCE_TABLE = RESULT_DIR / "same_trajectory_source_of_truth.csv"
CLAIM_METRICS_TABLE = RESULT_DIR / "claim_figures" / "key_claim_metrics.csv"
CLAIM_CHECKPOINT_TABLE = (
    RESULT_DIR / "claim_figures" / "checkpoint_locked_common_cohort.csv"
)

LEGACY_MODEL_PREDICTIONS = (
    OUT_DIR
    / "same_trajectory_monitoring_validation"
    / "sametraj_model_comparison_predictions.csv"
)
LEGACY_CHECKPOINT_RESULTS = (
    OUT_DIR / "checkpoint_curve_and_baselines" / "checkpoint_time_curve.csv"
)
LEGACY_COORDINATE_RESULTS = (
    OUT_DIR / "checkpoint_curve_and_baselines" / "baseline_comparison.csv"
)
LEGACY_REPLICA_RESULTS = (
    OUT_DIR
    / "same_traj_vs_independent_replica"
    / "controlled_comparison_table.csv"
)
LEGACY_VALIDATION2_USABLE = (
    OUT_DIR
    / "pocket_aware_validation2_100ns"
    / "validation2_100ns_outcomes_usable.csv"
)

LABEL = "late_pose_retained_locked"
PROTEIN = "protein"
CHECKPOINTS_NS = (5, 10, 15, 20, 30)
N_BOOTSTRAP = 5000
SEED = 7

MODEL_PARAMETERS = {
    "penalty": "l2",
    "C": 1.0,
    "fit_intercept": True,
    "class_weight": "balanced",
    "solver": "lbfgs",
    "max_iter": 100,
    "tol": 1e-4,
    "random_state": 7,
}

OUTCOME_DEFINITION = (
    "late_pose_retained_rmsd_v1: complete nominal (70,100] ns coverage "
    "and median PBC-aware protein-relative ligand RMSD < 3 A"
)
COHORT_DEFINITION = (
    "same 52 complete trajectories with a locked late outcome and all "
    "5/10/15/20/30 ns corrected-RMSD checkpoint measurements"
)
CHECKPOINT_LOCK_ID = "same_trajectory_checkpoint_20ns_v1"
ANALYSIS_LOCK_ID = "same_trajectory_triage_a3_a5_v1_2026-07-15"


def auc_from_scores(labels: np.ndarray, scores: np.ndarray) -> float:
    """Compute AUROC by pairwise ordering without relying on row independence."""
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    positive = scores[labels == 1]
    negative = scores[labels == 0]
    if len(positive) == 0 or len(negative) == 0:
        return float("nan")
    difference = positive[:, None] - negative[None, :]
    return float(
        (np.sum(difference > 0) + 0.5 * np.sum(difference == 0))
        / difference.size
    )


def classification_metrics(labels: np.ndarray, decisions: np.ndarray) -> dict[str, float]:
    labels = np.asarray(labels, dtype=int)
    decisions = np.asarray(decisions, dtype=int)
    tn = int(np.sum((decisions == 0) & (labels == 0)))
    fn = int(np.sum((decisions == 0) & (labels == 1)))
    tp = int(np.sum((decisions == 1) & (labels == 1)))
    fp = int(np.sum((decisions == 1) & (labels == 0)))
    return {
        "balanced_accuracy": float(balanced_accuracy_score(labels, decisions)),
        "sensitivity_retained": float(recall_score(labels, decisions, pos_label=1)),
        "specificity_nonretained": float(recall_score(labels, decisions, pos_label=0)),
        "correct_stop_tn": tn,
        "false_stop_fn": fn,
        "correct_continue_tp": tp,
        "nonretained_continue_fp": fp,
    }


def cluster_bootstrap_interval(
    table: pd.DataFrame,
    value_column: str,
    metric: str,
    sign: float = 1.0,
    seed: int = SEED,
) -> tuple[float, float]:
    """Resample proteins and evaluate a fixed score or fixed decision.

    The model is intentionally not refit in a bootstrap draw. This matches the
    uncertainty scope of script 26 and is stated explicitly in every output.
    """
    rng = np.random.default_rng(seed)
    proteins = table[PROTEIN].dropna().unique()
    blocks = [
        (
            table.loc[table[PROTEIN].eq(protein), LABEL].to_numpy(int),
            table.loc[table[PROTEIN].eq(protein), value_column].to_numpy(),
        )
        for protein in proteins
    ]
    draws: list[float] = []
    for _ in range(N_BOOTSTRAP):
        sampled = rng.integers(0, len(blocks), size=len(blocks))
        labels = np.concatenate([blocks[index][0] for index in sampled])
        values = np.concatenate([blocks[index][1] for index in sampled])
        if len(np.unique(labels)) < 2:
            continue
        if metric == "auroc":
            estimate = auc_from_scores(labels, sign * values.astype(float))
        elif metric == "balanced_accuracy":
            estimate = balanced_accuracy_score(labels, values.astype(int))
        else:
            raise ValueError(f"unsupported bootstrap metric: {metric}")
        if np.isfinite(estimate):
            draws.append(float(estimate))
    if not draws:
        return float("nan"), float("nan")
    lo, hi = np.percentile(draws, [2.5, 97.5])
    return float(lo), float(hi)


def lopo_probabilities(
    table: pd.DataFrame, feature_columns: list[str]
) -> tuple[np.ndarray, np.ndarray, list[str], np.ndarray, np.ndarray]:
    """Outer leave-one-protein-out probabilities with training-only scaling."""
    probability = np.full(len(table), np.nan)
    decision = np.full(len(table), -1, dtype=int)
    folds = [""] * len(table)
    train_n = np.zeros(len(table), dtype=int)
    train_proteins = np.zeros(len(table), dtype=int)

    for held_out_protein in sorted(table[PROTEIN].unique()):
        test_mask = table[PROTEIN].eq(held_out_protein).to_numpy()
        train = table.loc[~test_mask]
        if train[LABEL].nunique() != 2:
            raise AssertionError(
                f"training fold for {held_out_protein} does not contain both classes"
            )

        scaler = StandardScaler(with_mean=True, with_std=True)
        x_train = scaler.fit_transform(train[feature_columns].to_numpy())
        model = LogisticRegression(**MODEL_PARAMETERS)
        model.fit(x_train, train[LABEL].astype(int).to_numpy())

        x_test = scaler.transform(table.loc[test_mask, feature_columns].to_numpy())
        fold_probability = model.predict_proba(x_test)[:, 1]
        probability[test_mask] = fold_probability
        decision[test_mask] = (fold_probability >= 0.5).astype(int)
        for index in np.flatnonzero(test_mask):
            folds[index] = f"held_out_protein={held_out_protein}"
        train_n[test_mask] = len(train)
        train_proteins[test_mask] = train[PROTEIN].nunique()

    if np.isnan(probability).any() or np.any(decision < 0):
        raise AssertionError("at least one outer-LOPO row did not receive a prediction")
    return probability, decision, folds, train_n, train_proteins


def locked_cohort(source: pd.DataFrame) -> pd.DataFrame:
    required = [
        "complex_id",
        PROTEIN,
        LABEL,
        "eligible_locked_primary_pose_model",
        "same_traj_pose_rmsd_18_20_A",
        "same_traj_centroid_displacement_18_20_A",
        "same_traj_pocket_retention_18_20",
        "separate_launch_pose_rmsd_18_20_A",
        "separate_launch_pocket_retention_18_20",
        "baseline_no_pbc_pose_rmsd_18_20_A",
        "baseline_ligand_self_aligned_rmsd_18_20_A",
        "baseline_corrected_final_frame_pose_rmsd_20ns_A",
        "locked_lopo_probability_retained",
        "locked_lopo_decision_continue",
    ] + [f"same_traj_pose_rmsd_{checkpoint}ns_A" for checkpoint in CHECKPOINTS_NS]
    missing = [column for column in required if column not in source.columns]
    if missing:
        raise ValueError(f"source-of-truth table is missing columns: {missing}")

    cohort = source.loc[source["eligible_locked_primary_pose_model"].eq(True)].copy()
    cohort = cohort.dropna(subset=required[2:]).reset_index(drop=True)
    if len(cohort) != 52:
        raise AssertionError(f"expected locked common cohort n=52, found {len(cohort)}")
    if cohort[PROTEIN].nunique() != 15:
        raise AssertionError("expected all 15 protein folds in the locked cohort")
    if int(cohort[LABEL].sum()) != 13:
        raise AssertionError("expected 13 late-retained trajectories in locked cohort")
    if int((cohort[LABEL] == 0).sum()) != 39:
        raise AssertionError("expected 39 late-non-retained trajectories in locked cohort")
    return cohort


def build_checkpoint_results(
    cohort: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    prediction_frames: list[pd.DataFrame] = []
    result_rows: list[dict] = []

    for position, checkpoint in enumerate(CHECKPOINTS_NS):
        feature = f"same_traj_pose_rmsd_{checkpoint}ns_A"
        probability, decision, folds, train_n, train_proteins = lopo_probabilities(
            cohort, [feature]
        )
        checkpoint_predictions = pd.DataFrame(
            {
                "complex_id": cohort["complex_id"],
                "protein": cohort[PROTEIN],
                "checkpoint_ns": checkpoint,
                "feature_window_ns": f"({checkpoint - 2},{checkpoint}]",
                "corrected_pose_rmsd_A": cohort[feature],
                "late_pose_retained_locked": cohort[LABEL].astype(int),
                "outer_lopo_fold": folds,
                "outer_train_n": train_n,
                "outer_train_proteins": train_proteins,
                "oof_probability_retained": probability,
                "decision_threshold": 0.5,
                "oof_decision_continue": decision,
            }
        )
        prediction_frames.append(checkpoint_predictions)

        temp = cohort[[PROTEIN, LABEL, feature]].copy()
        temp["oof_probability_retained"] = probability
        temp["oof_decision_continue"] = decision
        continuous_auc = auc_from_scores(temp[LABEL], -temp[feature])
        model_auc = auc_from_scores(temp[LABEL], temp["oof_probability_retained"])
        class_metrics = classification_metrics(temp[LABEL], decision)
        continuous_ci = cluster_bootstrap_interval(
            temp, feature, "auroc", sign=-1.0, seed=SEED
        )
        model_auc_ci = cluster_bootstrap_interval(
            temp,
            "oof_probability_retained",
            "auroc",
            sign=1.0,
            seed=SEED + 20 + position,
        )
        ba_ci = cluster_bootstrap_interval(
            temp,
            "oof_decision_continue",
            "balanced_accuracy",
            seed=SEED + 40 + position,
        )
        result_rows.append(
            {
                "checkpoint_ns": checkpoint,
                "checkpoint_role": (
                    "locked primary checkpoint for future prospective validation"
                    if checkpoint == 20
                    else "secondary retrospective checkpoint exploration"
                ),
                "feature": feature,
                "feature_definition": (
                    "mean PBC-aware protein-CA-aligned ligand heavy-atom RMSD over "
                    f"({checkpoint - 2},{checkpoint}] ns; lower predicts retention"
                ),
                "n": len(temp),
                "n_proteins": temp[PROTEIN].nunique(),
                "n_retained": int(temp[LABEL].sum()),
                "n_nonretained": int((temp[LABEL] == 0).sum()),
                "continuous_feature_auroc": continuous_auc,
                "continuous_feature_auroc_ci_lo": continuous_ci[0],
                "continuous_feature_auroc_ci_hi": continuous_ci[1],
                "oof_probability_auroc": model_auc,
                "oof_probability_auroc_ci_lo": model_auc_ci[0],
                "oof_probability_auroc_ci_hi": model_auc_ci[1],
                "oof_balanced_accuracy_p05": class_metrics["balanced_accuracy"],
                "oof_balanced_accuracy_ci_lo": ba_ci[0],
                "oof_balanced_accuracy_ci_hi": ba_ci[1],
                **{key: value for key, value in class_metrics.items() if key != "balanced_accuracy"},
                "validation_design": (
                    "outer leave-one-protein-out; training-fold-only StandardScaler; "
                    "fixed L2 logistic C=1; fixed p>=0.5 continue decision; no tuning"
                ),
                "uncertainty_method": (
                    f"{N_BOOTSTRAP}-draw protein-clustered bootstrap of fixed scores/"
                    "OOF predictions; models not refit per draw"
                ),
                "outcome_definition": OUTCOME_DEFINITION,
                "cohort_definition": COHORT_DEFINITION,
            }
        )

    predictions = pd.concat(prediction_frames, ignore_index=True)
    results = pd.DataFrame(result_rows)
    if len(predictions) != 52 * len(CHECKPOINTS_NS):
        raise AssertionError("checkpoint prediction table has an unexpected row count")
    return results, predictions


def metric_registry_row(
    cohort: pd.DataFrame,
    result_id: str,
    task: str,
    feature_definition: str,
    feature_column: str,
    window: str,
    sign: float,
    seed: int,
    interpretation: str,
) -> dict:
    estimate = auc_from_scores(cohort[LABEL], sign * cohort[feature_column])
    ci_lo, ci_hi = cluster_bootstrap_interval(
        cohort, feature_column, "auroc", sign=sign, seed=seed
    )
    return {
        "result_id": result_id,
        "result_status": "locked",
        "task": task,
        "feature_definition": feature_definition,
        "feature_column": feature_column,
        "feature_window": window,
        "outcome_definition": OUTCOME_DEFINITION,
        "cohort_definition": COHORT_DEFINITION,
        "n": len(cohort),
        "n_proteins": cohort[PROTEIN].nunique(),
        "n_retained": int(cohort[LABEL].sum()),
        "n_nonretained": int((cohort[LABEL] == 0).sum()),
        "validation_design": "fixed continuous score on the locked common cohort",
        "metric": "AUROC",
        "estimate": estimate,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "ci_level": 0.95,
        "ci_method": (
            f"protein-clustered bootstrap, {N_BOOTSTRAP} draws; fixed score; "
            "models not applicable"
        ),
        "decision_threshold": np.nan,
        "model_specification": "none; univariable continuous score",
        "interpretation": interpretation,
    }


def build_headline_registry(cohort: pd.DataFrame) -> pd.DataFrame:
    rows = [
        metric_registry_row(
            cohort,
            "same_corrected_pose_rmsd_20ns_auc",
            "same-trajectory discrimination",
            "Mean corrected ligand heavy-atom RMSD after protein-CA alignment and PBC minimum-image correction; lower predicts retention",
            "same_traj_pose_rmsd_18_20_A",
            "(18,20] ns",
            -1.0,
            SEED,
            "Primary physical forecasting signal; measures within-run persistence",
        ),
        metric_registry_row(
            cohort,
            "same_no_pbc_pose_rmsd_20ns_auc",
            "coordinate control",
            "Mean protein-aligned ligand RMSD without PBC correction; lower predicts retention",
            "baseline_no_pbc_pose_rmsd_18_20_A",
            "(18,20] ns",
            -1.0,
            SEED + 4,
            "PBC-naive control; periodic wrapping can corrupt protein-relative motion",
        ),
        metric_registry_row(
            cohort,
            "same_ligand_self_aligned_rmsd_20ns_auc",
            "coordinate control",
            "Mean ligand-self-aligned RMSD; lower predicts retention",
            "baseline_ligand_self_aligned_rmsd_18_20_A",
            "(18,20] ns",
            -1.0,
            SEED + 5,
            "Internal-conformation control that discards protein-relative translation/rotation",
        ),
        metric_registry_row(
            cohort,
            "same_corrected_final_frame_20ns_auc",
            "aggregation control",
            "Corrected ligand RMSD at the single saved frame nearest 20 ns; lower predicts retention",
            "baseline_corrected_final_frame_pose_rmsd_20ns_A",
            "single frame nearest 20 ns",
            -1.0,
            SEED + 7,
            "Distinct final-frame estimand; must not be merged into a range with the 2 ns mean",
        ),
        metric_registry_row(
            cohort,
            "same_pocket_retention_20ns_auc",
            "same-trajectory supporting measurement",
            "Mean fraction of baseline pocket neighbors retained; higher predicts retention",
            "same_traj_pocket_retention_18_20",
            "(18,20] ns",
            1.0,
            SEED + 1,
            "Supporting pocket-neighborhood measurement, not the primary model",
        ),
        metric_registry_row(
            cohort,
            "separate_corrected_pose_rmsd_20ns_auc",
            "separate-launch transfer",
            "Mean corrected RMSD from a separately launched trajectory; lower predicts the target launch's late retention",
            "separate_launch_pose_rmsd_18_20_A",
            "(18,20] ns",
            -1.0,
            SEED + 2,
            "Tests transfer to one separately launched realization, not same-run forecasting",
        ),
        metric_registry_row(
            cohort,
            "separate_pocket_retention_20ns_auc",
            "separate-launch transfer",
            "Mean pocket-neighborhood retention from a separately launched trajectory; higher predicts the target launch's late retention",
            "separate_launch_pocket_retention_18_20",
            "(18,20] ns",
            1.0,
            SEED + 3,
            "Supporting separate-launch comparison; point estimate is near chance",
        ),
    ]

    model_temp = cohort[[PROTEIN, LABEL, "locked_lopo_probability_retained", "locked_lopo_decision_continue"]].copy()
    model_auc = auc_from_scores(
        model_temp[LABEL], model_temp["locked_lopo_probability_retained"]
    )
    model_auc_ci = cluster_bootstrap_interval(
        model_temp,
        "locked_lopo_probability_retained",
        "auroc",
        seed=SEED + 6,
    )
    model_ba = classification_metrics(
        model_temp[LABEL], model_temp["locked_lopo_decision_continue"]
    )["balanced_accuracy"]
    model_ba_ci = cluster_bootstrap_interval(
        model_temp,
        "locked_lopo_decision_continue",
        "balanced_accuracy",
        seed=SEED + 7,
    )
    common_model = {
        "result_status": "locked",
        "task": "same-trajectory outer-LOPO model",
        "feature_definition": (
            "Mean corrected ligand RMSD plus corrected ligand-centroid displacement; "
            "both PBC-aware and protein-CA-aligned"
        ),
        "feature_column": (
            "same_traj_pose_rmsd_18_20_A;"
            "same_traj_centroid_displacement_18_20_A"
        ),
        "feature_window": "(18,20] ns",
        "outcome_definition": OUTCOME_DEFINITION,
        "cohort_definition": COHORT_DEFINITION,
        "n": len(cohort),
        "n_proteins": cohort[PROTEIN].nunique(),
        "n_retained": int(cohort[LABEL].sum()),
        "n_nonretained": int((cohort[LABEL] == 0).sum()),
        "validation_design": (
            "outer leave-one-protein-out; training-fold-only StandardScaler; "
            "no tuning"
        ),
        "ci_level": 0.95,
        "ci_method": (
            f"protein-clustered bootstrap, {N_BOOTSTRAP} draws; fixed outer-LOPO "
            "predictions; models not refit per draw"
        ),
        "decision_threshold": 0.5,
        "model_specification": (
            "L2 logistic regression; C=1; class_weight=balanced; lbfgs; "
            "max_iter=100; tol=1e-4; predict_proba output is uncalibrated and "
            "is interpreted as a model score"
        ),
        "interpretation": (
            "Evaluation model only: 15 fold-specific models; no final deployment "
            "model has been trained or prospectively validated"
        ),
    }
    rows.extend(
        [
            {
                "result_id": "same_two_feature_lopo_oof_probability_auc",
                **common_model,
                "metric": "AUROC",
                "estimate": model_auc,
                "ci_lo": model_auc_ci[0],
                "ci_hi": model_auc_ci[1],
            },
            {
                "result_id": "same_two_feature_lopo_balanced_accuracy_p05",
                **common_model,
                "metric": "balanced accuracy",
                "estimate": model_ba,
                "ci_lo": model_ba_ci[0],
                "ci_hi": model_ba_ci[1],
            },
        ]
    )
    return pd.DataFrame(rows)


def read_legacy_artifacts() -> dict[str, float]:
    for path in [
        LEGACY_MODEL_PREDICTIONS,
        LEGACY_CHECKPOINT_RESULTS,
        LEGACY_COORDINATE_RESULTS,
        LEGACY_REPLICA_RESULTS,
        LEGACY_VALIDATION2_USABLE,
    ]:
        if not path.exists():
            raise FileNotFoundError(f"missing legacy provenance artifact: {path}")

    model = pd.read_csv(LEGACY_MODEL_PREDICTIONS)
    checkpoint = pd.read_csv(LEGACY_CHECKPOINT_RESULTS)
    coordinate = pd.read_csv(LEGACY_COORDINATE_RESULTS)
    replica = pd.read_csv(LEGACY_REPLICA_RESULTS)
    validation2 = pd.read_csv(LEGACY_VALIDATION2_USABLE)
    checkpoint20 = checkpoint.loc[np.isclose(checkpoint["checkpoint_ns"], 20.0)].iloc[0]
    corrected_window = coordinate.loc[
        coordinate["method"].str.startswith("validated (protein-CA")
    ].iloc[0]
    corrected_final = coordinate.loc[
        coordinate["method"].str.startswith("final-frame-only")
    ].iloc[0]
    no_pbc = coordinate.loc[coordinate["method"].str.contains("NO PBC")].iloc[0]
    ligand_self = coordinate.loc[
        coordinate["method"].str.contains("ligand-self")
    ].iloc[0]
    retention_replica = replica.loc[
        replica["metric"].str.startswith("Pocket retention")
    ].iloc[0]
    pose_replica = replica.loc[
        replica["metric"].str.startswith("Pose perturbation")
    ].iloc[0]
    return {
        "legacy_two_feature_ba": float(
            balanced_accuracy_score(model["y_stable100"], model["pred_pose_only"])
        ),
        "legacy_checkpoint20_rmsd_only_ba": float(checkpoint20["lopo_balanced_accuracy"]),
        "legacy_corrected_window_auc": float(corrected_window["auroc"]),
        "legacy_corrected_final_auc": float(corrected_final["auroc"]),
        "legacy_no_pbc_auc": float(no_pbc["auroc"]),
        "legacy_ligand_self_auc": float(ligand_self["auroc"]),
        "legacy_separate_retention_auc": float(
            retention_replica["auroc_independent_replica"]
        ),
        "legacy_separate_pose_auc": float(pose_replica["auroc_independent_replica"]),
        "legacy_subset_retention_auc": auc_from_scores(
            validation2["y_stable100"], validation2["r_pocket"]
        ),
        "legacy_subset_retention_n": len(validation2),
        "legacy_subset_retained_n": int(validation2["y_stable100"].sum()),
        "legacy_subset_nonretained_n": int(
            (validation2["y_stable100"] == 0).sum()
        ),
    }


def build_reconciliation(
    registry: pd.DataFrame, legacy: dict[str, float]
) -> pd.DataFrame:
    estimates = registry.set_index("result_id")["estimate"]
    return pd.DataFrame(
        [
            {
                "conflict_id": "balanced_accuracy_0.778_vs_0.783",
                "legacy_values": (
                    f"{legacy['legacy_two_feature_ba']:.9f} vs "
                    f"{legacy['legacy_checkpoint20_rmsd_only_ba']:.9f}"
                ),
                "legacy_provenance": (
                    "script 22 two-feature RMSD+centroid LOPO on n=58 vs script 24 "
                    "one-feature RMSD-only LOPO on n=58"
                ),
                "resolution": (
                    "They are different models, not rounding variants. The locked primary "
                    "two-feature model has balanced accuracy "
                    f"{estimates['same_two_feature_lopo_balanced_accuracy_p05']:.9f} "
                    "on n=52. Checkpoint-specific one-feature values remain in "
                    "checkpoint_model_results.csv."
                ),
                "locked_replacement_result_id": "same_two_feature_lopo_balanced_accuracy_p05",
            },
            {
                "conflict_id": "corrected_pose_auc_0.877_to_0.886_range",
                "legacy_values": (
                    f"{legacy['legacy_corrected_window_auc']:.9f} to "
                    f"{legacy['legacy_corrected_final_auc']:.9f}"
                ),
                "legacy_provenance": (
                    "mean over (18,20] ns on n=58 vs a single frame nearest 20 ns "
                    "on n=59"
                ),
                "resolution": (
                    "The range mixed different aggregation rules and cohorts. Use exactly "
                    f"{estimates['same_corrected_pose_rmsd_20ns_auc']:.9f} for the "
                    "locked (18,20] mean-RMSD signal on n=52; report the final-frame "
                    "control separately."
                ),
                "locked_replacement_result_id": "same_corrected_pose_rmsd_20ns_auc",
            },
            {
                "conflict_id": "coordinate_controls_shared_0.549_to_0.585_range",
                "legacy_values": (
                    f"no-PBC={legacy['legacy_no_pbc_auc']:.9f}; "
                    f"ligand-self={legacy['legacy_ligand_self_auc']:.9f}"
                ),
                "legacy_provenance": (
                    "two distinct measurements were collapsed into one shared range on the "
                    "legacy n=58 cohort"
                ),
                "resolution": (
                    "Report exact locked n=52 values separately: no-PBC="
                    f"{estimates['same_no_pbc_pose_rmsd_20ns_auc']:.9f}; "
                    "ligand-self="
                    f"{estimates['same_ligand_self_aligned_rmsd_20ns_auc']:.9f}."
                ),
                "locked_replacement_result_id": (
                    "same_no_pbc_pose_rmsd_20ns_auc;"
                    "same_ligand_self_aligned_rmsd_20ns_auc"
                ),
            },
            {
                "conflict_id": "separate_launch_pocket_retention_0.54_vs_0.603",
                "legacy_values": (
                    f"current controlled artifact={legacy['legacy_separate_retention_auc']:.9f}; "
                    f"selective P20-perturbed subset={legacy['legacy_subset_retention_auc']:.9f}"
                ),
                "legacy_provenance": (
                    "0.602978 is direct r_pocket AUROC in validation2_100ns_outcomes_usable.csv "
                    f"on a selective n={legacy['legacy_subset_retention_n']} P20-perturbed "
                    f"subset ({legacy['legacy_subset_retained_n']}/"
                    f"{legacy['legacy_subset_nonretained_n']} classes); 0.538020 is the "
                    "script-21 paired same-versus-separate n=58 intersection"
                ),
                "resolution": (
                    "The estimates answer different cohort questions. Do not use the "
                    "selective subset estimate as full-cohort transfer performance. The "
                    "exact locked identical-cohort estimate is "
                    f"{estimates['separate_pocket_retention_20ns_auc']:.9f} on n=52."
                ),
                "locked_replacement_result_id": "separate_pocket_retention_20ns_auc",
            },
            {
                "conflict_id": "separate_launch_pose_approximate_range",
                "legacy_values": f"current controlled artifact={legacy['legacy_separate_pose_auc']:.9f}",
                "legacy_provenance": "script-21 controlled comparison on legacy n=58 cohort",
                "resolution": (
                    "Use the exact locked identical-cohort estimate "
                    f"{estimates['separate_corrected_pose_rmsd_20ns_auc']:.9f} on n=52."
                ),
                "locked_replacement_result_id": "separate_corrected_pose_rmsd_20ns_auc",
            },
        ]
    )


def build_validation_audit() -> pd.DataFrame:
    rows = [
        (
            "outer_fold_unit",
            "pass",
            "All rows from one protein are held out together in each of 15 outer folds.",
            "25_build_same_trajectory_source_of_truth.py: recompute_lopo_pose_model; this script: lopo_probabilities",
            "Evaluation is protein-held-out within this 15-protein panel, not external validation.",
        ),
        (
            "scaling_scope",
            "pass",
            "StandardScaler is fit only on outer-training proteins and applied unchanged to the held-out protein.",
            "scaler.fit_transform(train); scaler.transform(test)",
            "No preprocessing leakage was found for the locked models.",
        ),
        (
            "regularization",
            "documented_fixed_default",
            "L2 logistic regression with C=1.0, class_weight=balanced, lbfgs, max_iter=100, tol=1e-4.",
            "MODEL_PARAMETERS and source-of-truth model_audit_metrics.csv",
            "C was inherited/fixed, not selected by evidence from an inner loop.",
        ),
        (
            "hyperparameter_tuning",
            "none",
            "No model hyperparameters are tuned; therefore no inner tuning loop is needed for the reported fixed model.",
            "No GridSearchCV, RandomizedSearchCV, or training-fold selection in scripts 22/24/25/27.",
            "Future tuning must occur entirely inside an inner training-only loop.",
        ),
        (
            "decision_threshold",
            "fixed_retrospective_operating_point",
            "Continue when outer-LOPO probability is >=0.5; stop otherwise.",
            "25_build_same_trajectory_source_of_truth.py and 26_plot_locked_claim_results.py",
            "The threshold was not selected to guarantee a false-stop ceiling; do not call it false-stop-controlled.",
        ),
        (
            "policy_prediction_source",
            "pass",
            "Every locked p=0.5 policy decision is derived from an outer-LOPO probability for that row.",
            "locked_lopo_fold and locked_lopo_probability_retained are nonmissing for all 52 policy rows.",
            "The threshold sweep remains retrospective and does not constitute prospective control.",
        ),
        (
            "feature_and_checkpoint_selection",
            "exploratory_not_nested",
            "The 20 ns checkpoint was historically motivated by the existing short-MD screen, but coordinate controls and 5/10/15/30 ns comparisons were examined retrospectively after outcome access.",
            "scripts 22 then 24; repository timestamps support ordering, but no preregistration or committed history proves prospective specification.",
            "Do not claim that 20 ns or the feature set was selected by nested validation. Treat alternatives as secondary exploration.",
        ),
        (
            "score_calibration",
            "not_calibrated",
            "The class-weighted logistic predict_proba output has no Platt, isotonic, or other calibration step and is treated as a model score.",
            "No calibration estimator is present in scripts 22/24/25/27.",
            "Do not interpret the score as an absolute probability of late retention.",
        ),
        (
            "evaluation_vs_deployment_model",
            "deployment_model_absent",
            "LOPO produces 15 fold-specific evaluation models. No final all-development-protein deployment model is generated here.",
            "One held-out-protein model per fold; no serialized final model artifact.",
            "After the procedure is frozen, train a separate final candidate and test it in a prospective shadow study.",
        ),
        (
            "bootstrap_scope",
            "limited_but_explicit",
            f"Intervals use {N_BOOTSTRAP} protein-clustered draws of fixed continuous scores or saved OOF predictions; models are not refit per draw.",
            "cluster_bootstrap_interval and script 26",
            "Intervals quantify clustered sampling variation conditional on the fitted OOF scores, not full pipeline-selection uncertainty.",
        ),
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "audit_item",
            "status",
            "implementation",
            "evidence",
            "required_interpretation_or_action",
        ],
    )


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def build_primary_fold_audit(cohort: pd.DataFrame) -> pd.DataFrame:
    """Persist scaler/model details for every primary two-feature LOPO fold."""
    features = [
        "same_traj_pose_rmsd_18_20_A",
        "same_traj_centroid_displacement_18_20_A",
    ]
    rows: list[dict] = []
    regenerated = np.full(len(cohort), np.nan)
    for held_out_protein in sorted(cohort[PROTEIN].unique()):
        test_mask = cohort[PROTEIN].eq(held_out_protein).to_numpy()
        train = cohort.loc[~test_mask]
        test = cohort.loc[test_mask]
        scaler = StandardScaler(with_mean=True, with_std=True)
        x_train = scaler.fit_transform(train[features].to_numpy())
        model = LogisticRegression(**MODEL_PARAMETERS)
        model.fit(x_train, train[LABEL].astype(int).to_numpy())
        regenerated[test_mask] = model.predict_proba(
            scaler.transform(test[features].to_numpy())
        )[:, 1]
        train_labels = train[LABEL].astype(int)
        class_weight_0 = len(train) / (2.0 * int((train_labels == 0).sum()))
        class_weight_1 = len(train) / (2.0 * int((train_labels == 1).sum()))
        n_iter = int(model.n_iter_[0])
        rows.append(
            {
                "analysis_id": ANALYSIS_LOCK_ID,
                "fold_id": f"held_out_protein={held_out_protein}",
                "held_out_protein": held_out_protein,
                "train_n": len(train),
                "test_n": len(test),
                "train_proteins": train[PROTEIN].nunique(),
                "train_retained": int(train_labels.sum()),
                "train_nonretained": int((train_labels == 0).sum()),
                "test_retained": int(test[LABEL].sum()),
                "test_nonretained": int((test[LABEL] == 0).sum()),
                "scaler_mean_rmsd": float(scaler.mean_[0]),
                "scaler_scale_rmsd": float(scaler.scale_[0]),
                "scaler_mean_centroid": float(scaler.mean_[1]),
                "scaler_scale_centroid": float(scaler.scale_[1]),
                "effective_class_weight_0": class_weight_0,
                "effective_class_weight_1": class_weight_1,
                "coef_standardized_rmsd": float(model.coef_[0, 0]),
                "coef_standardized_centroid": float(model.coef_[0, 1]),
                "intercept": float(model.intercept_[0]),
                "n_iter": n_iter,
                "max_iter": MODEL_PARAMETERS["max_iter"],
                "convergence_status": (
                    "did_not_hit_iteration_limit"
                    if n_iter < MODEL_PARAMETERS["max_iter"]
                    else "hit_iteration_limit_review_warning_log"
                ),
                "test_complex_ids_sha256": sha256_text(
                    "\n".join(sorted(test["complex_id"].astype(str)))
                ),
            }
        )
    if not np.allclose(
        regenerated,
        cohort["locked_lopo_probability_retained"].to_numpy(float),
        atol=1e-12,
        rtol=1e-12,
    ):
        raise AssertionError("fold audit refit does not reproduce saved locked OOF scores")
    return pd.DataFrame(rows)


def build_model_specification(cohort: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "analysis_id": ANALYSIS_LOCK_ID,
                "analysis_role": "locked retrospective evaluation",
                "outcome_id": "late_pose_retained_rmsd_v1",
                "cohort_id": "locked_complete_70_100ns_common_checkpoint_n52",
                "n": len(cohort),
                "n_proteins": cohort[PROTEIN].nunique(),
                "n_retained": int(cohort[LABEL].sum()),
                "n_nonretained": int((cohort[LABEL] == 0).sum()),
                "checkpoint_ns": 20,
                "feature_window": "(18,20] ns",
                "features": (
                    "same_traj_pose_rmsd_18_20_A;"
                    "same_traj_centroid_displacement_18_20_A"
                ),
                "outer_split_unit": "protein",
                "outer_folds": cohort[PROTEIN].nunique(),
                "scaler": "StandardScaler(with_mean=True,with_std=True)",
                "scaler_fit_scope": "outer_training_fold_only",
                "estimator": "LogisticRegression",
                "penalty": MODEL_PARAMETERS["penalty"],
                "C": MODEL_PARAMETERS["C"],
                "solver": MODEL_PARAMETERS["solver"],
                "max_iter": MODEL_PARAMETERS["max_iter"],
                "tol": MODEL_PARAMETERS["tol"],
                "class_weight": MODEL_PARAMETERS["class_weight"],
                "hyperparameter_selection": "none_inherited_fixed_default",
                "inner_cv": "none",
                "threshold": 0.5,
                "threshold_selection": "fixed_default_not_safety_tuned",
                "score_calibrated": False,
                "score_interpretation": "uncalibrated ranking/decision score",
                "prediction_provenance": "outer_LOPO_only",
                "feature_selection_timing": "retrospective_not_nested",
                "checkpoint_selection_timing": (
                    "20ns_inherited_before_exploratory_scan_formal_prespecification_unverified"
                ),
                "ci_method": (
                    f"{N_BOOTSTRAP}-draw protein-clustered bootstrap of fixed scores"
                ),
                "bootstrap_refit": False,
                "deployment_model_status": "not_trained",
                "prospective_validation": False,
                "scikit_learn_version": sklearn.__version__,
                "source_script": Path(__file__).name,
            }
        ]
    )


def build_analysis_lock(cohort: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "analysis_lock_id": ANALYSIS_LOCK_ID,
                "outcome_definition_id": "late_pose_retained_rmsd_v1",
                "checkpoint_lock_id": CHECKPOINT_LOCK_ID,
                "checkpoint_ns": 20,
                "checkpoint_window": "(18,20] ns",
                "checkpoint_rationale": (
                    "inherited operational checkpoint from the existing 20 ns screening "
                    "workflow and the first same-trajectory analysis; not chosen because "
                    "it maximized the retrospective checkpoint curve"
                ),
                "historical_prespecification_status": (
                    "historically motivated before the multi-checkpoint curve, but not "
                    "prospectively preregistered before outcome inspection"
                ),
                "alternative_checkpoint_status": (
                    "5/10/15/30 ns are secondary retrospective exploration"
                ),
                "locked_features": (
                    "same_traj_pose_rmsd_18_20_A;"
                    "same_traj_centroid_displacement_18_20_A"
                ),
                "model": (
                    "training-only StandardScaler + L2 logistic regression, C=1, "
                    "class_weight=balanced"
                ),
                "decision_rule": "continue if probability_retained >= 0.5; stop otherwise",
                "validation": "outer leave-one-protein-out, 15 protein folds",
                "cohort_n": len(cohort),
                "cohort_proteins": cohort[PROTEIN].nunique(),
                "cohort_retained": int(cohort[LABEL].sum()),
                "cohort_nonretained": int((cohort[LABEL] == 0).sum()),
                "deployment_status": (
                    "frozen for future prospective shadow validation; no deployable model "
                    "or prospective performance claim"
                ),
            }
        ]
    )


def validate_against_existing_outputs(
    registry: pd.DataFrame, checkpoint_results: pd.DataFrame
) -> None:
    if CLAIM_METRICS_TABLE.exists():
        claim = pd.read_csv(CLAIM_METRICS_TABLE).set_index("metric_id")
        registry_estimate = registry.set_index("result_id")["estimate"]
        mappings = {
            "same_corrected_pose_rmsd_20ns_auc": "same_corrected_pose",
            "same_no_pbc_pose_rmsd_20ns_auc": "same_no_pbc_pose",
            "same_ligand_self_aligned_rmsd_20ns_auc": "same_ligand_self_rmsd",
            "same_pocket_retention_20ns_auc": "same_pocket_retention",
            "separate_corrected_pose_rmsd_20ns_auc": "separate_corrected_pose",
            "separate_pocket_retention_20ns_auc": "separate_pocket_retention",
            "same_two_feature_lopo_oof_probability_auc": "locked_oof_probability",
        }
        for registry_id, claim_id in mappings.items():
            if not np.isclose(
                registry_estimate[registry_id], claim.loc[claim_id, "auroc"], atol=1e-12
            ):
                raise AssertionError(
                    f"locked result mismatch between scripts 26 and 27: {registry_id}"
                )

    if CLAIM_CHECKPOINT_TABLE.exists():
        claim_checkpoint = pd.read_csv(CLAIM_CHECKPOINT_TABLE).set_index("checkpoint_ns")
        current = checkpoint_results.set_index("checkpoint_ns")
        if not np.allclose(
            current.loc[list(CHECKPOINTS_NS), "continuous_feature_auroc"],
            claim_checkpoint.loc[list(CHECKPOINTS_NS), "auroc"],
            atol=1e-12,
        ):
            raise AssertionError("checkpoint AUROCs disagree between scripts 26 and 27")


def write_results_lock(
    cohort: pd.DataFrame,
    registry: pd.DataFrame,
    checkpoint: pd.DataFrame,
    reconciliation: pd.DataFrame,
    validation: pd.DataFrame,
) -> None:
    values = registry.set_index("result_id")
    cp = checkpoint.set_index("checkpoint_ns")
    policy = classification_metrics(
        cohort[LABEL], cohort["locked_lopo_decision_continue"]
    )
    n_stopped = policy["correct_stop_tn"] + policy["false_stop_fn"]
    compute_saved_fraction = n_stopped * 80.0 / (len(cohort) * 100.0)
    text = f"""# Results lock for manuscript TODO A3-A5

Generated by `27_audit_headlines_validation_checkpoints.py` from the locked
`same_trajectory_source_of_truth.csv`. Do not hand-edit generated CSV files.

## Locked cohort and outcome

- Cohort: 52 trajectories from 15 proteins; 13 late retained and 39 late
  non-retained.
- Outcome: `{OUTCOME_DEFINITION}`.
- Every checkpoint comparison uses these same 52 rows. The former expected
  common cohort of 57 was based on the legacy outcome and is obsolete after A1.

## A3: exact headline values

- Corrected same-trajectory RMSD, mean over `(18,20]` ns: AUROC
  {values.loc['same_corrected_pose_rmsd_20ns_auc', 'estimate']:.6f}, 95% CI
  [{values.loc['same_corrected_pose_rmsd_20ns_auc', 'ci_lo']:.6f},
  {values.loc['same_corrected_pose_rmsd_20ns_auc', 'ci_hi']:.6f}].
- Two-feature outer-LOPO model: OOF-probability AUROC
  {values.loc['same_two_feature_lopo_oof_probability_auc', 'estimate']:.6f}, 95% CI
  [{values.loc['same_two_feature_lopo_oof_probability_auc', 'ci_lo']:.6f},
  {values.loc['same_two_feature_lopo_oof_probability_auc', 'ci_hi']:.6f}]; balanced
  accuracy at p=0.5
  {values.loc['same_two_feature_lopo_balanced_accuracy_p05', 'estimate']:.6f}, 95% CI
  [{values.loc['same_two_feature_lopo_balanced_accuracy_p05', 'ci_lo']:.6f},
  {values.loc['same_two_feature_lopo_balanced_accuracy_p05', 'ci_hi']:.6f}].
- Exact coordinate controls: no-PBC AUROC
  {values.loc['same_no_pbc_pose_rmsd_20ns_auc', 'estimate']:.6f}; ligand-self-aligned
  AUROC {values.loc['same_ligand_self_aligned_rmsd_20ns_auc', 'estimate']:.6f}.
- Exact separate-launch values: corrected pose AUROC
  {values.loc['separate_corrected_pose_rmsd_20ns_auc', 'estimate']:.6f}; pocket
  retention AUROC {values.loc['separate_pocket_retention_20ns_auc', 'estimate']:.6f}.

The stored `predict_proba` output is uncalibrated because the model uses balanced
class weights and no calibration step. It is therefore interpreted as a model
score, not an absolute probability of retention. The raw corrected-RMSD AUROC
and the two-feature model-score AUROC answer different
questions. The former measures discrimination by one physical score; the latter
is based on outer-LOPO probabilities from a fitted two-feature model. They must
be labeled separately.

`legacy_headline_reconciliation.csv` records why `0.778` and `0.783` were
different models, why `0.877-0.886` mixed window and final-frame estimands, why
the coordinate controls require separate exact rows, and why the older `0.603`
estimate came from a selective 44-row P20-perturbed subset and is
superseded for the full identical-cohort comparison.

## A4: validation and threshold lock

- All preprocessing is fit only on outer-training proteins.
- Model: L2 logistic regression, C=1.0, class weight `balanced`, `lbfgs`.
- Hyperparameter tuning: none.
- Decision threshold: fixed p=0.5. It was not selected to guarantee a
  false-stop ceiling.
- Policy calculations use only outer-LOPO probabilities, but the threshold
  sweep and operating-point assessment are retrospective.
- LOPO creates 15 evaluation models. No final all-development-protein deployment
  model is trained or claimed here.
- Feature-set and alternative-checkpoint comparisons are retrospective and not
  nested. The clustered intervals condition on fixed scores/OOF predictions and
  do not refit the entire modeling/selection procedure per draw.

At the fixed p=0.5 retrospective operating point, the exact confusion counts are
TN={policy['correct_stop_tn']}, FN={policy['false_stop_fn']},
TP={policy['correct_continue_tp']}, and FP={policy['nonretained_continue_fp']}.
Thus {n_stopped} of {len(cohort)} trajectories would stop at 20 ns, saving
{100 * compute_saved_fraction:.1f}% of the nominal simulated time; this is an
out-of-fold development-cohort estimate, not prospective false-stop control.

See `validation_procedure_audit.csv` for the itemized audit.
The exact one-row specification is in `primary_model_validation_specification.csv`,
and `lopo_fold_model_audit.csv` records every scaler, coefficient, intercept,
iteration count, and held-out-protein fold.

## A5: role of 20 ns

Twenty nanoseconds was inherited from the pre-existing short-MD screening
workflow (version-controlled evidence reaches back to commits `66a17805` and
`9f9df8cb`) and was the checkpoint used by the first same-trajectory analysis
before the multi-checkpoint curve was generated. However, the repository has no
preregistration or committed history proving that it was prospectively specified
before any outcome inspection. The defensible description is therefore:

> The 20 ns checkpoint was historically and operationally motivated, not selected
> because it maximized the retrospective checkpoint curve. Analyses at 5, 10, 15,
> and 30 ns are secondary retrospective exploration.

The frozen prospective checkpoint is now `{CHECKPOINT_LOCK_ID}`: use the mean
corrected features over `(18,20]` ns and do not change the checkpoint after
examining a future validation batch.

On the common 52-row cohort, continuous corrected-RMSD AUROCs are:

| Checkpoint | AUROC | 95% clustered CI | One-feature LOPO BA at p=0.5 |
|---:|---:|---:|---:|
| 5 ns | {cp.loc[5, 'continuous_feature_auroc']:.3f} | [{cp.loc[5, 'continuous_feature_auroc_ci_lo']:.3f}, {cp.loc[5, 'continuous_feature_auroc_ci_hi']:.3f}] | {cp.loc[5, 'oof_balanced_accuracy_p05']:.3f} |
| 10 ns | {cp.loc[10, 'continuous_feature_auroc']:.3f} | [{cp.loc[10, 'continuous_feature_auroc_ci_lo']:.3f}, {cp.loc[10, 'continuous_feature_auroc_ci_hi']:.3f}] | {cp.loc[10, 'oof_balanced_accuracy_p05']:.3f} |
| 15 ns | {cp.loc[15, 'continuous_feature_auroc']:.3f} | [{cp.loc[15, 'continuous_feature_auroc_ci_lo']:.3f}, {cp.loc[15, 'continuous_feature_auroc_ci_hi']:.3f}] | {cp.loc[15, 'oof_balanced_accuracy_p05']:.3f} |
| 20 ns | {cp.loc[20, 'continuous_feature_auroc']:.3f} | [{cp.loc[20, 'continuous_feature_auroc_ci_lo']:.3f}, {cp.loc[20, 'continuous_feature_auroc_ci_hi']:.3f}] | {cp.loc[20, 'oof_balanced_accuracy_p05']:.3f} |
| 30 ns | {cp.loc[30, 'continuous_feature_auroc']:.3f} | [{cp.loc[30, 'continuous_feature_auroc_ci_lo']:.3f}, {cp.loc[30, 'continuous_feature_auroc_ci_hi']:.3f}] | {cp.loc[30, 'oof_balanced_accuracy_p05']:.3f} |

The 10 ns point estimate is numerically highest, which is additional evidence
that 20 ns should be justified by operational history and frozen protocol rather
than by retrospective performance maximization.

## Generated audit files

- `headline_result_registry.csv`: manuscript-ready definitions and exact values.
- `legacy_headline_reconciliation.csv`: one row per legacy discrepancy.
- `validation_procedure_audit.csv`: preprocessing/model/threshold/selection audit.
- `primary_model_validation_specification.csv`: exact locked model specification.
- `lopo_fold_model_audit.csv`: one row per held-out-protein evaluation model.
- `checkpoint_model_results.csv`: common-cohort continuous and one-feature LOPO
  results at every checkpoint.
- `checkpoint_oof_predictions.csv`: all 260 row-checkpoint OOF predictions.
- `analysis_lock_a3_a5.csv`: the prospective 20 ns lock and deployment status.
"""
    (RESULT_DIR / "RESULTS_LOCK_A3_A5.md").write_text(text)


def main() -> None:
    if not SOURCE_TABLE.exists():
        raise FileNotFoundError(
            f"missing {SOURCE_TABLE}; run 25_build_same_trajectory_source_of_truth.py first"
        )
    source = pd.read_csv(SOURCE_TABLE)
    cohort = locked_cohort(source)

    checkpoint_results, checkpoint_predictions = build_checkpoint_results(cohort)
    registry = build_headline_registry(cohort)
    legacy = read_legacy_artifacts()
    reconciliation = build_reconciliation(registry, legacy)
    validation = build_validation_audit()
    model_specification = build_model_specification(cohort)
    fold_audit = build_primary_fold_audit(cohort)
    lock = build_analysis_lock(cohort)

    validate_against_existing_outputs(registry, checkpoint_results)

    registry.to_csv(RESULT_DIR / "headline_result_registry.csv", index=False)
    reconciliation.to_csv(
        RESULT_DIR / "legacy_headline_reconciliation.csv", index=False
    )
    validation.to_csv(RESULT_DIR / "validation_procedure_audit.csv", index=False)
    model_specification.to_csv(
        RESULT_DIR / "primary_model_validation_specification.csv", index=False
    )
    fold_audit.to_csv(RESULT_DIR / "lopo_fold_model_audit.csv", index=False)
    checkpoint_results.to_csv(
        RESULT_DIR / "checkpoint_model_results.csv", index=False
    )
    checkpoint_predictions.to_csv(
        RESULT_DIR / "checkpoint_oof_predictions.csv", index=False
    )
    lock.to_csv(RESULT_DIR / "analysis_lock_a3_a5.csv", index=False)
    write_results_lock(
        cohort, registry, checkpoint_results, reconciliation, validation
    )

    print(f"Wrote A3-A5 audit outputs to {RESULT_DIR}")
    print(
        registry[["result_id", "metric", "estimate", "ci_lo", "ci_hi"]]
        .to_string(index=False)
    )
    print("\nCommon-cohort checkpoint results:")
    print(
        checkpoint_results[
            [
                "checkpoint_ns",
                "continuous_feature_auroc",
                "oof_probability_auroc",
                "oof_balanced_accuracy_p05",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
