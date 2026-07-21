#!/usr/bin/env python3
"""PROTOTYPE: frozen leave-one-protein-out [T, Theta] joint model.

Targeted test of one physical hypothesis: that a learned rigid-body
combination of translation (T) and reorientation (Theta) -- with internal
deformation (I) excluded, since script 33 showed it is chance-level --
predicts late pose retention at least as well as corrected RMSD, using a
combination method that (unlike the naive per-trajectory-MAD-normalized
Z_pose in script 33) does not fail.

Frozen model spec (do not tune):
  p_hat = sigmoid(b0 + b_T * z(T) + b_Theta * z(Theta))
  - features: T and Theta only (I excluded from the primary model; a
    T+Theta+I variant is included ONLY as a negative control to confirm I
    adds nothing on top of T+Theta, not as a candidate to adopt)
  - L2 logistic regression, C=1.0 (sklearn default), class_weight="balanced"
  - no hyperparameter search
  - outer leave-one-protein-out (LOPO): every fold's StandardScaler and
    LogisticRegression are fit on the training proteins only, then applied
    to the held-out protein -- never fit on the full cohort. This differs
    from script 33's Z_pose, which normalized by each trajectory's OWN
    baseline MAD (a per-trajectory, label-blind but fold-blind scale); here
    the scale is a single training-fold statistic, applied uniformly.
  - checkpoints: 20 ns (primary, inherited from the locked pipeline),
    5 ns and 10 ns (exploratory -- the 10 ns rotation result was observed
    after inspecting the checkpoint curve, so it is reported as
    exploratory, not confirmatory)

Six models compared at each checkpoint, all evaluated through the SAME
frozen LOPO logistic procedure (including the single-feature ones) so
balanced accuracy / compute-saved use a consistently fold-calibrated
decision threshold rather than an arbitrary raw-value cutoff:
  1. corrected RMSD only         -- current benchmark
  2. T only                      -- translation baseline
  3. Theta only                  -- rotation baseline
  4. T + Theta                   -- proposed learned rigid-body combination
  5. T + Theta + I                -- negative control (I should add nothing)
  6. Z_pose (from script 33)     -- known-failed hand-built combination,
                                     included for a like-for-like comparison
                                     against a properly fit model

For each model/checkpoint: AUROC (protein-clustered bootstrap CI), paired
AUROC difference vs RMSD, balanced accuracy at p=0.5 (protein-clustered
bootstrap CI), max compute-saved fraction subject to false-stop rate
<= 10% (threshold swept independently per model), AUROC restricted to the
near-boundary subset (same definition as the locked pipeline's Table
non-obvious: |RMSD_checkpoint - 3.0 A| <= 1.5 A), and for the T+Theta
model, per-fold standardized-coefficient stability across the 16 LOPO
folds.

Usage:
    python 34_prpd_frozen_lopo_joint_model.py \\
        --expanded-table ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/master/expanded_trajectory_source_of_truth.csv \\
        --prpd-table ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/prpd_candidate_measurements.csv \\
        --output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/prpd_lopo_joint_model_results.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

plot_mod = importlib.import_module("26_plot_locked_claim_results")

LABEL = "late_pose_retained"
PROTEIN = "protein"
# 26_plot_locked_claim_results.py's metric_row/cluster_bootstrap_auc/
# paired_cluster_bootstrap_difference read these as module-level globals
# (they default to the locked 52-cohort's column names) -- override to
# match the expanded-cohort table used here.
plot_mod.LABEL = LABEL
plot_mod.PROTEIN = PROTEIN
N_BOOTSTRAP = 5000
SEED = 7
FALSE_STOP_BUDGET = 0.10
NEAR_BOUNDARY_HALFWIDTH_A = 1.5
LABEL_CUTOFF_A = 3.0

CHECKPOINTS = [
    ("5ns", "3_5", True),
    ("10ns", "8_10", True),
    ("20ns", "18_20", False),
]


def lopo_logistic_proba(
    df: pd.DataFrame, feature_cols: list[str], protein_col: str = PROTEIN, y_col: str = LABEL
) -> tuple[np.ndarray, list[dict]]:
    """Out-of-fold probabilities, fit/scaled on training proteins only per fold.

    Returns (oof_probability, fold_records). fold_records carries the
    standardized coefficients for each fold, for the coefficient-stability
    report -- irrelevant for single-feature models but harmless to collect.
    """
    proba = np.full(len(df), np.nan)
    fold_records = []
    for held_out in df[protein_col].unique():
        test_mask = (df[protein_col] == held_out).to_numpy()
        train_df = df.loc[~test_mask]
        if train_df[y_col].nunique() < 2:
            majority = float(train_df[y_col].mode().iloc[0])
            proba[test_mask] = majority
            continue
        scaler = StandardScaler()
        x_train = scaler.fit_transform(train_df[feature_cols].to_numpy())
        clf = LogisticRegression(C=1.0, penalty="l2", class_weight="balanced")
        clf.fit(x_train, train_df[y_col].to_numpy())
        x_test = scaler.transform(df.loc[test_mask, feature_cols].to_numpy())
        proba[test_mask] = clf.predict_proba(x_test)[:, 1]
        fold_records.append(
            {
                "held_out_protein": held_out,
                "n_train": len(train_df),
                "intercept": float(clf.intercept_[0]),
                **{
                    f"coef_{name}": float(coef)
                    for name, coef in zip(feature_cols, clf.coef_[0])
                },
            }
        )
    return proba, fold_records


def cluster_bootstrap_ba(
    table: pd.DataFrame, pred_col: str, n: int = N_BOOTSTRAP, seed: int = SEED
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    proteins = table[PROTEIN].unique()
    values = []
    for _ in range(n):
        sampled = rng.choice(proteins, size=len(proteins), replace=True)
        boot = pd.concat([table.loc[table[PROTEIN] == p] for p in sampled], ignore_index=True)
        y = boot[LABEL].to_numpy()
        if len(np.unique(y)) < 2:
            continue
        values.append(balanced_accuracy_score(y, boot[pred_col].to_numpy()))
    return tuple(np.percentile(values, [2.5, 97.5])) if values else (float("nan"), float("nan"))


def compute_saved_at_false_stop_budget(
    proba: np.ndarray, labels: np.ndarray, checkpoint_ns: float, budget: float = FALSE_STOP_BUDGET
) -> dict:
    thresholds = np.linspace(0.0, 1.0, 501)
    n_retained = int(np.sum(labels == 1))
    n_total = len(labels)
    best = {"threshold": float("nan"), "compute_saved_fraction": 0.0, "false_stop_rate": 0.0, "n_stopped": 0}
    for threshold in thresholds:
        stopped = proba < threshold
        n_stopped = int(stopped.sum())
        false_stop_rate = (
            float(np.sum(stopped & (labels == 1)) / n_retained) if n_retained else float("nan")
        )
        if false_stop_rate <= budget:
            compute_saved = n_stopped * (100.0 - checkpoint_ns) / (n_total * 100.0)
            if compute_saved > best["compute_saved_fraction"]:
                best = {
                    "threshold": float(threshold),
                    "compute_saved_fraction": float(compute_saved),
                    "false_stop_rate": false_stop_rate,
                    "n_stopped": n_stopped,
                }
    return best


def evaluate_model(
    table: pd.DataFrame,
    model_name: str,
    feature_cols: list[str],
    checkpoint_label: str,
    checkpoint_ns: float,
    rmsd_proba_by_complex_id: pd.Series | None,
) -> tuple[dict, pd.Series, list[dict]]:
    sub = table.dropna(subset=[LABEL] + feature_cols).copy()
    proba, fold_records = lopo_logistic_proba(sub, feature_cols)
    sub["_proba"] = proba
    labels = sub[LABEL].to_numpy(int)
    proba_by_complex_id = pd.Series(proba, index=sub["complex_id"].to_numpy())

    row = plot_mod.metric_row(sub, model_name, checkpoint_label, "_proba", 1.0)

    if rmsd_proba_by_complex_id is not None:
        # Join by complex_id, not position -- feature sets can in principle
        # drop different rows to NaN, so positional alignment would be a
        # silent correctness bug.
        sub["_rmsd_proba"] = sub["complex_id"].map(rmsd_proba_by_complex_id)
        diff, lo, hi = plot_mod.paired_cluster_bootstrap_difference(
            sub.dropna(subset=["_rmsd_proba"]), "_proba", 1.0, "_rmsd_proba", 1.0
        )
    else:
        diff, lo, hi = 0.0, 0.0, 0.0

    ba = balanced_accuracy_score(labels, (proba >= 0.5).astype(int))
    ba_lo, ba_hi = cluster_bootstrap_ba(sub.assign(_pred=(proba >= 0.5).astype(int)), "_pred")

    savings = compute_saved_at_false_stop_budget(proba, labels, checkpoint_ns)

    row.update(
        {
            "model": model_name,
            "checkpoint": checkpoint_label,
            "auroc_ci_lo": row["ci_lo"],
            "auroc_ci_hi": row["ci_hi"],
            "paired_diff_vs_rmsd": diff,
            "paired_diff_ci_lo": lo,
            "paired_diff_ci_hi": hi,
            "balanced_accuracy": ba,
            "balanced_accuracy_ci_lo": ba_lo,
            "balanced_accuracy_ci_hi": ba_hi,
            "compute_saved_at_fsr10": savings["compute_saved_fraction"],
            "achieved_false_stop_rate": savings["false_stop_rate"],
            "decision_threshold": savings["threshold"],
        }
    )
    return row, proba_by_complex_id, fold_records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expanded-table", type=Path, required=True)
    parser.add_argument("--prpd-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    expanded = pd.read_csv(args.expanded_table)
    expanded = expanded[expanded["eligible_primary_cohort"].astype(bool)].copy()
    prpd = pd.read_csv(args.prpd_table)
    prpd_cols = ["complex_id"] + [
        c for c in prpd.columns if c.startswith("prpd_") or c.startswith("baseline_mad")
    ]
    table = expanded.merge(prpd[prpd_cols], on="complex_id", how="inner", validate="one_to_one")
    print(f"Merged cohort: {len(table)} rows, {table[PROTEIN].nunique()} proteins")

    all_rows = []
    coefficient_rows = []

    for checkpoint_label, suffix, exploratory in CHECKPOINTS:
        rmsd_col = f"corrected_pose_rmsd_mean_{suffix}_A"
        t_col = f"prpd_T_mean_{suffix}"
        theta_col = f"prpd_Theta_mean_{suffix}"
        i_col = f"prpd_I_mean_{suffix}"
        zpose_col = f"prpd_Zpose_mean_{suffix}"

        print(f"\n=== checkpoint {checkpoint_label} ({'exploratory' if exploratory else 'primary'}) ===")

        rmsd_row, rmsd_proba, _ = evaluate_model(table, "rmsd_only", [rmsd_col], checkpoint_label, float(checkpoint_label.rstrip("ns")), None)
        rmsd_row["exploratory"] = exploratory
        all_rows.append(rmsd_row)

        models = [
            ("T_only", [t_col]),
            ("Theta_only", [theta_col]),
            ("T_plus_Theta", [t_col, theta_col]),
            ("T_plus_Theta_plus_I", [t_col, theta_col, i_col]),
            ("Zpose_lopo_calibrated", [zpose_col]),
        ]
        for model_name, feature_cols in models:
            row, proba, fold_records = evaluate_model(
                table, model_name, feature_cols, checkpoint_label,
                float(checkpoint_label.rstrip("ns")), rmsd_proba,
            )
            row["exploratory"] = exploratory
            all_rows.append(row)
            print(
                f"{model_name:>22}  AUROC={row['auroc']:.3f} [{row['auroc_ci_lo']:.3f},{row['auroc_ci_hi']:.3f}]"
                f"  diff_vs_RMSD={row['paired_diff_vs_rmsd']:+.3f} [{row['paired_diff_ci_lo']:+.3f},{row['paired_diff_ci_hi']:+.3f}]"
                f"  BA={row['balanced_accuracy']:.3f}"
                f"  saved@10%FSR={row['compute_saved_at_fsr10']:.3f}"
            )
            if model_name == "T_plus_Theta":
                for record in fold_records:
                    record["checkpoint"] = checkpoint_label
                    coefficient_rows.append(record)

        print(
            f"{'rmsd_only':>22}  AUROC={rmsd_row['auroc']:.3f} [{rmsd_row['auroc_ci_lo']:.3f},{rmsd_row['auroc_ci_hi']:.3f}]"
            f"  BA={rmsd_row['balanced_accuracy']:.3f}"
            f"  saved@10%FSR={rmsd_row['compute_saved_at_fsr10']:.3f}"
        )

        # Near-boundary subset: same definition as the locked pipeline's
        # Table non-obvious, restricted to this checkpoint's RMSD distance
        # from the 3 A label cutoff.
        near = table[(table[rmsd_col] - LABEL_CUTOFF_A).abs() <= NEAR_BOUNDARY_HALFWIDTH_A].copy()
        print(f"  near-boundary subset: n={len(near)}")
        if near[LABEL].nunique() == 2 and len(near) >= 10:
            for model_name, feature_cols in [("rmsd_only", [rmsd_col])] + models:
                sub = table.dropna(subset=[LABEL] + feature_cols).copy()
                proba, _ = lopo_logistic_proba(sub, feature_cols)
                sub["_proba"] = proba
                near_sub = sub.loc[sub["complex_id"].isin(near["complex_id"])]
                if near_sub[LABEL].nunique() < 2:
                    continue
                auc = plot_mod.auc_from_scores(near_sub[LABEL], near_sub["_proba"])
                print(f"    {model_name:>22}  near-boundary AUROC={auc:.3f}  n={len(near_sub)}")
                all_rows.append(
                    {
                        "model": model_name,
                        "checkpoint": checkpoint_label,
                        "metric_id": f"near_boundary_{model_name}_{checkpoint_label}",
                        "n": len(near_sub),
                        "auroc": auc,
                        "exploratory": exploratory,
                        "subset": "near_boundary",
                    }
                )

    results = pd.DataFrame(all_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")

    coeff_table = pd.DataFrame(coefficient_rows)
    coeff_path = args.output.parent / "prpd_t_theta_fold_coefficients.csv"
    coeff_table.to_csv(coeff_path, index=False)
    print(f"Wrote {coeff_path}")
    print("\nT+Theta coefficient stability across LOPO folds (standardized units):")
    for checkpoint_label, suffix, _ in CHECKPOINTS:
        sub = coeff_table[coeff_table["checkpoint"] == checkpoint_label]
        if sub.empty:
            continue
        coef_t_col = f"coef_prpd_T_mean_{suffix}"
        coef_theta_col = f"coef_prpd_Theta_mean_{suffix}"
        print(
            f"  {checkpoint_label} (n_folds={len(sub)}): "
            f"beta_T mean={sub[coef_t_col].mean():+.3f} std={sub[coef_t_col].std():.3f}   "
            f"beta_Theta mean={sub[coef_theta_col].mean():+.3f} std={sub[coef_theta_col].std():.3f}"
        )


if __name__ == "__main__":
    main()
