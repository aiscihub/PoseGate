#!/usr/bin/env python3
"""Evaluation for the local-pocket-alignment control (script 35).

Answers exactly one question: does Theta (rotation) survive when measured
relative to a frozen local pocket frame instead of the whole-protein Cα
frame? Reuses the frozen LOPO logistic spec from script 34
(lopo_logistic_proba / evaluate_model: training-fold-only standardization,
L2, C=1, balanced class weights, no tuning) so the comparison is on
identical footing to the whole-frame results already reported.

No new models beyond what the control calls for:
  Theta_whole, Theta_pocket(8A), T_whole, T_pocket(8A),
  frozen LOPO T_whole+Theta_whole, frozen LOPO T_pocket+Theta_pocket,
  RMSD benchmark.
6A/10A pocket radii are reported as standalone-AUROC sensitivity checks
only (not run through the full LOPO/compute-saved battery), per the
control's own instruction not to choose or over-explore radius.

Usage:
    python 36_local_pocket_control_evaluation.py \\
        --expanded-table ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/master/expanded_trajectory_source_of_truth.csv \\
        --pocket-table ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/local_pocket_control_measurements.csv \\
        --output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/local_pocket_control_results.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

plot_mod = importlib.import_module("26_plot_locked_claim_results")
model_mod = importlib.import_module("34_prpd_frozen_lopo_joint_model")

LABEL = "late_pose_retained"
PROTEIN = "protein"
plot_mod.LABEL = LABEL
plot_mod.PROTEIN = PROTEIN

CHECKPOINTS = [("5ns", "3_5"), ("10ns", "8_10"), ("20ns", "18_20")]
LABEL_CUTOFF_A = 3.0
NEAR_BOUNDARY_HALFWIDTH_A = 1.5


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--expanded-table", type=Path, required=True)
    parser.add_argument("--pocket-table", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    expanded = pd.read_csv(args.expanded_table)
    expanded = expanded[expanded["eligible_primary_cohort"].astype(bool)].copy()
    pocket = pd.read_csv(args.pocket_table)
    keep_cols = ["complex_id"] + [
        c for c in pocket.columns if c.startswith(("T_whole", "Theta_whole", "I_mean", "T_pocket", "Theta_pocket"))
    ]
    table = expanded.merge(pocket[keep_cols], on="complex_id", how="inner", validate="one_to_one")
    print(f"Merged cohort: {len(table)} rows, {table[PROTEIN].nunique()} proteins")

    all_rows = []
    coeff_rows = []
    rank_change_rows = []

    for checkpoint_label, suffix in CHECKPOINTS:
        checkpoint_ns = float(checkpoint_label.rstrip("ns"))
        rmsd_col = f"corrected_pose_rmsd_mean_{suffix}_A"
        t_whole_col = f"T_whole_mean_{suffix}"
        theta_whole_col = f"Theta_whole_mean_{suffix}"
        t_pocket_col = f"T_pocket_8A_mean_{suffix}"
        theta_pocket_col = f"Theta_pocket_8A_mean_{suffix}"

        print(f"\n=== checkpoint {checkpoint_label} ===")

        rmsd_row, rmsd_proba, _ = model_mod.evaluate_model(
            table, "rmsd_only", [rmsd_col], checkpoint_label, checkpoint_ns, None
        )
        all_rows.append(rmsd_row)
        print(
            f"{'rmsd_only':>28}  AUROC={rmsd_row['auroc']:.3f} [{rmsd_row['auroc_ci_lo']:.3f},{rmsd_row['auroc_ci_hi']:.3f}]"
            f"  BA={rmsd_row['balanced_accuracy']:.3f}  saved@10%FSR={rmsd_row['compute_saved_at_fsr10']:.3f}"
        )

        models = [
            ("T_whole_only", [t_whole_col]),
            ("Theta_whole_only", [theta_whole_col]),
            ("T_pocket_only", [t_pocket_col]),
            ("Theta_pocket_only", [theta_pocket_col]),
            ("T_whole_plus_Theta_whole", [t_whole_col, theta_whole_col]),
            ("T_pocket_plus_Theta_pocket", [t_pocket_col, theta_pocket_col]),
        ]
        for model_name, feature_cols in models:
            row, _, fold_records = model_mod.evaluate_model(
                table, model_name, feature_cols, checkpoint_label, checkpoint_ns, rmsd_proba
            )
            all_rows.append(row)
            print(
                f"{model_name:>28}  AUROC={row['auroc']:.3f} [{row['auroc_ci_lo']:.3f},{row['auroc_ci_hi']:.3f}]"
                f"  diff_vs_RMSD={row['paired_diff_vs_rmsd']:+.3f} [{row['paired_diff_ci_lo']:+.3f},{row['paired_diff_ci_hi']:+.3f}]"
                f"  BA={row['balanced_accuracy']:.3f}  saved@10%FSR={row['compute_saved_at_fsr10']:.3f}"
            )
            if model_name in ("T_whole_plus_Theta_whole", "T_pocket_plus_Theta_pocket"):
                for record in fold_records:
                    record["checkpoint"] = checkpoint_label
                    record["model"] = model_name
                    coeff_rows.append(record)

        # Near-boundary subset (same definition as the locked pipeline's Table non-obvious).
        near = table[(table[rmsd_col] - LABEL_CUTOFF_A).abs() <= NEAR_BOUNDARY_HALFWIDTH_A].copy()
        print(f"  near-boundary subset: n={len(near)}")
        if near[LABEL].nunique() == 2 and len(near) >= 10:
            for model_name, feature_cols in [("rmsd_only", [rmsd_col])] + models:
                sub = table.dropna(subset=[LABEL] + feature_cols).copy()
                proba, _ = model_mod.lopo_logistic_proba(sub, feature_cols)
                sub["_proba"] = proba
                near_sub = sub.loc[sub["complex_id"].isin(near["complex_id"])]
                if near_sub[LABEL].nunique() < 2:
                    continue
                auc = plot_mod.auc_from_scores(near_sub[LABEL], near_sub["_proba"])
                print(f"    {model_name:>28}  near-boundary AUROC={auc:.3f}  n={len(near_sub)}")
                all_rows.append(
                    {
                        "model": model_name,
                        "checkpoint": checkpoint_label,
                        "metric_id": f"near_boundary_{model_name}_{checkpoint_label}",
                        "n": len(near_sub),
                        "auroc": auc,
                        "subset": "near_boundary",
                    }
                )

        # Theta_whole vs Theta_pocket agreement.
        sub = table.dropna(subset=[theta_whole_col, theta_pocket_col]).copy()
        pearson_r = float(np.corrcoef(sub[theta_whole_col], sub[theta_pocket_col])[0, 1])
        spearman_r = float(spearmanr(sub[theta_whole_col], sub[theta_pocket_col]).correlation)
        rank_whole = sub[theta_whole_col].rank()
        rank_pocket = sub[theta_pocket_col].rank()
        rank_shift = (rank_whole - rank_pocket).abs()
        n_shifted_10 = int((rank_shift > 10).sum())
        print(
            f"  Theta_whole vs Theta_pocket(8A): pearson r={pearson_r:.3f}  spearman r={spearman_r:.3f}"
            f"  n(|rank shift|>10 of {len(sub)})={n_shifted_10}"
        )
        rank_change_rows.append(
            {
                "checkpoint": checkpoint_label,
                "n": len(sub),
                "pearson_r_theta_whole_vs_pocket": pearson_r,
                "spearman_r_theta_whole_vs_pocket": spearman_r,
                "n_rank_shifted_gt10": n_shifted_10,
            }
        )

        # 6A / 10A sensitivity check, standalone AUROC only.
        print("  pocket-radius sensitivity (Theta_pocket standalone AUROC):")
        for radius_tag in ("6A", "8A", "10A"):
            col = f"Theta_pocket_{radius_tag}_mean_{suffix}"
            if col not in table.columns:
                continue
            sub = table.dropna(subset=[LABEL, col])
            row = plot_mod.metric_row(sub, f"theta_pocket_{radius_tag}", checkpoint_label, col, -1.0)
            print(f"    {radius_tag:>4}  AUROC={row['auroc']:.3f}  n={row['n']}")
            all_rows.append(row | {"model": f"Theta_pocket_{radius_tag}", "checkpoint": checkpoint_label, "subset": "radius_sensitivity"})

    results = pd.DataFrame(all_rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(args.output, index=False)
    print(f"\nWrote {args.output}")

    coeff_table = pd.DataFrame(coeff_rows)
    coeff_path = args.output.parent / "local_pocket_control_fold_coefficients.csv"
    coeff_table.to_csv(coeff_path, index=False)
    print(f"Wrote {coeff_path}")

    rank_table = pd.DataFrame(rank_change_rows)
    rank_path = args.output.parent / "local_pocket_control_theta_agreement.csv"
    rank_table.to_csv(rank_path, index=False)
    print(f"Wrote {rank_path}")

    print("\nCoefficient stability across LOPO folds (standardized units):")
    for model_name in ("T_whole_plus_Theta_whole", "T_pocket_plus_Theta_pocket"):
        for checkpoint_label, _ in CHECKPOINTS:
            sub = coeff_table[(coeff_table["model"] == model_name) & (coeff_table["checkpoint"] == checkpoint_label)]
            if sub.empty:
                continue
            coef_cols = [c for c in sub.columns if c.startswith("coef_")]
            t_col, theta_col = coef_cols[0], coef_cols[1]
            n_t_negative = int((sub[t_col] < 0).sum())
            n_theta_negative = int((sub[theta_col] < 0).sum())
            print(
                f"  {model_name} @ {checkpoint_label} (n_folds={len(sub)}): "
                f"beta_T mean={sub[t_col].mean():+.3f} std={sub[t_col].std():.3f} (negative in {n_t_negative}/{len(sub)} folds)   "
                f"beta_Theta mean={sub[theta_col].mean():+.3f} std={sub[theta_col].std():.3f} (negative in {n_theta_negative}/{len(sub)} folds)"
            )


if __name__ == "__main__":
    main()
