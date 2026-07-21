#!/usr/bin/env python3
"""Agreement analysis for the local-pocket-alignment control (check 2 + 4).

Beyond pooled Pearson correlation (which can be inflated by between-
protein differences), reports: mean/median absolute difference, linear
regression slope/intercept, Lin's concordance correlation coefficient
(CCC), Bland-Altman limits of agreement, and within-protein-demeaned
correlation (residual agreement after removing each protein's own mean --
this is what tests whether agreement survives INSIDE each transporter
group, not just across the whole heterogeneous panel).

Also (check 4) compares agreement at three granularities using the
per-frame trace from script 37: individual saved frames, a single
nearest-checkpoint frame (no window averaging), and the 2 ns
window-averaged checkpoint value -- to test whether window averaging is
doing the work of suppressing pocket-fit noise, as hypothesized.

Produces a two-panel figure: Theta_whole vs Theta_pocket scatter (colored
by protein, with identity line) and the paired Bland-Altman plot.

Usage:
    python 38_pocket_frame_agreement_analysis.py \\
        --pocket-table ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/local_pocket_control_measurements.csv \\
        --frame-trace ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/pocket_frame_diagnostic_traces.csv \\
        --diagnostic-summary ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/pocket_frame_diagnostic_summary.csv \\
        --output-dir ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CHECKPOINTS = [("5ns", "3_5"), ("10ns", "8_10"), ("20ns", "18_20")]


def concordance_correlation_coefficient(x: np.ndarray, y: np.ndarray) -> float:
    mean_x, mean_y = x.mean(), y.mean()
    var_x, var_y = x.var(), y.var()
    covariance = np.mean((x - mean_x) * (y - mean_y))
    return float(2 * covariance / (var_x + var_y + (mean_x - mean_y) ** 2))


def agreement_stats(x: np.ndarray, y: np.ndarray, protein: np.ndarray) -> dict:
    diff = y - x
    slope, intercept = np.polyfit(x, y, 1)
    pooled_r = float(np.corrcoef(x, y)[0, 1])

    protein_series = pd.Series(protein)
    x_demeaned = x - protein_series.map(pd.Series(x).groupby(protein_series).mean()).to_numpy()
    y_demeaned = y - protein_series.map(pd.Series(y).groupby(protein_series).mean()).to_numpy()
    within_protein_r = float(np.corrcoef(x_demeaned, y_demeaned)[0, 1])

    return {
        "n": len(x),
        "pooled_pearson_r": pooled_r,
        "within_protein_pearson_r": within_protein_r,
        "mean_abs_diff": float(np.mean(np.abs(diff))),
        "median_abs_diff": float(np.median(np.abs(diff))),
        "regression_slope": float(slope),
        "regression_intercept": float(intercept),
        "ccc": concordance_correlation_coefficient(x, y),
        "bland_altman_mean_diff": float(diff.mean()),
        "bland_altman_sd_diff": float(diff.std()),
        "bland_altman_loa_lower": float(diff.mean() - 1.96 * diff.std()),
        "bland_altman_loa_upper": float(diff.mean() + 1.96 * diff.std()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pocket-table", type=Path, required=True)
    parser.add_argument("--frame-trace", type=Path, required=True)
    parser.add_argument("--diagnostic-summary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    pocket = pd.read_csv(args.pocket_table)
    pocket = pocket[pocket["status"].eq("ok")].copy()

    print("=== Checkpoint-level agreement: Theta_whole vs Theta_pocket(8A) ===")
    rows = []
    for checkpoint_label, suffix in CHECKPOINTS:
        whole_col = f"Theta_whole_mean_{suffix}"
        pocket_col = f"Theta_pocket_8A_mean_{suffix}"
        sub = pocket.dropna(subset=[whole_col, pocket_col])
        stats = agreement_stats(sub[whole_col].to_numpy(), sub[pocket_col].to_numpy(), sub["protein"].to_numpy())
        stats["checkpoint"] = checkpoint_label
        rows.append(stats)
        print(
            f"{checkpoint_label:>5}  n={stats['n']}  pooled_r={stats['pooled_pearson_r']:.4f}  "
            f"within_protein_r={stats['within_protein_pearson_r']:.4f}  "
            f"CCC={stats['ccc']:.4f}  mean|diff|={stats['mean_abs_diff']:.4f} rad  "
            f"slope={stats['regression_slope']:.3f}  intercept={stats['regression_intercept']:.3f}  "
            f"BA limits=[{stats['bland_altman_loa_lower']:+.4f}, {stats['bland_altman_loa_upper']:+.4f}] rad"
        )

    agreement_table = pd.DataFrame(rows)
    agreement_path = args.output_dir / "pocket_frame_agreement_stats.csv"
    agreement_table.to_csv(agreement_path, index=False)
    print(f"\nWrote {agreement_path}")

    # --- Check 4: frame-level vs checkpoint-window-averaged agreement ---
    print("\n=== Check 4: agreement at frame-level vs single-checkpoint-frame vs 2ns-window-averaged (20ns) ===")
    frames = pd.read_csv(args.frame_trace)
    frame_level = frames[frames["time_ns"] > 0.05]  # drop the trivial frame-0 zero row
    r_frame_level = float(np.corrcoef(frame_level["Theta_whole"], frame_level["Theta_pocket"])[0, 1])
    print(f"  individual-frame level (all {len(frame_level)} saved frames): pearson r = {r_frame_level:.4f}")

    nearest_20ns = frames.iloc[(frames["time_ns"] - 20.0).abs().groupby(frames["complex_id"]).idxmin()]
    r_single_checkpoint = float(np.corrcoef(nearest_20ns["Theta_whole"], nearest_20ns["Theta_pocket"])[0, 1])
    print(f"  single nearest-frame at 20ns (n={len(nearest_20ns)}): pearson r = {r_single_checkpoint:.4f}")

    window_row = agreement_table[agreement_table["checkpoint"] == "20ns"].iloc[0]
    print(f"  2ns window-averaged at 20ns (n={int(window_row['n'])}): pearson r = {window_row['pooled_pearson_r']:.4f}")

    # --- Outlier investigation: worst-divergence and worst-conditioned trajectories ---
    print("\n=== Outlier check: does agreement hold for the most stressed cases? ===")
    diag = pd.read_csv(args.diagnostic_summary)
    diag = diag[diag["status"].eq("ok")]
    worst_divergence = diag.loc[diag["divergence_max_deg"].idxmax()]
    worst_conditioning = diag.loc[diag["pocket_conditioning_ratio"].idxmin()]
    print(
        f"  worst angular divergence: {worst_divergence['complex_id']} "
        f"(max={worst_divergence['divergence_max_deg']:.1f} deg, median={worst_divergence['divergence_median_deg']:.1f} deg)"
    )
    print(
        f"  worst pocket conditioning: {worst_conditioning['complex_id']} "
        f"(ratio={worst_conditioning['pocket_conditioning_ratio']:.4f}, n_pocket_ca={worst_conditioning['n_pocket_ca']})"
    )
    for label, cid in [("worst-divergence", worst_divergence["complex_id"]), ("worst-conditioning", worst_conditioning["complex_id"])]:
        row = pocket[pocket["complex_id"] == cid]
        if row.empty:
            continue
        row = row.iloc[0]
        for checkpoint_label, suffix in CHECKPOINTS:
            wv = row.get(f"Theta_whole_mean_{suffix}")
            pv = row.get(f"Theta_pocket_8A_mean_{suffix}")
            if pd.notna(wv) and pd.notna(pv):
                print(
                    f"    [{label}] {cid} @ {checkpoint_label}: "
                    f"Theta_whole={np.degrees(wv):.2f} deg  Theta_pocket={np.degrees(pv):.2f} deg  "
                    f"abs_diff={np.degrees(abs(wv-pv)):.2f} deg"
                )

    # --- Two-panel figure ---
    suffix_20ns = "18_20"
    whole_col = f"Theta_whole_mean_{suffix_20ns}"
    pocket_col = f"Theta_pocket_8A_mean_{suffix_20ns}"
    sub = pocket.dropna(subset=[whole_col, pocket_col]).copy()
    sub["theta_whole_deg"] = np.degrees(sub[whole_col])
    sub["theta_pocket_deg"] = np.degrees(sub[pocket_col])
    proteins = sorted(sub["protein"].unique())
    colors = plt.cm.tab20(np.linspace(0, 1, len(proteins)))
    color_map = dict(zip(proteins, colors))

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    ax = axes[0]
    for protein in proteins:
        rows_p = sub[sub["protein"] == protein]
        ax.scatter(rows_p["theta_whole_deg"], rows_p["theta_pocket_deg"], color=color_map[protein], label=protein, s=30, alpha=0.85)
    lims = [0, max(sub["theta_whole_deg"].max(), sub["theta_pocket_deg"].max()) * 1.05]
    ax.plot(lims, lims, "k--", linewidth=1, label="identity")
    ax.set_xlabel("Theta_whole (deg, 18-20ns)")
    ax.set_ylabel("Theta_pocket 8A (deg, 18-20ns)")
    ax.set_title("Whole-frame vs pocket-frame reorientation")
    ax.legend(fontsize=6, ncol=2, loc="upper left")

    ax = axes[1]
    mean_val = (sub["theta_whole_deg"] + sub["theta_pocket_deg"]) / 2
    diff_val = sub["theta_pocket_deg"] - sub["theta_whole_deg"]
    for protein in proteins:
        rows_p = sub[sub["protein"] == protein]
        m = (rows_p["theta_whole_deg"] + rows_p["theta_pocket_deg"]) / 2
        d = rows_p["theta_pocket_deg"] - rows_p["theta_whole_deg"]
        ax.scatter(m, d, color=color_map[protein], s=30, alpha=0.85)
    mean_diff = diff_val.mean()
    sd_diff = diff_val.std()
    ax.axhline(mean_diff, color="k", linestyle="-", linewidth=1)
    ax.axhline(mean_diff + 1.96 * sd_diff, color="k", linestyle="--", linewidth=1)
    ax.axhline(mean_diff - 1.96 * sd_diff, color="k", linestyle="--", linewidth=1)
    ax.set_xlabel("mean(Theta_whole, Theta_pocket) (deg)")
    ax.set_ylabel("Theta_pocket - Theta_whole (deg)")
    ax.set_title("Bland-Altman: pocket vs whole agreement")

    fig.tight_layout()
    figure_path = args.output_dir / "pocket_frame_agreement_figure.png"
    fig.savefig(figure_path, dpi=150)
    print(f"\nWrote {figure_path}")


if __name__ == "__main__":
    main()
