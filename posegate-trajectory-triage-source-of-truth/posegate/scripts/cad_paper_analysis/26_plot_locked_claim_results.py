#!/usr/bin/env python3
"""Generate the key claim figures from the locked source-of-truth table.

The plots use one strict complete-case cohort throughout:

  - complete nominal (70,100] ns target coverage;
  - locked RMSD-only late outcome;
  - corrected same-trajectory 18--20 ns pose features available.

The script does not read aggregate legacy result tables. Every plotted value is
recomputed from `same_trajectory_source_of_truth.csv`. Confidence intervals
resample whole proteins with replacement while keeping each protein's complexes
together. Continuous-score intervals do not refit a model; the policy plots use
the already generated outer-LOPO probabilities stored in the audit table.

Outputs are written to:

pipeline/training/outputs/yau_recapture_screen/cad_states/
same_trajectory_source_of_truth/claim_figures/

Usage:
    MPLCONFIGDIR=/tmp/matplotlib-cache python \
        yau_competition/scripts/cad_paper_analysis/26_plot_locked_claim_results.py
"""

from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import OUT_DIR


SOURCE_TABLE = (
    OUT_DIR
    / "same_trajectory_source_of_truth"
    / "same_trajectory_source_of_truth.csv"
)
RESULT_DIR = OUT_DIR / "same_trajectory_source_of_truth" / "claim_figures"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

LABEL = "late_pose_retained_locked"
PROTEIN = "protein"
N_BOOTSTRAP = 5000
SEED = 7

NAVY = "#173F73"
BLUE = "#3D78B5"
GREEN = "#2E7D4A"
ORANGE = "#D97706"
RED = "#B6403A"
GRAY = "#6B7280"
LIGHT_GRAY = "#E5E7EB"


def auc_from_scores(labels: np.ndarray, scores: np.ndarray) -> float:
    """AUROC by pairwise ordering; higher score means more likely retained."""
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    positive = scores[labels == 1]
    negative = scores[labels == 0]
    if len(positive) == 0 or len(negative) == 0:
        return float("nan")
    comparisons = positive[:, None] - negative[None, :]
    return float((np.sum(comparisons > 0) + 0.5 * np.sum(comparisons == 0)) / comparisons.size)


def cluster_bootstrap_auc(
    table: pd.DataFrame,
    score_column: str,
    sign: float,
    n: int = N_BOOTSTRAP,
    seed: int = SEED,
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    proteins = table[PROTEIN].dropna().unique()
    blocks = [
        (
            table.loc[table[PROTEIN].eq(protein), LABEL].to_numpy(int),
            sign * table.loc[table[PROTEIN].eq(protein), score_column].to_numpy(float),
        )
        for protein in proteins
    ]
    values: list[float] = []
    for _ in range(n):
        sampled = rng.integers(0, len(blocks), size=len(blocks))
        labels = np.concatenate([blocks[index][0] for index in sampled])
        scores = np.concatenate([blocks[index][1] for index in sampled])
        auc = auc_from_scores(labels, scores)
        if np.isfinite(auc):
            values.append(auc)
    if not values:
        return float("nan"), float("nan")
    return tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))


def paired_cluster_bootstrap_difference(
    table: pd.DataFrame,
    score_a: str,
    sign_a: float,
    score_b: str,
    sign_b: float,
    n: int = N_BOOTSTRAP,
    seed: int = SEED,
) -> tuple[float, float, float]:
    observed = auc_from_scores(table[LABEL], sign_a * table[score_a]) - auc_from_scores(
        table[LABEL], sign_b * table[score_b]
    )
    rng = np.random.default_rng(seed)
    proteins = table[PROTEIN].dropna().unique()
    blocks = [
        (
            table.loc[table[PROTEIN].eq(protein), LABEL].to_numpy(int),
            sign_a * table.loc[table[PROTEIN].eq(protein), score_a].to_numpy(float),
            sign_b * table.loc[table[PROTEIN].eq(protein), score_b].to_numpy(float),
        )
        for protein in proteins
    ]
    values: list[float] = []
    for _ in range(n):
        sampled = rng.integers(0, len(blocks), size=len(blocks))
        labels = np.concatenate([blocks[index][0] for index in sampled])
        scores_a = np.concatenate([blocks[index][1] for index in sampled])
        scores_b = np.concatenate([blocks[index][2] for index in sampled])
        auc_a = auc_from_scores(labels, scores_a)
        auc_b = auc_from_scores(labels, scores_b)
        if np.isfinite(auc_a) and np.isfinite(auc_b):
            values.append(auc_a - auc_b)
    lo, hi = np.percentile(values, [2.5, 97.5]) if values else (np.nan, np.nan)
    return float(observed), float(lo), float(hi)


def metric_row(
    table: pd.DataFrame,
    metric_id: str,
    task: str,
    score_column: str,
    sign: float,
    seed_offset: int = 0,
) -> dict:
    complete = table.dropna(subset=[LABEL, score_column]).copy()
    auc = auc_from_scores(complete[LABEL], sign * complete[score_column])
    ci_lo, ci_hi = cluster_bootstrap_auc(
        complete, score_column, sign, seed=SEED + seed_offset
    )
    return {
        "metric_id": metric_id,
        "task": task,
        "score_column": score_column,
        "higher_score_means": "late retained" if sign > 0 else "late non-retained",
        "n": len(complete),
        "n_proteins": complete[PROTEIN].nunique(),
        "n_retained": int(complete[LABEL].sum()),
        "n_nonretained": int((complete[LABEL] == 0).sum()),
        "auroc": auc,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "ci_method": f"protein-clustered bootstrap, {N_BOOTSTRAP} draws; fixed score",
    }


def build_policy_frontier(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    probability = table["locked_lopo_probability_retained"].to_numpy(float)
    labels = table[LABEL].to_numpy(int)
    thresholds = np.unique(np.concatenate([np.linspace(0.0, 1.0, 201), [0.5]]))
    rows = []
    for threshold in thresholds:
        decision_continue = probability >= threshold
        stopped = ~decision_continue
        tn = int(np.sum(stopped & (labels == 0)))
        fn = int(np.sum(stopped & (labels == 1)))
        tp = int(np.sum(decision_continue & (labels == 1)))
        fp = int(np.sum(decision_continue & (labels == 0)))
        n_retained = int(np.sum(labels == 1))
        n_nonretained = int(np.sum(labels == 0))
        n_stopped = int(stopped.sum())
        rows.append(
            {
                "probability_threshold_continue": threshold,
                "n_stopped": n_stopped,
                "n_continued": len(table) - n_stopped,
                "correct_stop_tn": tn,
                "false_stop_fn": fn,
                "correct_continue_tp": tp,
                "nonretained_continue_fp": fp,
                "false_stop_rate": fn / n_retained,
                "nonretained_stop_rate": tn / n_nonretained,
                "compute_saved_fraction": n_stopped * 80.0 / (len(table) * 100.0),
            }
        )
    frontier = pd.DataFrame(rows)
    operating = frontier.loc[
        np.isclose(frontier["probability_threshold_continue"], 0.5)
    ].copy()
    return frontier, operating


def plot_key_claim(
    table: pd.DataFrame,
    metrics: pd.DataFrame,
    frontier: pd.DataFrame,
    operating: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(14, 10.5))
    fig.suptitle(
        "Locked same-trajectory claim audit: complete 70–100 ns cohort",
        fontsize=18,
        fontweight="bold",
        color=NAVY,
        y=0.98,
    )

    # A. Direct within-run persistence.
    ax = axes[0, 0]
    for label, name, color, marker in [
        (0, "Late non-retained", RED, "o"),
        (1, "Late retained", GREEN, "^"),
    ]:
        subset = table[table[LABEL].eq(label)]
        ax.scatter(
            subset["same_traj_pose_rmsd_18_20_A"],
            subset["late_pose_rmsd_median_last30pct_A"],
            s=52,
            alpha=0.82,
            color=color,
            edgecolor="white",
            linewidth=0.6,
            marker=marker,
            label=f"{name} (n={len(subset)})",
        )
    rho, rho_p = spearmanr(
        table["same_traj_pose_rmsd_18_20_A"],
        table["late_pose_rmsd_median_last30pct_A"],
    )
    same_auc = metrics.loc[metrics["metric_id"].eq("same_corrected_pose"), "auroc"].iloc[0]
    ax.axhline(3.0, color=GRAY, linestyle="--", linewidth=1.4, label="Locked late RMSD threshold")
    ax.set_xlabel("Corrected ligand RMSD, mean over (18,20] ns (Å)")
    ax.set_ylabel("Corrected ligand RMSD, median over (70,100] ns (Å)")
    ax.set_title("A. Early motion persists within the same run", loc="left", fontweight="bold")
    ax.text(
        0.98,
        0.04,
        f"Spearman ρ={rho:.2f} (p={rho_p:.2g})\nOutcome AUROC from early RMSD={same_auc:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=LIGHT_GRAY),
    )
    ax.legend(frameon=False, fontsize=9, loc="upper left")

    # B. Coordinate representations.
    ax = axes[0, 1]
    coordinate_ids = ["same_corrected_pose", "same_no_pbc_pose", "same_ligand_self_rmsd"]
    coordinate_names = ["Protein-aligned\n+ PBC-aware", "Protein-aligned\nPBC-naive", "Ligand-self\naligned"]
    coordinate_colors = [NAVY, ORANGE, GRAY]
    coord = metrics.set_index("metric_id").loc[coordinate_ids]
    x = np.arange(len(coord))
    yerr = np.vstack([coord["auroc"] - coord["ci_lo"], coord["ci_hi"] - coord["auroc"]])
    ax.bar(x, coord["auroc"], color=coordinate_colors, width=0.67, zorder=2)
    ax.errorbar(x, coord["auroc"], yerr=yerr, fmt="none", ecolor="#1F2937", capsize=5, zorder=3)
    ax.axhline(0.5, color=GRAY, linestyle="--", linewidth=1)
    ax.set_xticks(x, coordinate_names)
    ax.set_ylim(0.35, 1.02)
    ax.set_ylabel("AUROC for locked late retention")
    ax.set_title("B. Protein-relative, PBC-aware coordinates matter", loc="left", fontweight="bold")
    for xi, value in zip(x, coord["auroc"]):
        ax.text(xi, value + 0.025, f"{value:.3f}", ha="center", va="bottom", fontweight="bold")

    # C. Same trajectory versus separate launch.
    ax = axes[1, 0]
    comparison = [
        ("same_corrected_pose", "separate_corrected_pose", "Corrected pose RMSD"),
        ("same_pocket_retention", "separate_pocket_retention", "Pocket retention"),
    ]
    x = np.arange(len(comparison))
    width = 0.34
    same_rows = metrics.set_index("metric_id").loc[[item[0] for item in comparison]]
    separate_rows = metrics.set_index("metric_id").loc[[item[1] for item in comparison]]
    same_err = np.vstack(
        [same_rows["auroc"] - same_rows["ci_lo"], same_rows["ci_hi"] - same_rows["auroc"]]
    )
    separate_err = np.vstack(
        [
            separate_rows["auroc"] - separate_rows["ci_lo"],
            separate_rows["ci_hi"] - separate_rows["auroc"],
        ]
    )
    ax.bar(x - width / 2, same_rows["auroc"], width, color=NAVY, label="Same running trajectory")
    ax.bar(x + width / 2, separate_rows["auroc"], width, color=LIGHT_GRAY, edgecolor=GRAY, label="Separate launch")
    ax.errorbar(x - width / 2, same_rows["auroc"], yerr=same_err, fmt="none", ecolor="#111827", capsize=4)
    ax.errorbar(x + width / 2, separate_rows["auroc"], yerr=separate_err, fmt="none", ecolor="#111827", capsize=4)
    ax.axhline(0.5, color=GRAY, linestyle="--", linewidth=1)
    ax.set_xticks(x, [item[2] for item in comparison])
    ax.set_ylim(0.3, 1.02)
    ax.set_ylabel("AUROC for locked late retention")
    ax.set_title("C. Same-run advantage is strongest for corrected pose", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=9, loc="upper right")

    # D. Retrospective policy frontier.
    ax = axes[1, 1]
    ax.plot(
        100 * frontier["false_stop_rate"],
        100 * frontier["compute_saved_fraction"],
        color=NAVY,
        linewidth=2.4,
    )
    point = operating.iloc[0]
    ax.scatter(
        [100 * point["false_stop_rate"]],
        [100 * point["compute_saved_fraction"]],
        s=90,
        color=ORANGE,
        edgecolor="white",
        linewidth=1,
        zorder=4,
    )
    ax.annotate(
        (
            "p=0.5 retrospective point\n"
            f"saved={100*point['compute_saved_fraction']:.1f}%\n"
            f"false stops={int(point['false_stop_fn'])}/{int(table[LABEL].sum())}\n"
            f"non-retained stopped={int(point['correct_stop_tn'])}/{int((table[LABEL] == 0).sum())}"
        ),
        xy=(100 * point["false_stop_rate"], 100 * point["compute_saved_fraction"]),
        xytext=(8, -58),
        textcoords="offset points",
        fontsize=9,
        bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=LIGHT_GRAY),
        arrowprops=dict(arrowstyle="->", color=GRAY),
    )
    ax.set_xlim(-1, 101)
    ax.set_ylim(-1, 82)
    ax.set_xlabel("False-stop rate among late-retained trajectories (%)")
    ax.set_ylabel("Simulated time saved (%)")
    ax.set_title("D. Compute saving has an explicit error tradeoff", loc="left", fontweight="bold")

    for ax in axes.ravel():
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.7, alpha=0.7, zorder=0)

    fig.text(
        0.5,
        0.012,
        (
            f"Locked cohort: n={len(table)} trajectories from {table[PROTEIN].nunique()} proteins. "
            f"Intervals: {N_BOOTSTRAP}-draw protein-clustered bootstrap. "
            "Policy frontier is retrospective and not a prospective guarantee."
        ),
        ha="center",
        fontsize=9.5,
        color="#374151",
    )
    fig.tight_layout(rect=(0, 0.04, 1, 0.95), h_pad=2.0, w_pad=1.8)
    fig.savefig(RESULT_DIR / "figure_key_locked_claim.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULT_DIR / "figure_key_locked_claim.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_checkpoint_and_oof(
    table: pd.DataFrame,
    checkpoint_metrics: pd.DataFrame,
    operating: pd.DataFrame,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))

    ax = axes[0]
    checkpoints = checkpoint_metrics["checkpoint_ns"].to_numpy()
    aucs = checkpoint_metrics["auroc"].to_numpy()
    yerr = np.vstack(
        [
            aucs - checkpoint_metrics["ci_lo"].to_numpy(),
            checkpoint_metrics["ci_hi"].to_numpy() - aucs,
        ]
    )
    ax.errorbar(
        checkpoints,
        aucs,
        yerr=yerr,
        color=NAVY,
        marker="o",
        markersize=7,
        linewidth=2.2,
        capsize=4,
    )
    ax.axhline(0.5, color=GRAY, linestyle="--", linewidth=1)
    ax.set_xticks(checkpoints)
    ax.set_ylim(0.5, 1.02)
    ax.set_xlabel("Checkpoint (ns)")
    ax.set_ylabel("AUROC from corrected RMSD")
    ax.set_title(
        "A. Exploratory checkpoint comparison, locked common cohort",
        loc="left",
        fontweight="bold",
    )

    ax = axes[1]
    rng = np.random.default_rng(SEED)
    groups = []
    positions = [0, 1]
    colors = [RED, GREEN]
    labels = ["Late non-retained", "Late retained"]
    for outcome in [0, 1]:
        values = table.loc[table[LABEL].eq(outcome), "locked_lopo_probability_retained"].to_numpy()
        groups.append(values)
    box = ax.boxplot(groups, positions=positions, widths=0.5, patch_artist=True, showfliers=False)
    for patch, color in zip(box["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.22)
        patch.set_edgecolor(color)
    for position, values, color in zip(positions, groups, colors):
        jitter = rng.normal(0, 0.055, size=len(values))
        ax.scatter(position + jitter, values, color=color, s=38, alpha=0.78, edgecolor="white", linewidth=0.5)
    ax.axhline(0.5, color=GRAY, linestyle="--", linewidth=1.2, label="Retrospective decision threshold")
    ax.set_xticks(positions, [f"{label}\n(n={len(values)})" for label, values in zip(labels, groups)])
    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Outer-LOPO model score for late retention")
    ax.set_title("B. Uncalibrated OOF scores and the p=0.5 decision", loc="left", fontweight="bold")
    point = operating.iloc[0]
    ax.text(
        0.03,
        0.97,
        (
            f"At p=0.5: TN={int(point['correct_stop_tn'])}, FN={int(point['false_stop_fn'])}, "
            f"TP={int(point['correct_continue_tp'])}, FP={int(point['nonretained_continue_fp'])}"
        ),
        transform=ax.transAxes,
        va="top",
        fontsize=9,
    )

    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", color=LIGHT_GRAY, linewidth=0.7, alpha=0.7)
    fig.suptitle(
        "Supporting locked-cohort results",
        fontsize=16,
        fontweight="bold",
        color=NAVY,
        y=1.01,
    )
    fig.tight_layout()
    fig.savefig(RESULT_DIR / "figure_checkpoint_and_oof.png", dpi=300, bbox_inches="tight")
    fig.savefig(RESULT_DIR / "figure_checkpoint_and_oof.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    if not SOURCE_TABLE.exists():
        raise FileNotFoundError(
            f"missing {SOURCE_TABLE}; run 25_build_same_trajectory_source_of_truth.py first"
        )
    source = pd.read_csv(SOURCE_TABLE)
    table = source.loc[source["eligible_locked_primary_pose_model"].eq(True)].copy()
    if len(table) != 52:
        raise AssertionError(f"expected locked cohort n=52, found n={len(table)}")
    if table[PROTEIN].nunique() != 15:
        raise AssertionError("locked cohort no longer contains the expected 15 proteins")

    metric_specs = [
        ("same_corrected_pose", "same trajectory", "same_traj_pose_rmsd_18_20_A", -1.0),
        ("same_pocket_retention", "same trajectory", "same_traj_pocket_retention_18_20", 1.0),
        ("separate_corrected_pose", "separate launch", "separate_launch_pose_rmsd_18_20_A", -1.0),
        ("separate_pocket_retention", "separate launch", "separate_launch_pocket_retention_18_20", 1.0),
        ("same_no_pbc_pose", "coordinate control", "baseline_no_pbc_pose_rmsd_18_20_A", -1.0),
        ("same_ligand_self_rmsd", "coordinate control", "baseline_ligand_self_aligned_rmsd_18_20_A", -1.0),
        ("locked_oof_probability", "same trajectory model", "locked_lopo_probability_retained", 1.0),
    ]
    metric_rows = [
        metric_row(table, metric_id, task, column, sign, seed_offset=index)
        for index, (metric_id, task, column, sign) in enumerate(metric_specs)
    ]
    metrics = pd.DataFrame(metric_rows)
    metrics.to_csv(RESULT_DIR / "key_claim_metrics.csv", index=False)

    checkpoint_rows = []
    for index, checkpoint in enumerate([5, 10, 15, 20, 30]):
        column = f"same_traj_pose_rmsd_{checkpoint}ns_A"
        row = metric_row(
            table,
            f"same_corrected_pose_{checkpoint}ns",
            "same trajectory checkpoint",
            column,
            -1.0,
            # Reuse the same protein-bootstrap draws for every checkpoint.
            # This also makes the duplicated 20 ns interval identical to the
            # primary corrected-pose row in key_claim_metrics.csv.
            seed_offset=0,
        )
        row["checkpoint_ns"] = checkpoint
        checkpoint_rows.append(row)
    checkpoint_metrics = pd.DataFrame(checkpoint_rows)
    checkpoint_metrics.to_csv(RESULT_DIR / "checkpoint_locked_common_cohort.csv", index=False)

    paired_specs = [
        (
            "corrected_minus_no_pbc",
            "same_traj_pose_rmsd_18_20_A",
            -1.0,
            "baseline_no_pbc_pose_rmsd_18_20_A",
            -1.0,
        ),
        (
            "corrected_minus_ligand_self",
            "same_traj_pose_rmsd_18_20_A",
            -1.0,
            "baseline_ligand_self_aligned_rmsd_18_20_A",
            -1.0,
        ),
        (
            "same_minus_separate_pose",
            "same_traj_pose_rmsd_18_20_A",
            -1.0,
            "separate_launch_pose_rmsd_18_20_A",
            -1.0,
        ),
        (
            "same_minus_separate_retention",
            "same_traj_pocket_retention_18_20",
            1.0,
            "separate_launch_pocket_retention_18_20",
            1.0,
        ),
    ]
    paired_rows = []
    for index, (comparison, score_a, sign_a, score_b, sign_b) in enumerate(paired_specs):
        complete = table.dropna(subset=[LABEL, score_a, score_b]).copy()
        difference, ci_lo, ci_hi = paired_cluster_bootstrap_difference(
            complete,
            score_a,
            sign_a,
            score_b,
            sign_b,
            seed=SEED + 100 + index,
        )
        paired_rows.append(
            {
                "comparison": comparison,
                "score_a": score_a,
                "score_b": score_b,
                "n": len(complete),
                "n_proteins": complete[PROTEIN].nunique(),
                "auroc_difference_a_minus_b": difference,
                "ci_lo": ci_lo,
                "ci_hi": ci_hi,
                "ci_method": f"paired protein-clustered bootstrap, {N_BOOTSTRAP} draws",
            }
        )
    pd.DataFrame(paired_rows).to_csv(RESULT_DIR / "paired_auroc_differences.csv", index=False)

    frontier, operating = build_policy_frontier(table)
    frontier.to_csv(RESULT_DIR / "retrospective_policy_frontier.csv", index=False)
    operating.to_csv(RESULT_DIR / "retrospective_policy_operating_point_p05.csv", index=False)

    confusion = pd.DataFrame(
        [
            {
                "decision": "stop",
                "late_retained": int(operating.iloc[0]["false_stop_fn"]),
                "late_nonretained": int(operating.iloc[0]["correct_stop_tn"]),
            },
            {
                "decision": "continue",
                "late_retained": int(operating.iloc[0]["correct_continue_tp"]),
                "late_nonretained": int(operating.iloc[0]["nonretained_continue_fp"]),
            },
        ]
    )
    confusion.to_csv(RESULT_DIR / "retrospective_confusion_matrix_p05.csv", index=False)

    plot_key_claim(table, metrics, frontier, operating)
    plot_checkpoint_and_oof(table, checkpoint_metrics, operating)

    readme = f"""# Locked claim figures

Generated from `same_trajectory_source_of_truth.csv` by
`26_plot_locked_claim_results.py`.

All primary comparisons use the identical locked cohort: {len(table)} trajectories,
{table[PROTEIN].nunique()} proteins, {int(table[LABEL].sum())} late retained and
{int((table[LABEL] == 0).sum())} late non-retained. Seven truncated trajectories
from the legacy outcome are excluded.

## Main files

- `figure_key_locked_claim.png/.pdf`: within-run persistence, coordinate controls,
  same-versus-separate-launch comparison, and retrospective policy frontier.
- `figure_checkpoint_and_oof.png/.pdf`: common-cohort checkpoint curve and outer-LOPO
  probability distribution.
- `key_claim_metrics.csv`: exact AUROCs and protein-clustered intervals.
- `paired_auroc_differences.csv`: paired AUROC differences and clustered intervals.
- `checkpoint_locked_common_cohort.csv`: checkpoint results on the same 52 rows.
- `retrospective_policy_frontier.csv`: every plotted policy threshold.

Intervals resample whole proteins with replacement. They treat each continuous
score or saved OOF probability as fixed and do not refit the complete modeling
procedure within each bootstrap draw. The policy frontier is retrospective and
must not be described as prospective false-stop control. Because the
class-weighted logistic output was not calibrated, its saved `predict_proba`
value is interpreted as a model score rather than an absolute probability.
"""
    (RESULT_DIR / "README.md").write_text(readme)

    print(f"Wrote locked claim figures and result tables to {RESULT_DIR}")
    print(metrics[["metric_id", "n", "auroc", "ci_lo", "ci_hi"]].to_string(index=False))
    print("\nRetrospective p=0.5 operating point:")
    print(operating.to_string(index=False))


if __name__ == "__main__":
    main()
