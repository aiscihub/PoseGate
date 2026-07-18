#!/usr/bin/env python3
"""3 A late-outcome threshold sensitivity and non-obvious-cases robustness checks
for the locked same-trajectory cohort in yau_aj_trajectory_triage.tex.

Both analyses only need the raw continuous late-window RMSD
(`late_pose_rmsd_median_last30pct_A`) and the early per-checkpoint corrected
RMSD (`same_traj_pose_rmsd_{5,10,15,20,30}ns_A` / `same_traj_pose_rmsd_18_20_A`)
already present in the master same-trajectory source-of-truth table. Nothing
here recomputes coordinates from raw trajectories, so this runs anywhere the
CSV is available -- no remote host needed.

Reruns cleanly if the locked cohort ever changes (it should not, by
definition of "locked"), and is a template for repeating both checks on a
future expanded cohort as a secondary stability check, mirroring how the
post-lock model-complexity audit already treats the 73/80-trajectory cohort
as exploratory rather than replacing the n=52 headline.

Usage:
    python3 threshold_and_robustness_sensitivity.py \
        --input path/to/same_trajectory_source_of_truth.csv \
        --output-dir path/to/output/dir
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.stats import rankdata


def auc(score: np.ndarray, y: np.ndarray) -> float:
    y = np.asarray(y)
    score = np.asarray(score)
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    ranks = rankdata(score)
    u = ranks[y == 1].sum() - n1 * (n1 + 1) / 2
    return u / (n1 * n0)


def clustered_ci(
    score: np.ndarray, y: np.ndarray, protein: np.ndarray, draws: int = 5000, seed: int = 7
) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    uniq_proteins = np.unique(protein)
    vals = []
    for _ in range(draws):
        sampled = rng.choice(uniq_proteins, size=len(uniq_proteins), replace=True)
        idx = np.concatenate([np.where(protein == p)[0] for p in sampled])
        yy, ss = y[idx], score[idx]
        if len(set(yy)) < 2:
            continue
        vals.append(auc(ss, yy))
    return tuple(float(x) for x in np.percentile(vals, [2.5, 97.5]))


def load_locked_cohort(input_path: Path):
    rows = list(csv.DictReader(open(input_path)))
    locked = [r for r in rows if r["eligible_locked_primary_pose_model"] == "True"]
    return locked


def threshold_sensitivity(locked, cutoffs=(2.0, 2.5, 3.0, 3.5, 4.0, 5.0), primary=3.0):
    late_raw = np.array([float(r["late_pose_rmsd_median_last30pct_A"]) for r in locked])
    rmsd5 = np.array([float(r["same_traj_pose_rmsd_5ns_A"]) for r in locked])
    rmsd20 = np.array([float(r["same_traj_pose_rmsd_18_20_A"]) for r in locked])
    proteins = np.array([r["protein"] for r in locked])

    # sanity check: recomputing the primary cutoff must exactly reproduce the locked label
    locked_label = np.array([int(r["late_pose_retained_locked"]) for r in locked])
    recomputed = (late_raw < primary).astype(int)
    if not np.array_equal(recomputed, locked_label):
        raise RuntimeError(
            "Recomputed primary-cutoff label does not match late_pose_retained_locked -- "
            "the late_pose_rmsd_median_last30pct_A column or primary cutoff may have changed."
        )

    results = []
    for cutoff in cutoffs:
        label = (late_raw < cutoff).astype(int)
        n_ret, n_non = int(label.sum()), int(len(label) - label.sum())
        a5 = auc(-rmsd5, label)
        a20 = auc(-rmsd20, label)
        ci5 = clustered_ci(-rmsd5, label, proteins)
        ci20 = clustered_ci(-rmsd20, label, proteins)
        results.append(
            {
                "late_threshold_A": cutoff,
                "is_primary": cutoff == primary,
                "n_retained": n_ret,
                "n_non_retained": n_non,
                "auroc_5ns": a5,
                "ci5_lo": ci5[0],
                "ci5_hi": ci5[1],
                "auroc_20ns": a20,
                "ci20_lo": ci20[0],
                "ci20_hi": ci20[1],
            }
        )
    return results


def non_obvious_robustness(locked):
    rmsd20 = np.array([float(r["same_traj_pose_rmsd_18_20_A"]) for r in locked])
    label = np.array([int(r["late_pose_retained_locked"]) for r in locked])
    proteins = np.array([r["protein"] for r in locked])

    subsets = {}
    subsets["full_locked_cohort"] = np.ones(len(label), dtype=bool)

    hi = np.percentile(rmsd20, 90)
    subsets["excluding_top10pct_early_rmsd"] = rmsd20 <= hi

    lo80, hi80 = np.percentile(rmsd20, [10, 90])
    subsets["middle_80pct_early_rmsd"] = (rmsd20 >= lo80) & (rmsd20 <= hi80)

    subsets["near_decision_boundary_pm1.5A"] = np.abs(rmsd20 - 3.0) <= 1.5

    results = []
    for name, mask in subsets.items():
        sub_label, sub_score, sub_protein = label[mask], rmsd20[mask], proteins[mask]
        a = auc(-sub_score, sub_label)
        ci = clustered_ci(-sub_score, sub_label, sub_protein) if len(set(sub_label)) > 1 else (float("nan"),) * 2
        results.append(
            {
                "subset": name,
                "n": int(mask.sum()),
                "n_proteins": len(set(sub_protein)),
                "n_retained": int(sub_label.sum()),
                "auroc_20ns": a,
                "ci_lo": ci[0],
                "ci_hi": ci[1],
            }
        )
    return results


def write_csv(rows, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path}")


def main() -> None:
    here = Path(__file__).resolve()
    repo_root = here.parents[3]
    default_input = (
        repo_root
        / "pipeline/training/outputs/yau_recapture_screen/cad_states/same_trajectory_source_of_truth"
        / "same_trajectory_source_of_truth.csv"
    )
    default_output = (
        repo_root
        / "pipeline/training/outputs/yau_recapture_screen/cad_states/same_trajectory_source_of_truth"
        / "threshold_and_robustness_sensitivity"
    )

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=default_input)
    parser.add_argument("--output-dir", type=Path, default=default_output)
    args = parser.parse_args()

    locked = load_locked_cohort(args.input)
    print(f"locked cohort n={len(locked)}")

    thresh_rows = threshold_sensitivity(locked)
    write_csv(thresh_rows, args.output_dir / "threshold_sensitivity.csv")

    robust_rows = non_obvious_robustness(locked)
    write_csv(robust_rows, args.output_dir / "non_obvious_robustness.csv")

    print("\nThreshold sensitivity:")
    for r in thresh_rows:
        flag = " (primary)" if r["is_primary"] else ""
        print(
            f"  {r['late_threshold_A']:.1f} A{flag}: {r['n_retained']}/{r['n_non_retained']}  "
            f"5ns={r['auroc_5ns']:.3f} [{r['ci5_lo']:.3f},{r['ci5_hi']:.3f}]  "
            f"20ns={r['auroc_20ns']:.3f} [{r['ci20_lo']:.3f},{r['ci20_hi']:.3f}]"
        )
    print("\nNon-obvious-cases robustness:")
    for r in robust_rows:
        print(
            f"  {r['subset']}: n={r['n']} ({r['n_proteins']} proteins), retained={r['n_retained']}, "
            f"AUROC={r['auroc_20ns']:.3f} [{r['ci_lo']:.3f},{r['ci_hi']:.3f}]"
        )


if __name__ == "__main__":
    main()
