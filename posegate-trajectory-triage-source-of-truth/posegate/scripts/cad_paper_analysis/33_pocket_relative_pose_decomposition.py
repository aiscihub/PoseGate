#!/usr/bin/env python3
"""PROTOTYPE: Pocket-Relative Pose Decomposition (PRPD).

Not part of the locked or expanded source-of-truth pipelines. Tests
whether splitting the validated corrected-RMSD signal into translation,
rotation, and internal-deformation components is (a) an interpretation
tool that shows which physical motion carries the signal, or (b) actually
improves the forecast over the single coherent RMSD number.

This intentionally reuses the SAME whole-protein Cα Kabsch frame as the
validated method (script 24/28) -- no local pocket-frame or medoid
baseline yet. That keeps this test isolated to one question: does
decomposing the motion help, independent of the separate (and separately
testable) question of whether a different reference frame helps.

Per-frame, relative to frame 0, after protein-Cα Kabsch alignment +
minimum-image PBC correction of the ligand into the protein frame
(identical procedure to script 24's `rmsd_corrected` / `disp_corrected`):

  T(t): ligand heavy-atom centroid displacement in the protein frame.
        Identical quantity to `disp_corrected` in script 24/28.

  Theta(t): rotation angle relating the ligand's current orientation to
        its frame-0 orientation, from a SEPARATE ligand-only Kabsch fit
        (center ligand on its own centroid, fit rotation R_lig onto the
        frame-0 ligand shape -- same R_lig script 24 already computes for
        `rmsd_ligself`). Theta(t) = arccos(clip((trace(R_lig) - 1) / 2, -1, 1)),
        in radians.

  I(t): ligand internal-deformation RMSD after removing translation and
        rotation via that same ligand-only Kabsch fit (apply R_lig, then
        RMSD against the frame-0 ligand shape). Identical quantity to
        `rmsd_ligself` in script 24/28.

Baseline-normalized combined score, per the proposal: within each
trajectory's own (0.2, 1.0] ns baseline window, compute the median
absolute deviation (MAD) of T, Theta, I as that trajectory's own
fluctuation scale (s_T, s_Theta, s_I), then per frame:

  Z_pose(t) = sqrt( (T(t)/(s_T+eps))^2 + (Theta(t)/(s_Theta+eps))^2 + (I(t)/(s_I+eps))^2 )

Both raw (T, Theta, I) and Z-normalized (Z_T, Z_Theta, Z_I, Z_pose) are
window-averaged into the same trailing-2ns checkpoints and (70,100] ns
late window used by the locked pipeline, so every column here lines up
with `same_traj_pose_rmsd_18_20_A` / `corrected_pose_rmsd_mean_18_20_A`
for direct comparison.

Usage (against the 73-row expanded cohort):
    python 33_pocket_relative_pose_decomposition.py \\
        --root /media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md \\
        --output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/prpd_candidate_measurements.csv \\
        --complex-ids ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/master/expanded_trajectory_source_of_truth.csv \\
        --workers 10

Evaluate exactly as with CW-PGD: merge onto
`expanded_trajectory_source_of_truth.csv` by complex_id and run
`metric_row` / `paired_cluster_bootstrap_difference` from
`26_plot_locked_claim_results.py` against `late_pose_retained`, sign=-1.0
for every column here (higher T/Theta/I/Z_pose = more motion = more
likely non-retained).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import minimize_vectors
import numpy as np
import pandas as pd


LIGAND_RESNAME = "UNK"
BASELINE_WINDOW_NS = (0.2, 1.0)
MIN_BASELINE_FRAMES = 2
CHECKPOINTS_NS = (5, 10, 15, 20, 30)
WINDOW_NS = 2.0
LATE_WINDOW_NS = (70.0, 100.0)
COMPLETE_TOLERANCE_NS = 0.05
EPS_T_A = 1e-3
EPS_THETA_RAD = 1e-3
EPS_I_A = 1e-3


def discover(root: Path) -> list[tuple[str, str, str, Path, Path]]:
    """Identical to 28_generate_full_100ns_measurements.discover()."""
    runs = []
    for trajectory in sorted(root.glob("*/simulation_explicit/*/*/*_trajectory.dcd")):
        replica_dir = trajectory.parent
        replica = replica_dir.name
        pocket = replica_dir.parent.name
        protein = replica_dir.parents[2].name
        topology_candidates = sorted(
            path
            for path in replica_dir.glob("*_initial_frame.pdb")
            if not path.name.endswith("_stripped_initial_frame.pdb")
        )
        if len(topology_candidates) != 1:
            continue
        runs.append((protein, pocket, replica, topology_candidates[0], trajectory))
    return runs


def mad(values: np.ndarray) -> float:
    med = float(np.median(values))
    return float(np.median(np.abs(values - med)))


def measure(run: tuple[str, str, str, Path, Path]) -> dict:
    protein, pocket, replica, topology, trajectory = run
    result = {
        "complex_id": f"{protein}|Milbemycin|{pocket}|{replica}",
        "protein": protein,
        "ligand": "Milbemycin",
        "pocket": pocket,
        "replica": replica,
        "status": "error",
        "error": "",
    }
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            universe = mda.Universe(str(topology), str(trajectory))

        ca = universe.select_atoms("protein and name CA")
        ligand = universe.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
        if ca.n_atoms == 0 or ligand.n_atoms == 0:
            raise ValueError(f"missing selection: CA={ca.n_atoms}, ligand={ligand.n_atoms}")

        dt_ns = float(getattr(universe.trajectory, "dt", 200.0)) / 1000.0
        if not np.isfinite(dt_ns) or dt_ns <= 0:
            dt_ns = 0.2

        times, t_vals, theta_vals, i_vals = [], [], [], []
        ref_ca_centered = ref_ligand_rel = ref_ligand_centered = None

        for frame_i, ts in enumerate(universe.trajectory):
            time_ns = (frame_i + 1) * dt_ns
            box = ts.dimensions
            ca_positions = ca.positions
            ligand_positions = ligand.positions
            protein_center = ca_positions.mean(axis=0)
            ca_centered = ca_positions - protein_center
            ligand_rel = minimize_vectors(ligand_positions - protein_center, box=box)
            ligand_centroid = ligand_positions.mean(axis=0)
            ligand_centered = ligand_positions - ligand_centroid

            if frame_i == 0:
                ref_ca_centered = ca_centered.copy()
                ref_ligand_rel = ligand_rel.copy()
                ref_ligand_centered = ligand_centered.copy()
                times.append(time_ns)
                t_vals.append(0.0)
                theta_vals.append(0.0)
                i_vals.append(0.0)
                continue

            # Translation: identical procedure/quantity to disp_corrected.
            protein_rotation, _ = rotation_matrix(ca_centered, ref_ca_centered)
            ligand_aligned = ligand_rel @ protein_rotation.T
            t_val = float(
                np.linalg.norm(ligand_aligned.mean(axis=0) - ref_ligand_rel.mean(axis=0))
            )

            # Rotation + internal deformation: separate ligand-only Kabsch fit.
            ligand_rotation, _ = rotation_matrix(ligand_centered, ref_ligand_centered)
            trace = float(np.trace(ligand_rotation))
            theta_val = float(np.arccos(np.clip((trace - 1.0) / 2.0, -1.0, 1.0)))
            ligand_self_aligned = ligand_centered @ ligand_rotation.T
            i_val = float(
                np.sqrt(np.mean(np.sum((ligand_self_aligned - ref_ligand_centered) ** 2, axis=1)))
            )

            times.append(time_ns)
            t_vals.append(t_val)
            theta_vals.append(theta_val)
            i_vals.append(i_val)

        times = np.array(times)
        t_vals = np.array(t_vals)
        theta_vals = np.array(theta_vals)
        i_vals = np.array(i_vals)

        baseline_mask = (times > BASELINE_WINDOW_NS[0]) & (times <= BASELINE_WINDOW_NS[1] + 1e-4)
        n_baseline_frames = int(baseline_mask.sum())
        if n_baseline_frames < MIN_BASELINE_FRAMES:
            raise ValueError(f"insufficient baseline frames: {n_baseline_frames}")

        s_t = mad(t_vals[baseline_mask])
        s_theta = mad(theta_vals[baseline_mask])
        s_i = mad(i_vals[baseline_mask])

        z_t = t_vals / (s_t + EPS_T_A)
        z_theta = theta_vals / (s_theta + EPS_THETA_RAD)
        z_i = i_vals / (s_i + EPS_I_A)
        z_pose = np.sqrt(z_t ** 2 + z_theta ** 2 + z_i ** 2)

        last_time_ns = float(times[-1])
        late_mask = (times > LATE_WINDOW_NS[0]) & (times <= LATE_WINDOW_NS[1] + 1e-4)
        complete = bool(
            last_time_ns >= LATE_WINDOW_NS[1] - COMPLETE_TOLERANCE_NS and late_mask.sum() > 0
        )

        result.update(
            status="ok",
            n_frames=len(universe.trajectory),
            frame_interval_ns=dt_ns,
            trajectory_end_ns=last_time_ns,
            n_baseline_frames=n_baseline_frames,
            baseline_mad_T_A=s_t,
            baseline_mad_Theta_rad=s_theta,
            baseline_mad_I_A=s_i,
            complete_70_100_coverage=complete,
            late_T_median_70_100_A=float(np.median(t_vals[late_mask])) if late_mask.any() else np.nan,
            late_Theta_median_70_100_rad=float(np.median(theta_vals[late_mask])) if late_mask.any() else np.nan,
            late_I_median_70_100_A=float(np.median(i_vals[late_mask])) if late_mask.any() else np.nan,
            late_Zpose_median_70_100=float(np.median(z_pose[late_mask])) if late_mask.any() else np.nan,
        )

        for checkpoint in CHECKPOINTS_NS:
            window_mask = (times > checkpoint - WINDOW_NS) & (times <= checkpoint + 1e-4)
            n_window = int(window_mask.sum())
            suffix = f"{checkpoint - WINDOW_NS:g}_{checkpoint:g}"
            for name, series in (("T", t_vals), ("Theta", theta_vals), ("I", i_vals), ("Zpose", z_pose)):
                result[f"prpd_{name}_mean_{suffix}"] = (
                    float(np.mean(series[window_mask])) if n_window else np.nan
                )
            result[f"n_frames_{suffix}"] = n_window
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--complex-ids",
        type=Path,
        default=None,
        help=(
            "Optional CSV with a complex_id column to restrict the run to "
            "(e.g. expanded_trajectory_source_of_truth.csv). If it has an "
            "eligible_primary_cohort column, only True rows are kept."
        ),
    )
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    runs = discover(args.root)
    if not runs:
        raise SystemExit(f"No trajectories discovered under {args.root}")

    if args.complex_ids is not None:
        filter_table = pd.read_csv(args.complex_ids)
        if "eligible_primary_cohort" in filter_table.columns:
            filter_table = filter_table[filter_table["eligible_primary_cohort"].astype(bool)]
        keep = set(filter_table["complex_id"])
        runs = [run for run in runs if f"{run[0]}|Milbemycin|{run[1]}|{run[2]}" in keep]
        if not runs:
            raise SystemExit("No discovered trajectories matched --complex-ids")

    print(f"Measuring PRPD for {len(runs)} trajectories", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(measure, run): run for run in runs}
        for completed, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            rows.append(row)
            print(f"[{completed}/{len(runs)}] {row['complex_id']} {row['status']}", flush=True)

    table = pd.DataFrame(rows).sort_values(["protein", "pocket", "replica"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    ok = table[table["status"].eq("ok")]
    complete = ok[ok["complete_70_100_coverage"].fillna(False)]
    print(f"Wrote {args.output}")
    print(
        f"Measured {len(ok)}/{len(table)} trajectories; "
        f"complete 100 ns: {len(complete)} across {complete['protein'].nunique()} proteins"
    )


if __name__ == "__main__":
    main()
