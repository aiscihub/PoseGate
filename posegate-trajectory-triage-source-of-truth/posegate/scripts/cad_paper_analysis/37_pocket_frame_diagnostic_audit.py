#!/usr/bin/env python3
"""Full-cohort diagnostic audit for the local-pocket-alignment control.

Covers checks 1, 4, and 5 of the requested verification battery (check 2,
agreement statistics, and check 3, the synthetic unit test, are separate:
`36_local_pocket_control_evaluation.py`'s companion agreement script and
`test_synthetic_frame_relative_theta.py` respectively):

  1. Full 73-trajectory frame-divergence audit: per-frame angular
     difference between the whole-protein and pocket-only (8 A) fitted
     Kabsch rotations, aggregated to per-trajectory median/IQR/max, joined
     with pocket Cα count and fraction of the full protein.
  4. Window-averaging check: the same per-frame T/Theta traces are saved
     so agreement can be compared at the individual-frame level, a single
     nearest-checkpoint-frame level, and the 2 ns window-averaged level
     (the last of these duplicates script 35's checkpoint columns, saved
     here again for a self-contained per-trajectory comparison).
  5. Geometric conditioning: for the frozen reference (frame 0) pocket and
     whole-protein point clouds, the eigenvalues of their centered spatial
     covariance matrix (largest/mid/smallest) and the ratio of smallest to
     largest -- a near-zero ratio flags a poorly conditioned (nearly
     planar or collinear) point set, which would make the Kabsch fit
     numerically unstable.

Writes two files:
  - per-frame trace (long format, ~73*500 rows): time_ns, T_whole,
    Theta_whole, T_pocket, Theta_pocket, angular_divergence_deg,
    kabsch_residual_whole_A, kabsch_residual_pocket_A
  - per-trajectory summary (one row each): pocket CA count/fraction,
    covariance eigenvalues and conditioning ratio for both point sets,
    median/IQR/max angular divergence, median/max Kabsch residuals

Usage:
    python 37_pocket_frame_diagnostic_audit.py \\
        --root /media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md \\
        --frames-output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/pocket_frame_diagnostic_traces.csv \\
        --summary-output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/pocket_frame_diagnostic_summary.csv \\
        --complex-ids ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/master/expanded_trajectory_source_of_truth.csv \\
        --workers 10
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import sys
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import distance_array, minimize_vectors
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

control_mod = importlib.import_module("35_local_pocket_alignment_control")
discover = control_mod.discover
frame_relative_theta = control_mod.frame_relative_theta

LIGAND_RESNAME = "UNK"
PRIMARY_RADIUS_A = 8.0
MIN_POCKET_CA = 4
CHECKPOINTS_NS = (5, 10, 20)
WINDOW_NS = 2.0


def kabsch_residual(current_centered: np.ndarray, ref_centered: np.ndarray) -> float:
    q, _ = rotation_matrix(current_centered, ref_centered)
    aligned = current_centered @ q.T
    return float(np.sqrt(np.mean(np.sum((aligned - ref_centered) ** 2, axis=1))))


def conditioning(points_centered: np.ndarray) -> tuple[float, float, float, float]:
    cov = np.cov(points_centered.T)
    eigenvalues = np.sort(np.linalg.eigvalsh(cov))[::-1]  # descending: lambda1 >= lambda2 >= lambda3
    ratio = float(eigenvalues[2] / eigenvalues[0]) if eigenvalues[0] > 0 else float("nan")
    return float(eigenvalues[0]), float(eigenvalues[1]), float(eigenvalues[2]), ratio


def measure(run: tuple[str, str, str, Path, Path]) -> tuple[dict, list[dict]]:
    protein, pocket, replica, topology, trajectory = run
    complex_id = f"{protein}|Milbemycin|{pocket}|{replica}"
    summary = {
        "complex_id": complex_id, "protein": protein, "pocket": pocket, "replica": replica,
        "status": "error", "error": "",
    }
    frame_rows: list[dict] = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            universe = mda.Universe(str(topology), str(trajectory))

        ca = universe.select_atoms("protein and name CA")
        protein_heavy = universe.select_atoms("protein and not name H*")
        ligand = universe.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
        if ca.n_atoms == 0 or ligand.n_atoms == 0 or protein_heavy.n_atoms == 0:
            raise ValueError("missing selection")

        dt_ns = float(getattr(universe.trajectory, "dt", 200.0)) / 1000.0
        if not np.isfinite(dt_ns) or dt_ns <= 0:
            dt_ns = 0.2

        universe.trajectory[0]
        box0 = universe.trajectory.ts.dimensions
        d0 = distance_array(protein_heavy.positions, ligand.positions, box=box0)
        min_dist = d0.min(axis=1)
        pocket_resindex_set = set(np.unique(protein_heavy.resindices[min_dist <= PRIMARY_RADIUS_A]))
        mask = np.isin(ca.resindices, list(pocket_resindex_set))
        n_pocket_ca = int(mask.sum())
        if n_pocket_ca < MIN_POCKET_CA:
            raise ValueError(f"insufficient pocket CA atoms (n={n_pocket_ca})")
        pocket_ca = ca[mask]

        ref_ca_pos0 = ca.positions.copy()
        ref_pocket_pos0 = pocket_ca.positions.copy()
        ref_ca_centered0 = ref_ca_pos0 - ref_ca_pos0.mean(axis=0)
        ref_pocket_centered0 = ref_pocket_pos0 - ref_pocket_pos0.mean(axis=0)

        eig_whole = conditioning(ref_ca_centered0)
        eig_pocket = conditioning(ref_pocket_centered0)
        summary.update(
            n_pocket_ca=n_pocket_ca,
            n_total_ca=int(ca.n_atoms),
            pocket_fraction=n_pocket_ca / ca.n_atoms,
            whole_eig1=eig_whole[0], whole_eig2=eig_whole[1], whole_eig3=eig_whole[2], whole_conditioning_ratio=eig_whole[3],
            pocket_eig1=eig_pocket[0], pocket_eig2=eig_pocket[1], pocket_eig3=eig_pocket[2], pocket_conditioning_ratio=eig_pocket[3],
        )

        ref_ca_centered = ref_ligand_rel_whole = None
        ref_pocket_centered = ref_ligand_rel_pocket = None

        universe.trajectory.rewind()
        for frame_i, ts in enumerate(universe.trajectory):
            time_ns = (frame_i + 1) * dt_ns
            box = ts.dimensions
            ca_pos = ca.positions
            ligand_pos = ligand.positions
            pocket_pos = pocket_ca.positions

            protein_center = ca_pos.mean(axis=0)
            ca_centered = ca_pos - protein_center
            ligand_rel_whole = minimize_vectors(ligand_pos - protein_center, box=box)

            pocket_center = pocket_pos.mean(axis=0)
            pocket_centered = pocket_pos - pocket_center
            ligand_rel_pocket = minimize_vectors(ligand_pos - pocket_center, box=box)

            if frame_i == 0:
                ref_ca_centered = ca_centered.copy()
                ref_ligand_rel_whole = ligand_rel_whole.copy()
                ref_pocket_centered = pocket_centered.copy()
                ref_ligand_rel_pocket = ligand_rel_pocket.copy()
                frame_rows.append(
                    dict(complex_id=complex_id, protein=protein, time_ns=time_ns,
                         T_whole=0.0, Theta_whole=0.0, T_pocket=0.0, Theta_pocket=0.0,
                         angular_divergence_deg=0.0, kabsch_residual_whole_A=0.0, kabsch_residual_pocket_A=0.0)
                )
                continue

            t_whole, theta_whole, _ = frame_relative_theta(
                ca_centered, ref_ca_centered, ligand_rel_whole, ref_ligand_rel_whole
            )
            t_pocket, theta_pocket, _ = frame_relative_theta(
                pocket_centered, ref_pocket_centered, ligand_rel_pocket, ref_ligand_rel_pocket
            )

            q_whole, _ = rotation_matrix(ca_centered, ref_ca_centered)
            q_pocket, _ = rotation_matrix(pocket_centered, ref_pocket_centered)
            delta = q_pocket @ q_whole.T
            angular_divergence = float(
                np.degrees(np.arccos(np.clip((np.trace(delta) - 1.0) / 2.0, -1.0, 1.0)))
            )
            residual_whole = kabsch_residual(ca_centered, ref_ca_centered)
            residual_pocket = kabsch_residual(pocket_centered, ref_pocket_centered)

            frame_rows.append(
                dict(
                    complex_id=complex_id, protein=protein, time_ns=time_ns,
                    T_whole=t_whole, Theta_whole=np.degrees(theta_whole),
                    T_pocket=t_pocket, Theta_pocket=np.degrees(theta_pocket),
                    angular_divergence_deg=angular_divergence,
                    kabsch_residual_whole_A=residual_whole,
                    kabsch_residual_pocket_A=residual_pocket,
                )
            )

        divergence = np.array([row["angular_divergence_deg"] for row in frame_rows[1:]])
        residual_pocket_series = np.array([row["kabsch_residual_pocket_A"] for row in frame_rows[1:]])
        summary.update(
            status="ok",
            n_frames=len(frame_rows),
            divergence_median_deg=float(np.median(divergence)),
            divergence_q25_deg=float(np.percentile(divergence, 25)),
            divergence_q75_deg=float(np.percentile(divergence, 75)),
            divergence_max_deg=float(np.max(divergence)),
            kabsch_residual_pocket_median_A=float(np.median(residual_pocket_series)),
            kabsch_residual_pocket_max_A=float(np.max(residual_pocket_series)),
        )
    except Exception as exc:
        summary["error"] = f"{type(exc).__name__}: {exc}"
        frame_rows = []
    return summary, frame_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--frames-output", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--complex-ids", type=Path, default=None)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    runs = discover(args.root)
    if args.complex_ids is not None:
        filter_table = pd.read_csv(args.complex_ids)
        if "eligible_primary_cohort" in filter_table.columns:
            filter_table = filter_table[filter_table["eligible_primary_cohort"].astype(bool)]
        keep = set(filter_table["complex_id"])
        runs = [run for run in runs if f"{run[0]}|Milbemycin|{run[1]}|{run[2]}" in keep]

    print(f"Auditing {len(runs)} trajectories", flush=True)

    summaries, all_frame_rows = [], []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(measure, run): run for run in runs}
        for completed, future in enumerate(as_completed(futures), start=1):
            summary, frame_rows = future.result()
            summaries.append(summary)
            all_frame_rows.extend(frame_rows)
            print(f"[{completed}/{len(runs)}] {summary['complex_id']} {summary['status']}", flush=True)

    summary_table = pd.DataFrame(summaries).sort_values(["protein", "pocket", "replica"])
    args.summary_output.parent.mkdir(parents=True, exist_ok=True)
    summary_table.to_csv(args.summary_output, index=False)
    print(f"Wrote {args.summary_output}")

    frames_table = pd.DataFrame(all_frame_rows)
    frames_table.to_csv(args.frames_output, index=False)
    print(f"Wrote {args.frames_output} ({len(frames_table)} rows)")

    ok = summary_table[summary_table["status"].eq("ok")]
    print(f"\n{len(ok)}/{len(summary_table)} trajectories succeeded")
    print("Pocket fraction of total protein: median={:.3f}  range=[{:.3f}, {:.3f}]".format(
        ok["pocket_fraction"].median(), ok["pocket_fraction"].min(), ok["pocket_fraction"].max()
    ))
    print("Angular divergence (median per trajectory): overall median={:.2f} deg  max={:.2f} deg".format(
        ok["divergence_median_deg"].median(), ok["divergence_max_deg"].max()
    ))
    print("Pocket conditioning ratio (lambda3/lambda1): median={:.4f}  min={:.4f}".format(
        ok["pocket_conditioning_ratio"].median(), ok["pocket_conditioning_ratio"].min()
    ))


if __name__ == "__main__":
    main()
