#!/usr/bin/env python3
"""PROTOTYPE: cross-ligand coordinate-ablation and early-forecast check.

Tests whether this paper's two central findings generalize beyond
milbemycin, using real 20 ns trajectories for two other ligands
(Beauvericin, Verapamil) recovered from an external drive after the
original raw DCDs were deleted from the compute host (only cached
fingerprint/label CSVs had survived there -- see
`cohort_expansion_80_status.md` for that recovery precedent). These
trajectories are NOT part of the primary 100 ns milbemycin cohort and do
not touch it; this is a standalone cross-ligand corroboration.

Two questions, both testable on the same single-pass measurement:
  1. Coordinate ablation: does protein-relative, PBC-corrected ligand RMSD
     beat PBC-naive and ligand-self-aligned processing for these ligands
     too (mirrors Section "Coordinate-consistent processing reveals the
     signal")?
  2. Short-gap forecasting: does an early (2,7] ns window (matching the
     companion paper's own validated checkpoint) forecast a late (14,20]
     ns window (the final-30%-of-trajectory analog to this paper's
     (70,100] ns rule on the 100 ns cohort)?

Per-frame procedure is identical to script 24's `pose_traces`: protein-CA
Kabsch alignment to frame 0, minimum-image PBC correction of the ligand,
then three pose-measurement variants (corrected / PBC-naive /
ligand-self-aligned) plus centroid displacement, all in one trajectory
pass. Adapted here for 20 ns (1002 frames, 20 ps/frame) instead of 100 ns,
and for a ligand-agnostic root/ligand-name pair instead of the
Milbemycin-only discovery used elsewhere in this pipeline.

Usage (run locally against a mounted external drive; MDAnalysis is not
required on the compute host for this script):
    python 40_cross_ligand_20ns_coordinate_check.py \\
        --root /Volumes/Research/simulation_20ns_md_beau --ligand Beauvericin \\
        --output /Volumes/Research/cross_ligand_analysis/beauvericin_20ns_measurements.csv \\
        --workers 8
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
EARLY_WINDOW_NS = (2.0, 7.0)
LATE_WINDOW_NS = (14.0, 20.0)
NOMINAL_END_NS = 20.0
COMPLETE_TOLERANCE_NS = 0.1
RETENTION_THRESHOLD_A = 3.0  # reused as-is from the milbemycin-calibrated label; not re-derived for these ligands


def discover(root: Path, ligand: str) -> list[tuple[str, str, str, Path, Path]]:
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


def measure(run: tuple[str, str, str, Path, Path], ligand: str) -> dict:
    protein, pocket, replica, topology, trajectory = run
    result = {
        "complex_id": f"{protein}|{ligand}|{pocket}|{replica}",
        "protein": protein,
        "ligand": ligand,
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
        ligand_sel = universe.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
        if ca.n_atoms == 0 or ligand_sel.n_atoms == 0:
            raise ValueError(f"missing selection: CA={ca.n_atoms}, ligand={ligand_sel.n_atoms}")

        dt_ns = float(getattr(universe.trajectory, "dt", 20.0)) / 1000.0
        if not np.isfinite(dt_ns) or dt_ns <= 0:
            dt_ns = 0.02

        times, rmsd_corrected, rmsd_nopbc, rmsd_ligself, disp_corrected = [], [], [], [], []
        ref_ca_centered = ref_lig_rel_pbc = ref_lig_rel_raw = ref_lig_centered = None

        for frame_i, ts in enumerate(universe.trajectory):
            t = (frame_i + 1) * dt_ns
            box = ts.dimensions
            ca_pos = ca.positions
            lig_pos = ligand_sel.positions
            c_t = ca_pos.mean(axis=0)

            lig_rel_raw = lig_pos - c_t
            lig_rel_pbc = minimize_vectors(lig_rel_raw, box=box)
            lig_centroid = lig_pos.mean(axis=0)
            lig_centered = lig_pos - lig_centroid

            if frame_i == 0:
                ref_ca_centered = ca_pos - c_t
                ref_lig_rel_pbc = lig_rel_pbc.copy()
                ref_lig_rel_raw = lig_rel_raw.copy()
                ref_lig_centered = lig_centered.copy()
                times.append(t)
                rmsd_corrected.append(0.0)
                rmsd_nopbc.append(0.0)
                rmsd_ligself.append(0.0)
                disp_corrected.append(0.0)
                continue

            ca_centered = ca_pos - c_t
            R, _ = rotation_matrix(ca_centered, ref_ca_centered)

            aligned_pbc = lig_rel_pbc @ R.T
            rmsd_corrected.append(float(np.sqrt(np.mean(np.sum((aligned_pbc - ref_lig_rel_pbc) ** 2, axis=1)))))
            disp_corrected.append(float(np.linalg.norm(aligned_pbc.mean(axis=0) - ref_lig_rel_pbc.mean(axis=0))))

            aligned_raw = lig_rel_raw @ R.T
            rmsd_nopbc.append(float(np.sqrt(np.mean(np.sum((aligned_raw - ref_lig_rel_raw) ** 2, axis=1)))))

            R_lig, _ = rotation_matrix(lig_centered, ref_lig_centered)
            aligned_lig = lig_centered @ R_lig.T
            rmsd_ligself.append(float(np.sqrt(np.mean(np.sum((aligned_lig - ref_lig_centered) ** 2, axis=1)))))

            times.append(t)

        times = np.array(times)
        rmsd_corrected = np.array(rmsd_corrected)
        last_time_ns = float(times[-1])

        early_mask = (times > EARLY_WINDOW_NS[0]) & (times <= EARLY_WINDOW_NS[1] + 1e-4)
        late_mask = (times > LATE_WINDOW_NS[0]) & (times <= LATE_WINDOW_NS[1] + 1e-4)
        complete = bool(
            last_time_ns >= NOMINAL_END_NS - COMPLETE_TOLERANCE_NS and late_mask.sum() > 0
        )

        result.update(
            status="ok",
            n_frames=len(times),
            frame_interval_ns=dt_ns,
            trajectory_end_ns=last_time_ns,
            complete_20ns_coverage=complete,
            corrected_rmsd_early_2_7_A=float(np.mean(rmsd_corrected[early_mask])) if early_mask.any() else np.nan,
            nopbc_rmsd_early_2_7_A=float(np.mean(np.array(rmsd_nopbc)[early_mask])) if early_mask.any() else np.nan,
            ligself_rmsd_early_2_7_A=float(np.mean(np.array(rmsd_ligself)[early_mask])) if early_mask.any() else np.nan,
            disp_early_2_7_A=float(np.mean(np.array(disp_corrected)[early_mask])) if early_mask.any() else np.nan,
            late_rmsd_median_14_20_A=float(np.median(rmsd_corrected[late_mask])) if late_mask.any() else np.nan,
            late_pose_retained=(
                bool(np.median(rmsd_corrected[late_mask]) < RETENTION_THRESHOLD_A)
                if (complete and late_mask.any())
                else None
            ),
        )
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--ligand", type=str, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    runs = discover(args.root, args.ligand)
    if not runs:
        raise SystemExit(f"No trajectories discovered under {args.root}")
    print(f"Measuring {args.ligand}: {len(runs)} trajectories", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(measure, run, args.ligand): run for run in runs}
        for completed, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            rows.append(row)
            print(f"[{completed}/{len(runs)}] {row['complex_id']} {row['status']}", flush=True)

    table = pd.DataFrame(rows).sort_values(["protein", "pocket", "replica"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output, index=False)
    ok = table[table["status"].eq("ok")]
    complete = ok[ok["complete_20ns_coverage"].fillna(False)]
    print(f"Wrote {args.output}")
    print(
        f"Measured {len(ok)}/{len(table)} trajectories; "
        f"complete 20 ns coverage: {len(complete)} across {complete['protein'].nunique()} proteins; "
        f"retained: {int(complete['late_pose_retained'].sum())}/{len(complete)}"
    )


if __name__ == "__main__":
    main()
