#!/usr/bin/env python3
"""Generate early and late corrected-pose measurements for every raw 100 ns run.

Unlike the locked source-of-truth builder, this script discovers trajectories
directly from the raw simulation tree. It is intended to bring newly completed
runs into an auditable inventory before any headline cohort is re-locked.
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


CHECKPOINTS_NS = (5, 10, 15, 20, 30)
WINDOW_NS = 2.0
LATE_WINDOW_NS = (70.0, 100.0)
LIGAND_RESNAME = "UNK"
COMPLETE_TOLERANCE_NS = 0.05


def discover(root: Path) -> list[tuple[str, str, str, Path, Path]]:
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


def measure(run: tuple[str, str, str, Path, Path]) -> dict:
    protein, pocket, replica, topology, trajectory = run
    result = {
        "complex_id": f"{protein}|Milbemycin|{pocket}|{replica}",
        "protein": protein,
        "ligand": "Milbemycin",
        "pocket": pocket,
        "replica": replica,
        "topology_path": str(topology),
        "trajectory_path": str(trajectory),
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

        rmsd_windows = {checkpoint: [] for checkpoint in CHECKPOINTS_NS}
        displacement_windows = {checkpoint: [] for checkpoint in CHECKPOINTS_NS}
        late_rmsd = []
        ref_ca_centered = None
        ref_ligand_rel = None
        last_time_ns = np.nan

        for frame_i, ts in enumerate(universe.trajectory):
            time_ns = (frame_i + 1) * dt_ns
            last_time_ns = time_ns
            box = ts.dimensions
            if box is None or len(box) < 3 or np.any(np.asarray(box[:3]) <= 0):
                raise ValueError(f"invalid periodic box at frame {frame_i}")

            ca_positions = ca.positions
            protein_center = ca_positions.mean(axis=0)
            ligand_rel = minimize_vectors(ligand.positions - protein_center, box=box)
            ca_centered = ca_positions - protein_center

            if frame_i == 0:
                ref_ca_centered = ca_centered.copy()
                ref_ligand_rel = ligand_rel.copy()
                continue

            rotation, _ = rotation_matrix(ca_centered, ref_ca_centered)
            ligand_aligned = ligand_rel @ rotation.T
            delta = ligand_aligned - ref_ligand_rel
            pose_rmsd = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
            centroid_displacement = float(np.linalg.norm(delta.mean(axis=0)))

            for checkpoint in CHECKPOINTS_NS:
                if checkpoint - WINDOW_NS < time_ns <= checkpoint + 1e-4:
                    rmsd_windows[checkpoint].append(pose_rmsd)
                    displacement_windows[checkpoint].append(centroid_displacement)
            if LATE_WINDOW_NS[0] < time_ns <= LATE_WINDOW_NS[1] + 1e-4:
                late_rmsd.append(pose_rmsd)

        result.update(
            status="ok",
            n_frames=len(universe.trajectory),
            frame_interval_ns=dt_ns,
            trajectory_end_ns=last_time_ns,
            complete_70_100_coverage=bool(
                np.isfinite(last_time_ns)
                and last_time_ns >= LATE_WINDOW_NS[1] - COMPLETE_TOLERANCE_NS
                and len(late_rmsd) > 0
            ),
            late_pose_rmsd_median_70_100_A=(
                float(np.median(late_rmsd)) if late_rmsd else np.nan
            ),
        )
        for checkpoint in CHECKPOINTS_NS:
            result[f"pose_rmsd_mean_{checkpoint - 2}_{checkpoint}_A"] = (
                float(np.mean(rmsd_windows[checkpoint]))
                if rmsd_windows[checkpoint]
                else np.nan
            )
            result[f"centroid_displacement_mean_{checkpoint - 2}_{checkpoint}_A"] = (
                float(np.mean(displacement_windows[checkpoint]))
                if displacement_windows[checkpoint]
                else np.nan
            )
            result[f"n_frames_{checkpoint - 2}_{checkpoint}"] = len(rmsd_windows[checkpoint])
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    runs = discover(args.root)
    if not runs:
        raise SystemExit(f"No trajectories discovered under {args.root}")
    print(f"Discovered {len(runs)} trajectories", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(measure, run): run for run in runs}
        for completed, future in enumerate(as_completed(futures), start=1):
            row = future.result()
            rows.append(row)
            print(
                f"[{completed}/{len(runs)}] {row['complex_id']} {row['status']}",
                flush=True,
            )

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
