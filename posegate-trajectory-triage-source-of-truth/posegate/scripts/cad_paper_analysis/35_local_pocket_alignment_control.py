#!/usr/bin/env python3
"""PROTOTYPE: local-pocket-alignment control for PRPD's rotation component.

The one remaining analysis before interpreting Theta (script 33/34) as
"ligand reorientation inside the binding pocket": does Theta survive when
measured relative to a frozen LOCAL pocket frame instead of the whole
protein's Cα frame? If it collapses, the whole-frame signal was partly or
wholly transporter-domain motion, not ligand reorientation.

CORRECTNESS FIX FOLDED IN (not a new feature -- required for this control
to test the intended hypothesis): script 33 computed Theta via a ligand-
only Kabsch fit in raw lab coordinates, never referencing the protein at
all. That does not measure "rotation relative to the protein frame" -- it
measures the ligand's absolute rotation in the simulation box, which is
contaminated by whatever rigid-body tumbling the whole complex undergoes
in an unrestrained NPT box. Theta here is instead computed properly
relative to each reference frame:
  1. Fit the reference frame's own alignment rotation Q_frame(t) (whole-
     protein Cα, or pocket-only Cα) exactly as for T(t).
  2. Apply Q_frame(t) to the (PBC-corrected, frame-centroid-relative)
     ligand coordinates -- this expresses the ligand's current position
     AS SEEN FROM the frame-0 reference orientation, removing the
     reference's own rotation.
  3. Fit a second, ligand-only Kabsch between that frame-corrected ligand
     shape and the frame-0 ligand shape, centered on the ligand's own
     centroid. Theta_frame(t) = arccos(clip((trace(Q_lig)-1)/2, -1, 1))
     from that second fit.
Internal deformation I(t) is unaffected by this fix -- it was already a
pure ligand-self-superposition with no frame reference, which is correct
for I by definition (present in every model as a fixed identity, kept
here only as a sanity check that it stays ~unchanged between arms).

Pocket definition (frozen from the reference/frame-0 structure only, no
outcome dependence, never updated over the trajectory):
  P_0 = { residue i : min_{a in i, b in ligand} ||r_a - r_b|| <= radius,
          measured over protein heavy atoms in the reference frame }
Alignment atoms from P_0 are Cα only (matches the existing whole-protein
convention; the user's own spec allows this simplification explicitly, so
the comparison isolates ONE variable -- spatial extent of the alignment
atom set -- rather than conflating it with a second variable, atom type).
Radius 8 A is primary; 6 and 10 A are sensitivity checks. The radius is
NOT chosen by AUROC.

PBC handling: identical minimum-image (`minimize_vectors`) convention used
throughout the rest of this pipeline (scripts 24/28/32/33), applied
consistently to both the whole-protein and pocket arms so the comparison
is controlled. Full bond-connectivity molecule-whole reconstruction was
NOT implemented (that is a separate, lower-priority correctness item);
this is a scoping choice, not an oversight -- flagged explicitly because
milbemycin is small relative to the >=1 nm box padding and stays pocket-
bound (not edge-adjacent) through the runs used here, so the existing
per-atom minimum-image correction should already be adequate for this
cohort.

Usage:
    python 35_local_pocket_alignment_control.py \\
        --root /media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md \\
        --output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/local_pocket_control_measurements.csv \\
        --complex-ids ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/master/expanded_trajectory_source_of_truth.csv \\
        --workers 10
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import warnings

import MDAnalysis as mda
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import distance_array, minimize_vectors
import numpy as np
import pandas as pd


LIGAND_RESNAME = "UNK"
POCKET_RADII_A = (6.0, 8.0, 10.0)
PRIMARY_RADIUS_A = 8.0
MIN_POCKET_CA = 4
CHECKPOINTS_NS = (5, 10, 20)
WINDOW_NS = 2.0
LATE_WINDOW_NS = (70.0, 100.0)
COMPLETE_TOLERANCE_NS = 0.05


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


def frame_relative_theta(
    frame_centered: np.ndarray,
    ref_frame_centered: np.ndarray,
    ligand_rel: np.ndarray,
    ref_ligand_rel: np.ndarray,
) -> tuple[float, float, np.ndarray]:
    """T and Theta of the ligand relative to one reference frame (whole or pocket).

    `frame_centered` / `ref_frame_centered`: current/reference alignment-atom
    coordinates, centered on their own centroid.
    `ligand_rel` / `ref_ligand_rel`: current/reference ligand coordinates,
    PBC-corrected relative to the alignment-atom centroid (NOT yet rotated).
    Returns (T, Theta, ligand_aligned) -- ligand_aligned is returned so the
    caller can reuse it if needed.
    """
    q_frame, _ = rotation_matrix(frame_centered, ref_frame_centered)
    ligand_aligned = ligand_rel @ q_frame.T
    t_val = float(np.linalg.norm(ligand_aligned.mean(axis=0) - ref_ligand_rel.mean(axis=0)))

    cur_lig_centered = ligand_aligned - ligand_aligned.mean(axis=0)
    ref_lig_centered = ref_ligand_rel - ref_ligand_rel.mean(axis=0)
    q_lig, _ = rotation_matrix(cur_lig_centered, ref_lig_centered)
    theta_val = float(np.arccos(np.clip((np.trace(q_lig) - 1.0) / 2.0, -1.0, 1.0)))
    return t_val, theta_val, ligand_aligned


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
        protein_heavy = universe.select_atoms("protein and not name H*")
        ligand = universe.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
        if ca.n_atoms == 0 or ligand.n_atoms == 0 or protein_heavy.n_atoms == 0:
            raise ValueError(
                f"missing selection: CA={ca.n_atoms}, protein_heavy={protein_heavy.n_atoms}, ligand={ligand.n_atoms}"
            )

        dt_ns = float(getattr(universe.trajectory, "dt", 200.0)) / 1000.0
        if not np.isfinite(dt_ns) or dt_ns <= 0:
            dt_ns = 0.2

        # --- Frame 0: freeze the pocket residue set at each radius ---
        universe.trajectory[0]
        box0 = universe.trajectory.ts.dimensions
        d0 = distance_array(protein_heavy.positions, ligand.positions, box=box0)
        min_dist_per_heavy_atom = d0.min(axis=1)
        heavy_resindices = protein_heavy.resindices
        ca_resindices = ca.resindices

        pocket_ca_by_radius: dict[float, object] = {}
        n_pocket_residues = {}
        for radius in POCKET_RADII_A:
            pocket_resindex_set = set(np.unique(heavy_resindices[min_dist_per_heavy_atom <= radius]))
            mask = np.isin(ca_resindices, list(pocket_resindex_set))
            n_pocket_residues[radius] = int(mask.sum())
            pocket_ca_by_radius[radius] = ca[mask] if mask.sum() >= MIN_POCKET_CA else None

        if pocket_ca_by_radius[PRIMARY_RADIUS_A] is None:
            raise ValueError(
                f"insufficient pocket CA atoms at primary radius {PRIMARY_RADIUS_A} A "
                f"(n={n_pocket_residues[PRIMARY_RADIUS_A]})"
            )

        # --- Single full-trajectory pass ---
        times: list[float] = []
        t_whole_vals: list[float] = []
        theta_whole_vals: list[float] = []
        i_vals: list[float] = []
        t_pocket_vals: dict[float, list[float]] = {r: [] for r in POCKET_RADII_A}
        theta_pocket_vals: dict[float, list[float]] = {r: [] for r in POCKET_RADII_A}

        ref_ca_centered = ref_ligand_rel_whole = ref_ligand_self_centered = None
        ref_pocket_centered: dict[float, np.ndarray | None] = {r: None for r in POCKET_RADII_A}
        ref_ligand_rel_pocket: dict[float, np.ndarray | None] = {r: None for r in POCKET_RADII_A}

        universe.trajectory.rewind()
        for frame_i, ts in enumerate(universe.trajectory):
            time_ns = (frame_i + 1) * dt_ns
            box = ts.dimensions
            ca_pos = ca.positions
            ligand_pos = ligand.positions

            protein_center = ca_pos.mean(axis=0)
            ca_centered = ca_pos - protein_center
            ligand_rel_whole = minimize_vectors(ligand_pos - protein_center, box=box)

            ligand_centroid = ligand_pos.mean(axis=0)
            ligand_self_centered = ligand_pos - ligand_centroid

            pocket_frame_data = {}
            for radius in POCKET_RADII_A:
                pocket_ca = pocket_ca_by_radius[radius]
                if pocket_ca is None:
                    pocket_frame_data[radius] = None
                    continue
                pocket_pos = pocket_ca.positions
                pocket_center = pocket_pos.mean(axis=0)
                pocket_centered = pocket_pos - pocket_center
                ligand_rel_pocket = minimize_vectors(ligand_pos - pocket_center, box=box)
                pocket_frame_data[radius] = (pocket_centered, ligand_rel_pocket)

            if frame_i == 0:
                ref_ca_centered = ca_centered.copy()
                ref_ligand_rel_whole = ligand_rel_whole.copy()
                ref_ligand_self_centered = ligand_self_centered.copy()
                for radius in POCKET_RADII_A:
                    if pocket_frame_data[radius] is not None:
                        pocket_centered, ligand_rel_pocket = pocket_frame_data[radius]
                        ref_pocket_centered[radius] = pocket_centered.copy()
                        ref_ligand_rel_pocket[radius] = ligand_rel_pocket.copy()
                times.append(time_ns)
                t_whole_vals.append(0.0)
                theta_whole_vals.append(0.0)
                i_vals.append(0.0)
                for radius in POCKET_RADII_A:
                    t_pocket_vals[radius].append(0.0)
                    theta_pocket_vals[radius].append(0.0)
                continue

            t_whole, theta_whole, _ = frame_relative_theta(
                ca_centered, ref_ca_centered, ligand_rel_whole, ref_ligand_rel_whole
            )
            t_whole_vals.append(t_whole)
            theta_whole_vals.append(theta_whole)

            # I(t): pure ligand self-superposition, no frame reference (unchanged from script 33).
            q_selfshape, _ = rotation_matrix(ligand_self_centered, ref_ligand_self_centered)
            ligand_self_aligned = ligand_self_centered @ q_selfshape.T
            i_val = float(
                np.sqrt(np.mean(np.sum((ligand_self_aligned - ref_ligand_self_centered) ** 2, axis=1)))
            )
            i_vals.append(i_val)

            for radius in POCKET_RADII_A:
                if pocket_frame_data[radius] is None or ref_pocket_centered[radius] is None:
                    t_pocket_vals[radius].append(np.nan)
                    theta_pocket_vals[radius].append(np.nan)
                    continue
                pocket_centered, ligand_rel_pocket = pocket_frame_data[radius]
                t_pocket, theta_pocket, _ = frame_relative_theta(
                    pocket_centered, ref_pocket_centered[radius], ligand_rel_pocket, ref_ligand_rel_pocket[radius]
                )
                t_pocket_vals[radius].append(t_pocket)
                theta_pocket_vals[radius].append(theta_pocket)

            times.append(time_ns)

        times_arr = np.array(times)
        last_time_ns = float(times_arr[-1])
        late_mask = (times_arr > LATE_WINDOW_NS[0]) & (times_arr <= LATE_WINDOW_NS[1] + 1e-4)
        complete = bool(
            last_time_ns >= LATE_WINDOW_NS[1] - COMPLETE_TOLERANCE_NS and late_mask.sum() > 0
        )

        result.update(
            status="ok",
            n_frames=len(universe.trajectory),
            frame_interval_ns=dt_ns,
            trajectory_end_ns=last_time_ns,
            complete_70_100_coverage=complete,
        )
        for radius in POCKET_RADII_A:
            result[f"n_pocket_residues_{radius:g}A"] = n_pocket_residues[radius]

        def window_average(series: list[float], checkpoint: float) -> float:
            arr = np.array(series)
            mask = (times_arr > checkpoint - WINDOW_NS) & (times_arr <= checkpoint + 1e-4)
            return float(np.nanmean(arr[mask])) if mask.any() else np.nan

        for checkpoint in CHECKPOINTS_NS:
            suffix = f"{checkpoint - WINDOW_NS:g}_{checkpoint:g}"
            result[f"T_whole_mean_{suffix}"] = window_average(t_whole_vals, checkpoint)
            result[f"Theta_whole_mean_{suffix}"] = window_average(theta_whole_vals, checkpoint)
            result[f"I_mean_{suffix}"] = window_average(i_vals, checkpoint)
            for radius in POCKET_RADII_A:
                tag = f"{radius:g}A"
                result[f"T_pocket_{tag}_mean_{suffix}"] = window_average(t_pocket_vals[radius], checkpoint)
                result[f"Theta_pocket_{tag}_mean_{suffix}"] = window_average(
                    theta_pocket_vals[radius], checkpoint
                )
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
            "Optional CSV with a complex_id column to restrict the run to. "
            "If it has an eligible_primary_cohort column, only True rows are kept."
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

    print(f"Measuring local-pocket-alignment control for {len(runs)} trajectories", flush=True)

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
