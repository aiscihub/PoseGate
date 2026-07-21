#!/usr/bin/env python3
"""PROTOTYPE: Contact-Weighted Pocket Geometry Deviation (CW-PGD).

Not part of the locked or expanded source-of-truth pipelines. This is a
candidate early predictor, evaluated the same way as the existing
same-trajectory checkpoint features (script 24) but computed independently
here so it can be tested without touching locked artifacts.

Motivation: the validated corrected-RMSD measurement (script 24 / 28)
depends on a global protein-Cα Kabsch alignment -- it is sensitive to
domain motion unrelated to the ligand leaving its pocket, and it weights
every ligand heavy atom equally regardless of whether that atom actually
contacts the protein. CW-PGD instead asks how much the ligand's geometry
relative to its *own early contacts* has changed, weighting atom pairs by
how consistently they were in contact during a 0-1 ns baseline window.

Algorithm, per trajectory (two passes -- baseline contacts must be known
before the full-trajectory pass can be restricted to them):

  Pass 1 (frames with time_ns in (0, 1] ns only):
    1. Restrict candidate pocket atoms to protein heavy atoms within
       NEIGHBOR_SHELL_A of any ligand heavy atom at frame 0 (a generous
       shell so nothing that could plausibly contact is excluded).
    2. For every (ligand atom, pocket atom) pair, compute minimum-image
       distance d_ia(t) at each baseline frame via MDAnalysis
       `distance_array` (fully PBC-aware pairwise distance, not a
       centroid-relative approximation).
    3. Baseline occupancy p_ia = fraction of baseline frames with
       d_ia(t) <= CONTACT_CUTOFF_A.
    4. Baseline reference distance d_ia^0 = median_{t in baseline} d_ia(t)
       (median, not the frame-0 value alone, so one unusual starting
       snapshot cannot set the reference).
    5. Retain only pairs with p_ia > 0. Weight w_ia = p_ia^POWER,
       normalized so sum(w_ia) = 1 -- this makes G_pocket a proper
       weighted quadratic mean (its magnitude does not depend on how many
       contacts exist or how large the ligand/pocket are). POWER=1
       reproduces the original linear-occupancy-weight spec; POWER>1
       sharpens the preference for persistently-occupied contacts over
       intermittent ones (see --weight-power).

  Pass 2 (every frame, restricted to the retained pairs only -- cheap):
    G_pocket(t) = sqrt( sum_ia w_ia * (d_ia(t) - d_ia^0)^2 )

    averaged into the same trailing-2ns checkpoint windows and the same
    (70,100] ns late window used by the locked pipeline, so this column
    can be merged onto the existing cohort table and compared directly
    against `same_traj_pose_rmsd_18_20_A` with a paired bootstrap test.

Usage (against the 73-row expanded cohort, reusing raw trajectories
already discovered by script 28):
    python 32_contact_weighted_pocket_geometry_deviation.py \\
        --root /media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md \\
        --output ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/cwpgd_candidate_measurements.csv \\
        --complex-ids ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth/master/expanded_trajectory_source_of_truth.csv \\
        --workers 4

Evaluate the result by merging `cwpgd_mean_18_20_A` onto
`expanded_trajectory_source_of_truth.csv` (join key `complex_id`) and
running the same `auc_from_scores` / `cluster_bootstrap_auc` /
`paired_cluster_bootstrap_difference` machinery from
`26_plot_locked_claim_results.py` against `late_pose_retained_locked`,
with `sign=-1.0` (lower deviation -> more likely retained, same
convention as the RMSD features).
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import warnings

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
import numpy as np
import pandas as pd


LIGAND_RESNAME = "UNK"
CONTACT_CUTOFF_A = 4.5       # matches common.CONTACT_CUTOFF_A
NEIGHBOR_SHELL_A = 12.0      # generous frame-0 shell for candidate pocket atoms
BASELINE_WINDOW_NS = (0.0, 1.0)
MIN_BASELINE_FRAMES = 2
CHECKPOINTS_NS = (5, 10, 15, 20, 30)
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


def measure(run: tuple[str, str, str, Path, Path], weight_power: float) -> dict:
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

        ligand = universe.select_atoms(f"resname {LIGAND_RESNAME} and not name H*")
        protein_heavy = universe.select_atoms("protein and not name H*")
        if ligand.n_atoms == 0 or protein_heavy.n_atoms == 0:
            raise ValueError(
                f"missing selection: ligand={ligand.n_atoms}, protein_heavy={protein_heavy.n_atoms}"
            )

        dt_ns = float(getattr(universe.trajectory, "dt", 200.0)) / 1000.0
        if not np.isfinite(dt_ns) or dt_ns <= 0:
            dt_ns = 0.2

        # Candidate pocket shell from frame 0 only (cheap, generous cutoff).
        universe.trajectory[0]
        frame0_box = universe.trajectory.ts.dimensions
        d0_shell = distance_array(ligand.positions, protein_heavy.positions, box=frame0_box)
        shell_mask = np.any(d0_shell <= NEIGHBOR_SHELL_A, axis=0)
        pocket_atoms = protein_heavy[shell_mask]
        if pocket_atoms.n_atoms == 0:
            raise ValueError("no candidate pocket atoms within neighbor shell")

        # --- Pass 1: baseline occupancy + baseline reference distances ---
        baseline_distance_frames = []
        n_baseline_frames = 0
        for frame_i, ts in enumerate(universe.trajectory):
            time_ns = (frame_i + 1) * dt_ns
            if time_ns > BASELINE_WINDOW_NS[1] + 1e-4:
                break
            if BASELINE_WINDOW_NS[0] < time_ns <= BASELINE_WINDOW_NS[1] + 1e-4:
                box = ts.dimensions
                d = distance_array(ligand.positions, pocket_atoms.positions, box=box)
                baseline_distance_frames.append(d)
                n_baseline_frames += 1

        if n_baseline_frames < MIN_BASELINE_FRAMES:
            raise ValueError(f"insufficient baseline frames: {n_baseline_frames}")

        stacked = np.stack(baseline_distance_frames, axis=0)  # (n_baseline, n_lig, n_pocket)
        p_ia = np.mean(stacked <= CONTACT_CUTOFF_A, axis=0)
        d0_ia = np.median(stacked, axis=0)

        retained_mask = p_ia > 0.0
        n_retained = int(retained_mask.sum())
        if n_retained == 0:
            raise ValueError(f"no baseline contacts within {CONTACT_CUTOFF_A} A")

        li_idx, pk_idx = np.nonzero(retained_mask)
        raw_weight = p_ia[li_idx, pk_idx] ** weight_power
        w_flat = raw_weight / raw_weight.sum()
        d0_flat = d0_ia[li_idx, pk_idx]

        # --- Pass 2: full-trajectory G_pocket(t), restricted to retained pairs ---
        window_values = {checkpoint: [] for checkpoint in CHECKPOINTS_NS}
        late_values = []
        last_time_ns = np.nan

        universe.trajectory.rewind()
        for frame_i, ts in enumerate(universe.trajectory):
            time_ns = (frame_i + 1) * dt_ns
            last_time_ns = time_ns
            box = ts.dimensions
            d = distance_array(ligand.positions, pocket_atoms.positions, box=box)
            d_flat = d[li_idx, pk_idx]
            g = float(np.sqrt(np.sum(w_flat * (d_flat - d0_flat) ** 2)))

            for checkpoint in CHECKPOINTS_NS:
                if checkpoint - WINDOW_NS < time_ns <= checkpoint + 1e-4:
                    window_values[checkpoint].append(g)
            if LATE_WINDOW_NS[0] < time_ns <= LATE_WINDOW_NS[1] + 1e-4:
                late_values.append(g)

        result.update(
            status="ok",
            n_frames=len(universe.trajectory),
            frame_interval_ns=dt_ns,
            trajectory_end_ns=last_time_ns,
            n_pocket_atoms_candidate=int(pocket_atoms.n_atoms),
            n_baseline_frames=n_baseline_frames,
            n_retained_contacts=n_retained,
            weight_power=weight_power,
            complete_70_100_coverage=bool(
                np.isfinite(last_time_ns)
                and last_time_ns >= LATE_WINDOW_NS[1] - COMPLETE_TOLERANCE_NS
                and len(late_values) > 0
            ),
            late_cwpgd_median_70_100_A=(
                float(np.median(late_values)) if late_values else np.nan
            ),
        )
        for checkpoint in CHECKPOINTS_NS:
            result[f"cwpgd_mean_{checkpoint - 2}_{checkpoint}_A"] = (
                float(np.mean(window_values[checkpoint]))
                if window_values[checkpoint]
                else np.nan
            )
            result[f"n_frames_{checkpoint - 2}_{checkpoint}"] = len(window_values[checkpoint])
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
    parser.add_argument(
        "--weight-power",
        type=float,
        default=1.0,
        help="Exponent applied to baseline occupancy before normalizing to weights (default 1.0 = linear occupancy).",
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
        runs = [
            run
            for run in runs
            if f"{run[0]}|Milbemycin|{run[1]}|{run[2]}" in keep
        ]
        if not runs:
            raise SystemExit("No discovered trajectories matched --complex-ids")

    print(f"Measuring CW-PGD for {len(runs)} trajectories", flush=True)

    rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(measure, run, args.weight_power): run for run in runs
        }
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
