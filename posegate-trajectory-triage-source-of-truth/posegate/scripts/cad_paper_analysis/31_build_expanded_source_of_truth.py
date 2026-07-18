#!/usr/bin/env python3
"""Build the canonical expanded trajectory source-of-truth dataset.

The input is the direct raw-trajectory measurement table produced by
28_generate_full_100ns_measurements.py. The output retains every discovered
trajectory and makes all cohort decisions explicit rather than dropping rows.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import sys

import MDAnalysis
import numpy as np
import pandas as pd


CHECKPOINTS = (5, 10, 15, 20, 30)
PRIMARY_CHECKPOINT_NS = 20
LATE_WINDOW = "(70,100] ns"
POSE_RETENTION_THRESHOLD_A = 3.0
EXPECTED_FRAME_COUNT = 500
EXPECTED_FRAME_INTERVAL_NS = 0.2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exclusion_reason(row: pd.Series) -> str:
    reasons = []
    if row["measurement_status"] != "ok":
        reasons.append("measurement_error")
    if not bool(row["complete_100ns_coverage"]):
        reasons.append("incomplete_0_100ns_dcd")
    if not bool(row["has_all_required_early_features"]):
        reasons.append("missing_pre20ns_feature")
    if not bool(row["has_late_outcome_measurement"]):
        reasons.append("missing_70_100ns_outcome")
    return ";".join(reasons) if reasons else "eligible_primary_cohort"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--measurements", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--generator", type=Path)
    args = parser.parse_args()

    raw = pd.read_csv(args.measurements)
    if raw["complex_id"].duplicated().any():
        duplicated = raw.loc[raw["complex_id"].duplicated(False), "complex_id"].tolist()
        raise ValueError(f"Duplicate complex_id values: {duplicated}")

    early_columns = []
    rename = {
        "status": "measurement_status",
        "error": "measurement_error",
        "complete_70_100_coverage": "complete_100ns_coverage",
        "late_pose_rmsd_median_70_100_A": "late_pose_rmsd_median_70_100_A",
    }
    for checkpoint in CHECKPOINTS:
        old_rmsd = f"pose_rmsd_mean_{checkpoint - 2}_{checkpoint}_A"
        old_disp = f"centroid_displacement_mean_{checkpoint - 2}_{checkpoint}_A"
        new_rmsd = f"corrected_pose_rmsd_mean_{checkpoint - 2}_{checkpoint}_A"
        new_disp = f"corrected_centroid_displacement_mean_{checkpoint - 2}_{checkpoint}_A"
        rename[old_rmsd] = new_rmsd
        rename[old_disp] = new_disp
        early_columns.extend((new_rmsd, new_disp))

    table = raw.rename(columns=rename).copy()
    table.insert(0, "source_of_truth_version", "expanded_same_trajectory_v1")
    table.insert(1, "audit_row_id", np.arange(1, len(table) + 1))
    table["expected_frame_count"] = EXPECTED_FRAME_COUNT
    table["expected_frame_interval_ns"] = EXPECTED_FRAME_INTERVAL_NS
    table["primary_checkpoint_ns"] = PRIMARY_CHECKPOINT_NS
    table["primary_feature_window"] = "(18,20] ns"
    table["late_outcome_window"] = LATE_WINDOW
    table["pose_retention_threshold_A"] = POSE_RETENTION_THRESHOLD_A
    table["has_all_required_early_features"] = table[early_columns].notna().all(axis=1)
    table["has_late_outcome_measurement"] = table[
        "late_pose_rmsd_median_70_100_A"
    ].notna()
    table["eligible_primary_cohort"] = (
        table["measurement_status"].eq("ok")
        & table["complete_100ns_coverage"].fillna(False)
        & table["has_all_required_early_features"]
        & table["has_late_outcome_measurement"]
    )
    table["primary_exclusion_reason"] = table.apply(exclusion_reason, axis=1)
    table["late_pose_retained"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    eligible = table["eligible_primary_cohort"]
    table.loc[eligible, "late_pose_retained"] = (
        table.loc[eligible, "late_pose_rmsd_median_70_100_A"]
        < POSE_RETENTION_THRESHOLD_A
    ).astype(int)

    r5 = "corrected_pose_rmsd_mean_3_5_A"
    r10 = "corrected_pose_rmsd_mean_8_10_A"
    r15 = "corrected_pose_rmsd_mean_13_15_A"
    r20 = "corrected_pose_rmsd_mean_18_20_A"
    d5 = "corrected_centroid_displacement_mean_3_5_A"
    d10 = "corrected_centroid_displacement_mean_8_10_A"
    d15 = "corrected_centroid_displacement_mean_13_15_A"
    d20 = "corrected_centroid_displacement_mean_18_20_A"
    table["rmsd_change_5_20_A"] = table[r20] - table[r5]
    table["rmsd_change_15_20_A"] = table[r20] - table[r15]
    table["rmsd_acceleration_A"] = (table[r20] - table[r15]) - (table[r15] - table[r10])
    table["centroid_change_5_20_A"] = table[d20] - table[d5]
    table["centroid_change_15_20_A"] = table[d20] - table[d15]
    table["centroid_acceleration_A"] = (table[d20] - table[d15]) - (table[d15] - table[d10])

    table = table.sort_values(["protein", "pocket", "replica"]).reset_index(drop=True)
    table["audit_row_id"] = np.arange(1, len(table) + 1)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    master_path = args.output_dir / "expanded_trajectory_source_of_truth.csv"
    table.to_csv(master_path, index=False)

    flow = pd.DataFrame(
        [
            {"stage": "raw_dcd_discovered", "trajectory_n": len(table), "protein_n": table["protein"].nunique()},
            {"stage": "measurement_success", "trajectory_n": int(table["measurement_status"].eq("ok").sum()), "protein_n": table.loc[table["measurement_status"].eq("ok"), "protein"].nunique()},
            {"stage": "complete_0_100ns_dcd", "trajectory_n": int(table["complete_100ns_coverage"].sum()), "protein_n": table.loc[table["complete_100ns_coverage"], "protein"].nunique()},
            {"stage": "all_required_early_features", "trajectory_n": int(table["has_all_required_early_features"].sum()), "protein_n": table.loc[table["has_all_required_early_features"], "protein"].nunique()},
            {"stage": "primary_complete_cohort", "trajectory_n": int(table["eligible_primary_cohort"].sum()), "protein_n": table.loc[table["eligible_primary_cohort"], "protein"].nunique()},
        ]
    )
    flow.to_csv(args.output_dir / "cohort_flow.csv", index=False)
    table.loc[~table["eligible_primary_cohort"], [
        "complex_id", "protein", "pocket", "n_frames", "trajectory_end_ns",
        "primary_exclusion_reason", "trajectory_path"
    ]].to_csv(args.output_dir / "primary_cohort_exclusions.csv", index=False)

    retained = table.loc[eligible, "late_pose_retained"]
    code_dir = Path(__file__).resolve().parent
    code_files = [
        code_dir / "28_generate_full_100ns_measurements.py",
        code_dir / "29_analyze_expanded_100ns_cohort.py",
        code_dir / "30_physics_guided_ai_audit.py",
        code_dir / "31_build_expanded_source_of_truth.py",
        code_dir / "run_expanded_source_of_truth_pipeline.py",
    ]
    summary = {
        "source_of_truth_version": "expanded_same_trajectory_v1",
        "master_rows": int(len(table)),
        "master_proteins": int(table["protein"].nunique()),
        "primary_cohort_rows": int(eligible.sum()),
        "primary_cohort_proteins": int(table.loc[eligible, "protein"].nunique()),
        "late_retained": int((retained == 1).sum()),
        "late_non_retained": int((retained == 0).sum()),
        "measurements_sha256": sha256(args.measurements),
        "master_sha256": sha256(master_path),
        "generator_sha256": sha256(args.generator) if args.generator else None,
        "code_sha256": {
            path.name: sha256(path) for path in code_files if path.exists()
        },
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "MDAnalysis": MDAnalysis.__version__,
    }
    (args.output_dir / "build_manifest.json").write_text(json.dumps(summary, indent=2) + "\n")

    readme = f"""# Expanded same-trajectory source of truth

Canonical version: `expanded_same_trajectory_v1`.

- Discovered raw trajectories: {len(table)} across {table['protein'].nunique()} proteins.
- Primary complete cohort: {int(eligible.sum())} across {table.loc[eligible, 'protein'].nunique()} proteins.
- Late retained: {int((retained == 1).sum())}.
- Late non-retained: {int((retained == 0).sum())}.

The master CSV retains every discovered trajectory. Use
`eligible_primary_cohort == True` for the primary 18--20 ns to 70--100 ns
same-trajectory analysis. Never infer eligibility from non-null labels alone;
use the explicit flag and audit `primary_exclusion_reason`.

Outcome: complete nominal 0--100 ns DCD coverage and median corrected ligand
pose RMSD below {POSE_RETENTION_THRESHOLD_A:g} A over {LATE_WINDOW}.

This dataset is retrospective. It is the row-level numerical source of truth,
not evidence of prospective validation or experimental binding.
"""
    (args.output_dir / "README.md").write_text(readme)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
