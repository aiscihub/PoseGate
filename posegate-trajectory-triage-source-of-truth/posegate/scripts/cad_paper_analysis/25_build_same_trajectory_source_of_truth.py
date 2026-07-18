#!/usr/bin/env python3
"""Lock the same-trajectory outcome and build one auditable master table.

This script implements manuscript TODO items A1 and A2 without rescanning raw
trajectories.  It joins the current cached analysis artifacts, verifies their
join keys, audits the late-outcome redundancy, and regenerates the primary
two-feature leave-one-protein-out (LOPO) predictions from the locked outcome.

Locked late outcome (A1)
------------------------
The legacy label was

    pose_rmsd_late < 3 A AND com_displacement_late < 5 A.

For the current definitions, at each frame the ligand-centroid displacement is
no greater than ligand RMSD.  The same ordering holds for their paired median
late-window summaries.  Therefore the 5 A centroid condition is redundant when
RMSD is below 3 A.

The legacy code also used the final 30% of *whatever frames were available*.
That admitted truncated trajectories and, in one case, made the nominal late
window overlap the 18--20 ns predictor.  The locked audit outcome therefore
also requires complete nominal 70--100 ns coverage:

    late_pose_retained_locked = 1[complete 70--100 ns coverage]
                                * 1[pose_rmsd_late < 3 A].

The script keeps both the stored and recomputed legacy labels and asserts that
the locked label is identical on every usable row.  It does not overwrite old
artifacts.

Source-of-truth universe (A2)
-----------------------------
The audit universe is the complete 64-row discovered separate-launch cohort.
Explicit membership flags recover the 61-row usable-feature cohort, legacy
59-row matched cohort, and strict complete-coverage outcome cohort. The
resulting wide table records identity, provenance, outcome components,
same-trajectory checkpoint features, separate-launch features, QC fields,
explicit eligibility flags, fold assignment, and regenerated out-of-fold
probabilities/decisions.

Outputs
-------
pipeline/training/outputs/yau_recapture_screen/cad_states/
same_trajectory_source_of_truth/
    same_trajectory_source_of_truth.csv
    outcome_redundancy_audit.csv
    cohort_flow.csv
    model_audit_metrics.csv
    source_artifact_manifest.csv
    README.md

Usage:
    python 25_build_same_trajectory_source_of_truth.py
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, recall_score, roc_auc_score
from sklearn.preprocessing import StandardScaler

from common import OUT_DIR, REPO_ROOT


DISCOVERY_TABLE = OUT_DIR / "pocket_aware_state_20ns" / "pocket_aware_state_20ns_full.csv"
VALIDATION3_FULL_TABLE = (
    OUT_DIR / "pocket_aware_validation3_model_comparison" / "validation3_full_features.csv"
)
MATCHED_COHORT_TABLE = (
    OUT_DIR / "pocket_aware_validation3_model_comparison" / "validation3_matched_cohort.csv"
)
CHECKPOINT_TABLE = (
    OUT_DIR / "checkpoint_curve_and_baselines" / "checkpoint_curve_and_baselines_full.csv"
)
SAME_VS_SEPARATE_TABLE = (
    OUT_DIR / "same_traj_vs_independent_replica" / "same_traj_vs_independent_replica_full.csv"
)
LEGACY_MODEL_TABLE = (
    OUT_DIR / "same_trajectory_monitoring_validation" / "sametraj_model_comparison_predictions.csv"
)
REPLICA_AUDIT_TABLE = (
    OUT_DIR / "replica_independence_audit" / "replica_independence_audit_full.csv"
)
OUTCOME_DETAIL_TABLE = (
    OUT_DIR / "pocket_aware_validation2_100ns" / "validation2_100ns_outcomes_full.csv"
)
PBC_TRACE_TABLE = OUT_DIR / "rbe_traces_pbc_corrected.csv"

RESULT_DIR = OUT_DIR / "same_trajectory_source_of_truth"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_TABLE = RESULT_DIR / "same_trajectory_source_of_truth.csv"

RMSD_THRESHOLD_A = 3.0
LEGACY_CENTROID_THRESHOLD_A = 5.0
NUMERICAL_TOLERANCE_A = 1e-9
EXPECTED_FRAME_INTERVAL_NS = 0.2
REQUIRED_TRAJECTORY_END_NS = 100.0
COVERAGE_TOLERANCE_NS = 0.05
CHECKPOINTS_NS = (5, 10, 15, 20, 30)
PRIMARY_FEATURES = (
    "same_traj_pose_rmsd_18_20_A",
    "same_traj_centroid_displacement_18_20_A",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_unique(path: Path, keys: list[str]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"required source artifact is missing: {path}")
    table = pd.read_csv(path)
    missing = [key for key in keys if key not in table.columns]
    if missing:
        raise ValueError(f"{path} is missing join keys: {missing}")
    duplicated = table.duplicated(keys, keep=False)
    if duplicated.any():
        examples = table.loc[duplicated, keys].head().to_dict("records")
        raise ValueError(f"{path} has duplicate join keys; examples: {examples}")
    return table


def nullable_int(values: pd.Series) -> pd.Series:
    return pd.to_numeric(values, errors="coerce").round().astype("Int64")


def recompute_lopo_pose_model(table: pd.DataFrame) -> pd.DataFrame:
    """Regenerate the current two-feature LOPO model from locked labels.

    This intentionally mirrors 22_same_trajectory_monitoring_validation.py:
    training-only StandardScaler, an explicitly parameterized L2 logistic
    regression matching the former sklearn defaults, and probability >= 0.5
    as the continue decision. No tuning or inner cross-validation is
    performed; that limitation is explicit in the output and README.
    """
    result = pd.DataFrame(index=table.index)
    result["locked_lopo_fold"] = pd.Series(pd.NA, index=table.index, dtype="string")
    result["locked_lopo_train_n"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    result["locked_lopo_train_proteins"] = pd.Series(pd.NA, index=table.index, dtype="Int64")
    result["locked_lopo_probability_retained"] = np.nan
    result["locked_lopo_decision_continue"] = pd.Series(pd.NA, index=table.index, dtype="Int64")

    eligible = table["eligible_locked_primary_pose_model"].fillna(False)
    model_table = table.loc[eligible].copy()
    for held_out_protein in sorted(model_table["protein"].unique()):
        test_index = model_table.index[model_table["protein"].eq(held_out_protein)]
        train = model_table.loc[~model_table.index.isin(test_index)]
        if train["late_pose_retained_locked"].nunique() < 2:
            raise ValueError(f"LOPO training fold for {held_out_protein} has one outcome class")

        scaler = StandardScaler(with_mean=True, with_std=True)
        x_train = scaler.fit_transform(train[list(PRIMARY_FEATURES)].to_numpy())
        model = LogisticRegression(
            penalty="l2",
            C=1.0,
            fit_intercept=True,
            class_weight="balanced",
            solver="lbfgs",
            max_iter=100,
            tol=1e-4,
            random_state=7,
        )
        model.fit(x_train, train["late_pose_retained_locked"].astype(int).to_numpy())

        x_test = scaler.transform(model_table.loc[test_index, list(PRIMARY_FEATURES)].to_numpy())
        probability = model.predict_proba(x_test)[:, 1]
        decision = (probability >= 0.5).astype(int)

        result.loc[test_index, "locked_lopo_fold"] = f"held_out_protein={held_out_protein}"
        result.loc[test_index, "locked_lopo_train_n"] = len(train)
        result.loc[test_index, "locked_lopo_train_proteins"] = train["protein"].nunique()
        result.loc[test_index, "locked_lopo_probability_retained"] = probability
        result.loc[test_index, "locked_lopo_decision_continue"] = decision

    return result


def exclusion_reason(row: pd.Series) -> str:
    missing: list[str] = []
    if not bool(row["late_metrics_available"]):
        missing.append("late_metrics")
    if not bool(row["late_window_70_100_coverage_ok"]):
        missing.append("incomplete_70_100ns_coverage")
    if pd.isna(row[PRIMARY_FEATURES[0]]):
        missing.append("same_traj_pose_rmsd_18_20")
    if pd.isna(row[PRIMARY_FEATURES[1]]):
        missing.append("same_traj_centroid_displacement_18_20")
    return "eligible" if not missing else "missing:" + ";".join(missing)


def make_cohort_flow(table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []

    def add(stage: str, mask: pd.Series, definition: str) -> None:
        subset = table.loc[mask.fillna(False)]
        labels = subset["late_pose_retained_locked"].dropna().astype(int)
        rows.append(
            {
                "stage": stage,
                "definition": definition,
                "n_complexes": len(subset),
                "n_proteins": subset["protein"].nunique(),
                "n_locked_labels_available": len(labels),
                "n_late_retained": int(labels.sum()) if len(labels) else 0,
                "n_late_nonretained": int((labels == 0).sum()) if len(labels) else 0,
            }
        )

    all_rows = pd.Series(True, index=table.index)
    add("audit_universe", all_rows, "complete discovered separate-launch cohort")
    add(
        "separate_launch_usable_20ns",
        table["separate_launch_usable_20ns"].fillna(False),
        "separate-launch 18-20 ns pose and retention inputs usable",
    )
    add(
        "current_matched_cohort",
        table["in_current_59row_matched_cohort"],
        "membership in the legacy 59-row feature/outcome cohort",
    )
    add("legacy_late_metrics", table["late_metrics_available"], "cached final-30%-of-available-frame metrics")
    add(
        "locked_late_outcome",
        table["late_outcome_usable"],
        "late metrics plus confirmed complete nominal 70-100 ns coverage",
    )
    for checkpoint in CHECKPOINTS_NS:
        add(
            f"checkpoint_{checkpoint}ns_rmsd",
            table[f"eligible_checkpoint_{checkpoint}ns_rmsd"],
            f"locked late outcome plus corrected RMSD at {checkpoint} ns",
        )
        add(
            f"checkpoint_{checkpoint}ns_pose_pair",
            table[f"eligible_checkpoint_{checkpoint}ns_pose_pair"],
            f"locked late outcome plus corrected RMSD and centroid displacement at {checkpoint} ns",
        )
    add(
        "locked_checkpoint_common_cohort",
        table["eligible_locked_checkpoint_common_cohort"],
        "locked late outcome plus corrected RMSD at every 5/10/15/20/30 ns checkpoint",
    )
    add(
        "locked_primary_pose_model_20ns",
        table["eligible_locked_primary_pose_model"],
        "locked late outcome plus corrected RMSD and centroid displacement over (18,20] ns",
    )
    add(
        "legacy_saved_model_predictions",
        table["legacy_pose_prediction_available"],
        "rows carrying saved pred_pose_only from script 22",
    )
    return pd.DataFrame(rows)


def make_column_dictionary(table: pd.DataFrame) -> pd.DataFrame:
    """Create a compact, explicit dictionary for every master-table column."""
    exact = {
        "audit_row_id": ("identity", "Stable 1-based row number in the generated audit table."),
        "complex_id": ("identity", "Canonical protein|ligand|pocket|replica join key."),
        "protein": ("identity", "Protein identifier and LOPO grouping unit."),
        "ligand": ("identity", "Ligand identifier."),
        "pocket": ("identity", "Predicted pocket identifier."),
        "replica": ("identity", "Replica directory identifier."),
        "late_pose_rmsd_median_last30pct_A": (
            "legacy outcome",
            "Cached median corrected ligand heavy-atom RMSD over the final 30% of available frames; valid as nominal 70-100 ns only when coverage is complete.",
        ),
        "late_centroid_displacement_median_last30pct_A": (
            "legacy outcome diagnostic",
            "Cached median corrected unweighted ligand heavy-atom centroid displacement over the final 30% of available frames.",
        ),
        "late_pose_retained_locked": (
            "locked outcome",
            "1 if complete nominal 70-100 ns coverage and late corrected RMSD median is below 3 A; missing otherwise.",
        ),
        "locked_late_outcome_definition_id": (
            "locked outcome",
            "Versioned outcome identifier: late_pose_retained_rmsd_v1.",
        ),
        "late_window_70_100_coverage_ok": (
            "coverage QC",
            "True only when audited target-trajectory end time reaches nominal 100 ns.",
        ),
        "locked_lopo_probability_retained": (
            "model output",
            "Outer-LOPO out-of-fold probability of the locked late-retained class.",
        ),
        "locked_lopo_decision_continue": (
            "model output",
            "1=continue/predicted retained when OOF probability is at least 0.5; 0=stop.",
        ),
    }

    rows = []
    for column in table.columns:
        if column in exact:
            category, description = exact[column]
        elif column.startswith("same_traj_pose_rmsd_"):
            category = "same-trajectory feature"
            description = "Mean PBC-corrected, protein-CA-aligned ligand RMSD in the final 2 ns before the named checkpoint."
        elif column.startswith("same_traj_centroid_displacement_"):
            category = "same-trajectory feature"
            description = "Mean corrected unweighted ligand-centroid displacement in the final 2 ns before the named checkpoint."
        elif column.startswith("same_traj_pocket_retention_"):
            category = "same-trajectory state"
            description = "Mean baseline-neighborhood retention in the final 2 ns before the named checkpoint."
        elif column.startswith("separate_launch_"):
            category = "separate-launch source/feature"
            description = "Field from the separately launched 20 ns trajectory; exact quantity is encoded in the column name."
        elif column.startswith("baseline_"):
            category = "coordinate baseline"
            description = "Alternative coordinate-processing measurement used for controlled comparison."
        elif column.startswith("eligible_"):
            category = "cohort flag"
            description = "True when the row has the locked outcome and features required by the named analysis."
        elif column.startswith("audit_"):
            category = "audit assertion"
            description = "Row-level consistency or equivalence check named by the column."
        elif column.startswith("legacy_") or "stored_legacy" in column:
            category = "legacy provenance"
            description = "Value copied or recomputed from the pre-lock analysis for comparison only."
        elif column.startswith("source17_") or column.startswith("source21_"):
            category = "source provenance"
            description = "Cached field from the numbered upstream analysis script named in the prefix."
        elif column.startswith("trace_") or column.startswith("trajectory_"):
            category = "trajectory coverage"
            description = "Frame-count or time-coverage metadata used to audit the nominal 70-100 ns outcome window."
        elif column.startswith("late_") or column.startswith("locked_late_"):
            category = "outcome/QC"
            description = "Late-outcome component, definition, or coverage field named by the column."
        elif column.startswith("locked_lopo_"):
            category = "model output"
            description = "Explicit outer-LOPO fold, training-set, probability, or decision field."
        elif column.startswith("in_"):
            category = "cohort membership"
            description = "Membership in the named upstream cohort or artifact."
        else:
            category = "supporting provenance"
            description = "Supporting cached field retained for audit; see the generator mapping for its exact source."
        rows.append(
            {
                "column": column,
                "pandas_dtype": str(table[column].dtype),
                "category": category,
                "description": description,
            }
        )
    return pd.DataFrame(rows)


def build_table() -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    base = read_unique(DISCOVERY_TABLE, ["complex_id"])
    validation3 = read_unique(VALIDATION3_FULL_TABLE, ["complex_id"])
    matched = read_unique(MATCHED_COHORT_TABLE, ["complex_id"])
    checkpoint = read_unique(CHECKPOINT_TABLE, ["complex_id"])
    same = read_unique(SAME_VS_SEPARATE_TABLE, ["complex_id"])
    legacy = read_unique(LEGACY_MODEL_TABLE, ["complex_id"])
    replica = read_unique(REPLICA_AUDIT_TABLE, ["protein", "ligand", "pocket", "replica"])
    outcome_detail = read_unique(OUTCOME_DETAIL_TABLE, ["complex_id"])
    if not PBC_TRACE_TABLE.exists():
        raise FileNotFoundError(f"required source artifact is missing: {PBC_TRACE_TABLE}")
    traces = pd.read_csv(PBC_TRACE_TABLE, usecols=["complex_id", "time_ns"])
    trace_coverage = (
        traces.groupby("complex_id", as_index=False)
        .agg(
            trace_start_time_ns=("time_ns", "min"),
            trace_end_time_ns=("time_ns", "max"),
            trace_n_frames=("time_ns", "size"),
        )
    )

    validation3_by_id = validation3.set_index("complex_id").reindex(base["complex_id"])
    checkpoint_by_id = checkpoint.set_index("complex_id").reindex(base["complex_id"])

    # Start with the complete 64-row discovered separate-launch cohort. Rename
    # source columns immediately so separate-launch and same-trajectory
    # quantities cannot be confused in downstream auditing.
    table = pd.DataFrame(
        {
            "audit_row_id": np.arange(1, len(base) + 1),
            "complex_id": base["complex_id"],
            "protein": base["protein"],
            "ligand": base["ligand"],
            "pocket": base["pocket"],
            "replica": base["replica"],
            "audit_universe": "complete_discovered_separate_launch_cohort",
            "in_validation3_usable_feature_cohort": base["complex_id"].isin(
                validation3["complex_id"]
            ),
            "in_current_59row_matched_cohort": base["complex_id"].isin(matched["complex_id"]),
            "in_checkpoint_artifact": base["complex_id"].isin(checkpoint["complex_id"]),
            "separate_launch_topology_path": base["topology"],
            "separate_launch_trajectory_path": base["trajectory"],
            "separate_launch_qc_status": base["qc_status"],
            "separate_launch_usable_20ns": base["usable"],
            "separate_launch_notes": base["notes"],
            "separate_launch_n_terminal_pose_frames": base["n_terminal_frames_pose"],
            "separate_launch_pose_rmsd_18_20_A": base["pose_rmsd_18_20"],
            "separate_launch_centroid_displacement_18_20_A": base["com_disp_18_20"],
            "separate_launch_pocket_retention_18_20": base["r_pocket"],
            "separate_launch_n_baseline_neighbors": base["n_anchors"],
            "separate_launch_structural_state_20ns": base["s20_state"],
            "late_pose_rmsd_median_last30pct_A": validation3_by_id["pose_rmsd_late"].to_numpy(),
            "late_centroid_displacement_median_last30pct_A": validation3_by_id[
                "com_displacement_late"
            ].to_numpy(),
            "late_outcome_usable_stored": validation3_by_id["usable_100ns"].to_numpy(),
            "late_pose_retained_stored_legacy": nullable_int(
                pd.Series(validation3_by_id["y_stable100"].to_numpy(), index=base.index)
            ),
            "existing_contact_label_unstable": nullable_int(
                pd.Series(
                    validation3_by_id["existing_label_unstable"].to_numpy(), index=base.index
                )
            ),
            "existing_contact_label_usable": validation3_by_id["usable_existing"].to_numpy(),
            "existing_contact_fraction": validation3_by_id["existing_f_contact"].to_numpy(),
            "existing_rmsd_late": validation3_by_id["existing_rmsd_late"].to_numpy(),
            "existing_drift": validation3_by_id["existing_drift"].to_numpy(),
            "baseline_corrected_final_frame_pose_rmsd_20ns_A": checkpoint_by_id[
                "rmsd_corrected_final_frame_20ns"
            ].to_numpy(),
            "baseline_no_pbc_pose_rmsd_18_20_A": checkpoint_by_id[
                "rmsd_nopbc_20ns"
            ].to_numpy(),
            "baseline_no_pbc_centroid_displacement_18_20_A": checkpoint_by_id[
                "disp_nopbc_20ns"
            ].to_numpy(),
            "baseline_ligand_self_aligned_rmsd_18_20_A": checkpoint_by_id[
                "rmsd_ligself_20ns"
            ].to_numpy(),
        }
    )

    for checkpoint in CHECKPOINTS_NS:
        table[f"same_traj_pose_rmsd_{checkpoint}ns_A"] = checkpoint_by_id[
            f"rmsd_corrected_{checkpoint}ns"
        ].to_numpy()
        table[f"same_traj_centroid_displacement_{checkpoint}ns_A"] = checkpoint_by_id[
            f"disp_corrected_{checkpoint}ns"
        ].to_numpy()
        table[f"same_traj_pocket_retention_{checkpoint}ns"] = checkpoint_by_id[
            f"r_pocket_{checkpoint}ns"
        ].to_numpy()

    # Canonical aliases used by the primary model make the exact terminal
    # window visible in the column name.
    table[PRIMARY_FEATURES[0]] = table["same_traj_pose_rmsd_20ns_A"]
    table[PRIMARY_FEATURES[1]] = table["same_traj_centroid_displacement_20ns_A"]
    table["same_traj_pocket_retention_18_20"] = table["same_traj_pocket_retention_20ns"]

    same_select = same[
        [
            "complex_id",
            "r_pocket_sametraj",
            "pose_rmsd_18_20_sametraj",
            "com_disp_18_20_sametraj",
            "r_pocket_notes",
            "pose_notes",
        ]
    ].rename(
        columns={
            "r_pocket_sametraj": "source21_same_traj_pocket_retention_18_20",
            "pose_rmsd_18_20_sametraj": "source21_same_traj_pose_rmsd_18_20_A",
            "com_disp_18_20_sametraj": "source21_same_traj_centroid_displacement_18_20_A",
            "r_pocket_notes": "source21_pocket_notes",
            "pose_notes": "source21_pose_notes",
        }
    )
    table = table.merge(same_select, on="complex_id", how="left", validate="one_to_one")

    legacy_select = legacy[
        [
            "complex_id",
            "s20_state_sametraj",
            "pred_pose_only",
            "pred_retention_only",
            "pred_combined",
        ]
    ].rename(
        columns={
            "s20_state_sametraj": "legacy_same_traj_structural_state_20ns",
            "pred_pose_only": "legacy_saved_pose_decision_continue",
            "pred_retention_only": "legacy_saved_retention_decision_continue",
            "pred_combined": "legacy_saved_combined_decision_continue",
        }
    )
    table = table.merge(legacy_select, on="complex_id", how="left", validate="one_to_one")

    replica_select = replica.copy()
    replica_select["complex_id"] = (
        replica_select["protein"].astype(str)
        + "|"
        + replica_select["ligand"].astype(str)
        + "|"
        + replica_select["pocket"].astype(str)
        + "|"
        + replica_select["replica"].astype(str)
    )
    replica_columns = [
        "complex_id",
        "starting_structure_identical",
        "energy_check",
        "n_overlap_points",
        "raw_pe_corr",
        "fluctuation_corr",
        "mean_abs_pe_diff",
        "pe100_std",
    ]
    table = table.merge(
        replica_select[replica_columns], on="complex_id", how="left", validate="one_to_one"
    )

    coverage_detail = outcome_detail[
        ["complex_id", "n_frames_100ns", "n_late_frames", "notes_100ns"]
    ].rename(
        columns={
            "n_frames_100ns": "source17_n_frames_100ns",
            "n_late_frames": "source17_n_late_frames",
            "notes_100ns": "source17_outcome_notes",
        }
    )
    table = table.merge(coverage_detail, on="complex_id", how="left", validate="one_to_one")
    table = table.merge(trace_coverage, on="complex_id", how="left", validate="one_to_one")

    # Prefer the pose-outcome script's frame count where present; otherwise use
    # the independently cached PBC trace. The checkpoint pipeline uses a 0.2 ns
    # frame cadence, made explicit here for rows without a cached time trace.
    table["trajectory_frame_interval_ns"] = EXPECTED_FRAME_INTERVAL_NS
    table["trajectory_n_frames_audited"] = table["source17_n_frames_100ns"].fillna(
        table["trace_n_frames"]
    )
    table["trajectory_end_time_ns_audited"] = table["trace_end_time_ns"].fillna(
        table["trajectory_n_frames_audited"] * EXPECTED_FRAME_INTERVAL_NS
    )
    table["late_window_required_start_ns"] = 70.0
    table["late_window_required_end_ns"] = REQUIRED_TRAJECTORY_END_NS
    table["late_window_70_100_coverage_ok"] = (
        table["trajectory_end_time_ns_audited"]
        >= REQUIRED_TRAJECTORY_END_NS - COVERAGE_TOLERANCE_NS
    )
    table["late_window_coverage_status"] = np.where(
        table["late_window_70_100_coverage_ok"],
        "complete_70_100ns",
        "excluded_incomplete_before_100ns",
    )

    # A1: lock the simpler outcome and keep a row-level audit of every logical
    # implication and stored-label comparison.
    late_rmsd = table["late_pose_rmsd_median_last30pct_A"]
    late_centroid = table["late_centroid_displacement_median_last30pct_A"]
    table["late_metrics_available"] = late_rmsd.notna() & late_centroid.notna()
    table["late_outcome_usable"] = (
        table["late_metrics_available"] & table["late_window_70_100_coverage_ok"]
    )
    table["late_rmsd_minus_centroid_A"] = late_rmsd - late_centroid
    table["audit_centroid_le_rmsd"] = (
        table["late_rmsd_minus_centroid_A"] >= -NUMERICAL_TOLERANCE_A
    ).where(table["late_metrics_available"], pd.NA).astype("boolean")

    legacy_recomputed = (
        (late_rmsd < RMSD_THRESHOLD_A) & (late_centroid < LEGACY_CENTROID_THRESHOLD_A)
    )
    locked = late_rmsd < RMSD_THRESHOLD_A
    table["late_pose_retained_rmsd_only_all_metrics"] = nullable_int(
        locked.where(table["late_metrics_available"])
    )
    table["late_pose_retained_legacy_recomputed"] = nullable_int(
        legacy_recomputed.where(table["late_metrics_available"])
    )
    table["late_pose_retained_locked"] = nullable_int(locked.where(table["late_outcome_usable"]))
    table["audit_locked_equals_legacy_formula"] = (
        table["late_pose_retained_locked"].eq(table["late_pose_retained_legacy_recomputed"])
    ).where(table["late_outcome_usable"], pd.NA).astype("boolean")
    table["audit_locked_equals_stored_legacy"] = (
        table["late_pose_retained_locked"].eq(table["late_pose_retained_stored_legacy"])
    ).where(table["late_outcome_usable"], pd.NA).astype("boolean")
    table["locked_late_outcome_definition"] = (
        "complete nominal 70-100 ns coverage AND median corrected pose RMSD < 3 A"
    )
    table["locked_late_outcome_definition_id"] = "late_pose_retained_rmsd_v1"
    table["legacy_centroid_condition_status"] = "audited redundant; retained for provenance only"

    # Audit that independent cached computations of the same 20 ns quantities
    # agree.  Missing values remain explicit rather than being imputed.
    comparisons = [
        (
            "audit_source21_pose_matches_checkpoint",
            table["source21_same_traj_pose_rmsd_18_20_A"],
            table[PRIMARY_FEATURES[0]],
        ),
        (
            "audit_source21_centroid_matches_checkpoint",
            table["source21_same_traj_centroid_displacement_18_20_A"],
            table[PRIMARY_FEATURES[1]],
        ),
        (
            "audit_source21_retention_matches_checkpoint",
            table["source21_same_traj_pocket_retention_18_20"],
            table["same_traj_pocket_retention_18_20"],
        ),
    ]
    for name, left, right in comparisons:
        comparable = left.notna() & right.notna()
        table[name] = pd.Series(pd.NA, index=table.index, dtype="boolean")
        table.loc[comparable, name] = np.isclose(
            left[comparable], right[comparable], atol=1e-10, rtol=1e-10
        )

    for checkpoint in CHECKPOINTS_NS:
        rmsd_available = table[f"same_traj_pose_rmsd_{checkpoint}ns_A"].notna()
        centroid_available = table[
            f"same_traj_centroid_displacement_{checkpoint}ns_A"
        ].notna()
        table[f"eligible_checkpoint_{checkpoint}ns_rmsd"] = (
            table["late_outcome_usable"] & rmsd_available
        )
        table[f"eligible_checkpoint_{checkpoint}ns_pose_pair"] = (
            table["late_outcome_usable"] & rmsd_available & centroid_available
        )

    table["eligible_locked_checkpoint_common_cohort"] = table[
        [f"eligible_checkpoint_{checkpoint}ns_rmsd" for checkpoint in CHECKPOINTS_NS]
    ].all(axis=1)

    table["eligible_locked_primary_pose_model"] = (
        table["late_outcome_usable"]
        & table[PRIMARY_FEATURES[0]].notna()
        & table[PRIMARY_FEATURES[1]].notna()
    )
    table["locked_primary_pose_model_exclusion_reason"] = table.apply(exclusion_reason, axis=1)
    table["legacy_pose_prediction_available"] = table[
        "legacy_saved_pose_decision_continue"
    ].notna()

    lopo = recompute_lopo_pose_model(table)
    for column in lopo.columns:
        table[column] = lopo[column]

    table["audit_locked_decision_matches_legacy_saved"] = pd.Series(
        pd.NA, index=table.index, dtype="boolean"
    )
    comparable_decisions = (
        table["locked_lopo_decision_continue"].notna()
        & table["legacy_saved_pose_decision_continue"].notna()
    )
    table.loc[comparable_decisions, "audit_locked_decision_matches_legacy_saved"] = (
        table.loc[comparable_decisions, "locked_lopo_decision_continue"].astype(int).to_numpy()
        == table.loc[comparable_decisions, "legacy_saved_pose_decision_continue"]
        .astype(int)
        .to_numpy()
    )

    # Hard audit gates: generation fails instead of writing a table if the A1
    # equivalence or independent cache consistency is violated.
    usable = table["late_outcome_usable"]
    metrics_available = table["late_metrics_available"]
    if table.loc[metrics_available, "trajectory_end_time_ns_audited"].isna().any():
        raise AssertionError("at least one cached late outcome lacks trajectory coverage metadata")
    if not table.loc[metrics_available, "audit_centroid_le_rmsd"].all():
        raise AssertionError("stored late summaries violate centroid displacement <= RMSD")
    if not table.loc[usable, "audit_locked_equals_legacy_formula"].all():
        raise AssertionError("locked RMSD-only outcome changes at least one legacy formula label")
    if not table.loc[usable, "audit_locked_equals_stored_legacy"].all():
        raise AssertionError("locked outcome disagrees with at least one stored y_stable100 label")
    for name, _, _ in comparisons:
        comparable = table[name].notna()
        if comparable.any() and not table.loc[comparable, name].all():
            raise AssertionError(f"independent cached values disagree: {name}")
    # Regenerated decisions use the strict coverage cohort and are therefore
    # not expected to match a model trained on the legacy truncated cohort.

    outcome_audit = pd.DataFrame(
        [
            {
                "audit_universe_n": len(table),
                "late_metrics_available_n": int(metrics_available.sum()),
                "usable_late_outcome_n": int(usable.sum()),
                "excluded_incomplete_coverage_n": int(
                    (metrics_available & ~table["late_window_70_100_coverage_ok"]).sum()
                ),
                "locked_rmsd_threshold_A": RMSD_THRESHOLD_A,
                "legacy_centroid_threshold_A": LEGACY_CENTROID_THRESHOLD_A,
                "late_summary_statistic": "median over final 30% of frames",
                "min_stored_rmsd_minus_centroid_A": float(
                    table.loc[metrics_available, "late_rmsd_minus_centroid_A"].min()
                ),
                "centroid_gt_rmsd_violations": int(
                    (~table.loc[metrics_available, "audit_centroid_le_rmsd"]).sum()
                ),
                "locked_vs_legacy_formula_mismatches": int(
                    (~table.loc[usable, "audit_locked_equals_legacy_formula"]).sum()
                ),
                "locked_vs_stored_label_mismatches": int(
                    (~table.loc[usable, "audit_locked_equals_stored_legacy"]).sum()
                ),
                "rmsd_only_vs_legacy_formula_mismatches_all_metrics": int(
                    (
                        table.loc[metrics_available, "late_pose_retained_rmsd_only_all_metrics"]
                        != table.loc[
                            metrics_available, "late_pose_retained_legacy_recomputed"
                        ]
                    ).sum()
                ),
                "locked_late_retained_n": int(
                    table.loc[usable, "late_pose_retained_locked"].sum()
                ),
                "locked_late_nonretained_n": int(
                    (table.loc[usable, "late_pose_retained_locked"] == 0).sum()
                ),
                "conclusion": (
                    "centroid < 5 A condition is redundant under RMSD < 3 A; "
                    "locked outcome is RMSD-only on complete nominal 70-100 ns trajectories"
                ),
            }
        ]
    )

    model_rows = table["eligible_locked_primary_pose_model"]
    model_table = table.loc[model_rows].copy()
    y_true = model_table["late_pose_retained_locked"].astype(int).to_numpy()
    probability = model_table["locked_lopo_probability_retained"].to_numpy()
    decision = model_table["locked_lopo_decision_continue"].astype(int).to_numpy()
    tn = int(((decision == 0) & (y_true == 0)).sum())
    fn = int(((decision == 0) & (y_true == 1)).sum())
    tp = int(((decision == 1) & (y_true == 1)).sum())
    fp = int(((decision == 1) & (y_true == 0)).sum())
    model_audit = pd.DataFrame(
        [
            {
                "model": "two_feature_corrected_pose_lopo",
                "features": ";".join(PRIMARY_FEATURES),
                "outcome": "late_pose_retained_locked",
                "n": len(model_table),
                "n_proteins": model_table["protein"].nunique(),
                "n_retained": int(y_true.sum()),
                "n_nonretained": int((y_true == 0).sum()),
                "scaler": "StandardScaler(with_mean=True, with_std=True), fit within outer training fold",
                "penalty": "l2",
                "regularization_C": 1.0,
                "solver": "lbfgs",
                "max_iter": 100,
                "tolerance": 1e-4,
                "fit_intercept": True,
                "class_weight": "balanced",
                "decision_threshold": 0.5,
                "hyperparameter_tuning": "none",
                "pandas_version": pd.__version__,
                "scikit_learn_version": sklearn.__version__,
                "auroc_oof_probability": roc_auc_score(y_true, probability),
                "balanced_accuracy": balanced_accuracy_score(y_true, decision),
                "sensitivity_retained": recall_score(y_true, decision, pos_label=1),
                "specificity_nonretained": recall_score(y_true, decision, pos_label=0),
                "correct_stop_tn": tn,
                "false_stop_fn": fn,
                "correct_continue_tp": tp,
                "nonretained_continue_fp": fp,
            }
        ]
    )

    cohort_flow = make_cohort_flow(table)
    column_dictionary = make_column_dictionary(table)
    coverage_exclusions = table.loc[
        table["late_metrics_available"] & ~table["late_window_70_100_coverage_ok"],
        [
            "complex_id",
            "protein",
            "pocket",
            "trajectory_n_frames_audited",
            "trajectory_end_time_ns_audited",
            "late_pose_rmsd_median_last30pct_A",
            "late_centroid_displacement_median_last30pct_A",
            "late_pose_retained_stored_legacy",
            "late_window_coverage_status",
        ],
    ].sort_values("trajectory_end_time_ns_audited")

    source_paths = [
        DISCOVERY_TABLE,
        VALIDATION3_FULL_TABLE,
        MATCHED_COHORT_TABLE,
        CHECKPOINT_TABLE,
        SAME_VS_SEPARATE_TABLE,
        LEGACY_MODEL_TABLE,
        REPLICA_AUDIT_TABLE,
        OUTCOME_DETAIL_TABLE,
        PBC_TRACE_TABLE,
    ]
    manifest_rows = []
    for path in source_paths:
        source_table = pd.read_csv(path)
        manifest_rows.append(
            {
                "path_relative_to_repo": str(path.relative_to(REPO_ROOT)),
                "sha256": sha256(path),
                "rows": len(source_table),
                "columns": len(source_table.columns),
            }
        )
    manifest = pd.DataFrame(manifest_rows)

    return table, {
        "outcome_redundancy_audit.csv": outcome_audit,
        "cohort_flow.csv": cohort_flow,
        "coverage_exclusions.csv": coverage_exclusions,
        "column_dictionary.csv": column_dictionary,
        "model_audit_metrics.csv": model_audit,
        "source_artifact_manifest.csv": manifest,
    }


def write_readme(table: pd.DataFrame, outputs: dict[str, pd.DataFrame]) -> None:
    audit = outputs["outcome_redundancy_audit.csv"].iloc[0]
    metrics = outputs["model_audit_metrics.csv"].iloc[0]
    text = f"""# Same-trajectory source-of-truth audit

Generated by `yau_competition/scripts/cad_paper_analysis/25_build_same_trajectory_source_of_truth.py`.

## Locked A1 outcome

The late-window implementation stores the **median over the final 30% of trajectory frames**, not a mean. The legacy rule was:

```text
pose_rmsd_late < {RMSD_THRESHOLD_A:g} A AND centroid_displacement_late < {LEGACY_CENTROID_THRESHOLD_A:g} A
```

The locked rule is:

```text
late_pose_retained_locked = 1[
    trajectory covers the nominal 70-100 ns window
    AND
    median corrected pose_rmsd_late < {RMSD_THRESHOLD_A:g} A
]
```

For atom-displacement vectors `delta_r_a` and their centroid `delta_C`,

```text
RMSD^2 = ||delta_C||^2 + mean_a ||delta_r_a - delta_C||^2
```

so centroid displacement cannot exceed RMSD. The ordering also holds for their paired medians. In the current artifacts:

- cached late metrics: {int(audit['late_metrics_available_n'])};
- excluded for incomplete 70-100 ns coverage: {int(audit['excluded_incomplete_coverage_n'])};
- locked usable late outcomes: {int(audit['usable_late_outcome_n'])};
- minimum stored `RMSD - centroid displacement`: {audit['min_stored_rmsd_minus_centroid_A']:.6f} A;
- inequality violations: {int(audit['centroid_gt_rmsd_violations'])};
- RMSD-only versus legacy-formula changes across all cached metrics: {int(audit['rmsd_only_vs_legacy_formula_mismatches_all_metrics'])};
- locked-versus-legacy formula label changes: {int(audit['locked_vs_legacy_formula_mismatches'])};
- locked-versus-stored label changes: {int(audit['locked_vs_stored_label_mismatches'])}.

The centroid value remains in the table as a descriptive/predictor feature, but not as an independent late-label condition. The seven truncated trajectories retain their legacy values for provenance, but their locked outcomes are missing and they are excluded from locked modeling.

## A2 table scope

`same_trajectory_source_of_truth.csv` contains {len(table)} rows from the complete discovered separate-launch cohort. Explicit flags identify the 61-row usable-feature cohort, legacy 59-row matched cohort, and strict 52-row complete-coverage outcome cohort.

The table records:

- complex identity and source paths;
- stored and locked late outcomes;
- same-trajectory corrected RMSD, centroid displacement, and pocket retention at 5/10/15/20/30 ns;
- separate-launch 20 ns measurements;
- coordinate-processing baselines;
- QC and replica-audit fields;
- eligibility and exclusion reasons;
- LOPO fold, training counts, out-of-fold probability, and stop/continue decision;
- row-level agreement checks against independently cached artifacts.

## Regenerated primary model audit

- eligible rows: {int(metrics['n'])};
- proteins: {int(metrics['n_proteins'])};
- late retained/non-retained: {int(metrics['n_retained'])}/{int(metrics['n_nonretained'])};
- out-of-fold probability AUROC: {metrics['auroc_oof_probability']:.6f};
- balanced accuracy at probability 0.5: {metrics['balanced_accuracy']:.6f};
- sensitivity/specificity: {metrics['sensitivity_retained']:.6f}/{metrics['specificity_nonretained']:.6f};
- confusion counts `(TN, FN, TP, FP)`: ({int(metrics['correct_stop_tn'])}, {int(metrics['false_stop_fn'])}, {int(metrics['correct_continue_tp'])}, {int(metrics['nonretained_continue_fp'])}).

This reruns the script-22 two-feature procedure on the stricter complete-coverage cohort. It does not reproduce the legacy decisions exactly because six truncated trajectories formerly included in the 58-row headline model are now excluded. It uses no hyperparameter tuning or inner cross-validation; that limitation is explicit rather than hidden.

## Supporting files

- `outcome_redundancy_audit.csv`: one-row A1 audit.
- `cohort_flow.csv`: exact inclusion counts at each checkpoint and model stage.
- `coverage_exclusions.csv`: the seven cached outcomes excluded for ending before 100 ns.
- `column_dictionary.csv`: name, dtype, category, and definition for every master-table column.
- `model_audit_metrics.csv`: regenerated primary LOPO metrics and confusion counts.
- `source_artifact_manifest.csv`: input paths, sizes, and SHA-256 hashes.

## Downstream A3-A5 audit

After this table is generated, run
`27_audit_headlines_validation_checkpoints.py`. It writes the exact headline
registry, legacy-number reconciliation, validation/model/fold audit, common-cohort
checkpoint results, and the dated 20 ns prospective lock. Those files are
derived from this source-of-truth table; legacy aggregates are used only to
explain superseded numbers.

Do not hand-edit generated files. Change the script or upstream artifacts and rerun it.
"""
    (RESULT_DIR / "README.md").write_text(text)


def main() -> None:
    table, supporting = build_table()
    table.to_csv(OUTPUT_TABLE, index=False)
    for filename, output in supporting.items():
        output.to_csv(RESULT_DIR / filename, index=False)
    write_readme(table, supporting)

    audit = supporting["outcome_redundancy_audit.csv"].iloc[0]
    metrics = supporting["model_audit_metrics.csv"].iloc[0]
    print(f"Wrote {len(table)} audit rows to {OUTPUT_TABLE}")
    print(
        "A1 outcome audit: "
        f"violations={int(audit['centroid_gt_rmsd_violations'])}, "
        f"label_changes={int(audit['locked_vs_legacy_formula_mismatches'])}"
    )
    print(
        "Locked LOPO model: "
        f"n={int(metrics['n'])}, BA={metrics['balanced_accuracy']:.6f}, "
        f"AUROC(probability)={metrics['auroc_oof_probability']:.6f}"
    )
    print(f"pandas={pd.__version__}; scikit-learn={sklearn.__version__}")


if __name__ == "__main__":
    main()
