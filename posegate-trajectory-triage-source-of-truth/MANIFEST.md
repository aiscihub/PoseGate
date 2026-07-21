# Provenance manifest

Copied from `valleyfevermutation` at local HEAD commit `4c66c9c4a8023a9be8d4d4138e35f271c90c4eac`
(2026-06-26), copy performed 2026-07-18.

**Note:** in the source repo, `doc/YAU_competition_2026/yau_aj_trajectory_triage.tex`
and `pipeline/training/outputs/yau_recapture_screen/cad_states/` are
untracked working-tree content (`git status` reports `??`), not committed
to any `valleyfevermutation` commit. The HEAD hash above identifies the
repo state at copy time, not a commit that contains these files. If the
source repo later commits and revises these paths, this snapshot will
diverge and needs re-syncing by rerunning the copy.

## Verification performed at copy time

- `pipeline/training/outputs/yau_recapture_screen/cad_states/` — `diff -rq`
  against the source directory: **no differences** (full recursive
  byte-for-byte match, 15 MB, all subdirectories including
  `same_trajectory_source_of_truth/` and `expanded_source_of_truth/`).
- All 10 copied scripts (`common.py`, `25`–`31`, `run_expanded_source_of_truth_pipeline.py`,
  `threshold_and_robustness_sensitivity.py`) — `filecmp.cmp(..., shallow=False)`
  against the source files: **all OK** (byte-for-byte match).
- `doc/YAU_competition_2026/yau_aj_trajectory_triage.tex` — same check: **OK**.

## Copied files

### Scripts — `posegate/scripts/cad_paper_analysis/`
- `common.py` — shared path constants; defines `REPO_ROOT`/`OUT_DIR` by
  relative position, which is why this snapshot preserves the source repo's
  directory layout.
- `25_build_same_trajectory_source_of_truth.py` — builds the locked
  52-trajectory master table from `cad_states/` inputs.
- `26_plot_locked_claim_results.py` — regenerates the paper's headline
  figures and `key_claim_metrics.csv`.
- `27_audit_headlines_validation_checkpoints.py` — regenerates the headline
  registry, legacy reconciliation, and prospective-lock audit.
- `28_generate_full_100ns_measurements.py` — reads raw DCD trajectories via
  MDAnalysis. **Cannot run from this snapshot** (raw trajectories not
  included); output is already checked in under `cad_states/`.
- `29_analyze_expanded_100ns_cohort.py`, `30_physics_guided_ai_audit.py`,
  `31_build_expanded_source_of_truth.py` — build and audit the 73-trajectory
  expanded/exploratory cohort.
- `run_expanded_source_of_truth_pipeline.py` — single entrypoint chaining
  28→31→29→30; supports `--reuse-measurements` to skip the raw-trajectory
  step.
- `threshold_and_robustness_sensitivity.py` — threshold-sensitivity and
  non-obvious-subset robustness tables (paper Tables threshold-sensitivity,
  non-obvious).

### Data — `pipeline/training/outputs/yau_recapture_screen/cad_states/` (entire directory, 15 MB)

Includes every upstream input CSV plus both generated source-of-truth
result sets:
- `same_trajectory_source_of_truth/` — locked 52-trajectory cohort (paper's
  primary claims): master table, claim figures, headline registry, audits.
- `expanded_source_of_truth/` — 73-trajectory exploratory cohort (paper's
  model-complexity-audit claim only): master table, build manifest,
  primary/physics-guided-AI analysis outputs.
- Remaining top-level and subdirectory CSVs are the shared upstream inputs
  both cohorts are built from (`rbe_traces_pbc_corrected.csv`,
  `pocket_aware_state_20ns/`, `pocket_aware_validation2_100ns/`,
  `pocket_aware_validation3_model_comparison/`,
  `checkpoint_curve_and_baselines/`, `same_traj_vs_independent_replica/`,
  `same_trajectory_monitoring_validation/`, `replica_independence_audit/`,
  etc.).

### Report
- `doc/YAU_competition_2026/yau_aj_trajectory_triage.tex` — the manuscript
  itself, included for reference only (not needed to reproduce the data).

## 2026-07-20 addition: scripts 32-40 (code only)

Copied from `valleyfevermutation` branch `yau-trajectory-triage-code`
(commit `0cb16991`, which also holds `25`-`31` plus everything below —
that branch is a dedicated code-snapshot branch in the source repo, kept
separate from the paper/data working tree specifically so this kind of
code-only sync has a clean, citable commit to copy from). Verified
byte-for-byte identical to source via `filecmp.cmp(..., shallow=False)`
at copy time: all 10 files **OK**, no mismatches.

**No corresponding data was copied for this batch** (see README.md,
"Scripts 32-40" section, for exactly what each script needs and does not
have here). This is a deliberate scope choice, not an oversight: these
scripts' input CSVs (local-pocket-control measurements, PRPD
decomposition outputs, the cross-ligand 20 ns measurements) were not
part of the original 15 MB `cad_states/` snapshot and were not requested
to be added.

Files added: `32_contact_weighted_pocket_geometry_deviation.py`,
`33_pocket_relative_pose_decomposition.py`,
`34_prpd_frozen_lopo_joint_model.py`,
`35_local_pocket_alignment_control.py`,
`36_local_pocket_control_evaluation.py`,
`37_pocket_frame_diagnostic_audit.py`,
`38_pocket_frame_agreement_analysis.py`,
`39_plot_prpd_pocket_control_figure.py`,
`40_cross_ligand_20ns_coordinate_check.py`,
`test_synthetic_frame_relative_theta.py`.
