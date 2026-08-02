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

## 2026-08-02 addition: n=52 claim_figures refresh + n=81 expanded cohort (data only)

Re-synced `same_trajectory_source_of_truth/claim_figures/` from the source
repo (additive only — the 12 original files were unchanged; added 13 files
covering later analyses: `spearman_early_late_clustered.csv`,
`table2_5ns_paired_differences.csv`, `table2_5ns_coordinate_and_transfer.csv`,
`harmonized_cross_ligand_3_5ns_stats.csv`, `forecast_horizon_curve.csv`,
`protein_influence_jackknife.csv`, `strict_survivor_subset_auroc.csv`,
`same_vs_separate_5ns.csv`, `separate_launch_pose_rmsd_5ns.csv`,
`paired_auroc_difference_5ns_vs_20ns.csv`,
`consolidated_100ns_extension.csv`, `contact_retention_outcome.csv`,
`table2_5ns_recalculation_manifest.json`). All top-level n=52 files
(`same_trajectory_source_of_truth.csv`, `README.md`, `RESULTS_LOCK_A3_A5.md`,
etc.) were verified unchanged since the 2026-07-18 copy — this cohort is
locked/frozen, as expected.

Added `expanded_source_of_truth_frozen_n81/` in full (not previously
present here) — the current locked n=81 expanded stress-test cohort, which
supersedes the n=73 snapshot (`expanded_source_of_truth/`, still present
below for history) as the paper's expanded-cohort numbers as of 2026-07-30.
Includes `master/expanded_trajectory_source_of_truth.csv` (84 rows,
`eligible_primary_cohort` sums to 81), `checkpoint_and_robustness/` (5/10/
15/20/30ns AUROC+CI, matches Figure 2B and the Results-section prose), and
`physics_guided_ai_audit/` (the model-complexity audit behind Table 3,
`tab:model-audit-main`).

Data only — no additional scripts copied this round (scripts 25-31 already
here remain current for the n=52 build; the expanded-cohort/model-audit
scripts, and everything numbered above 31, were not requested and are not
included).
