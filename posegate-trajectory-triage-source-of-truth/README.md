# PoseGate trajectory-triage source-of-truth (reproducibility snapshot)

Self-contained snapshot of the data and scripts backing
`doc/YAU_competition_2026/yau_aj_trajectory_triage.tex` (in
`valleyfevermutation`), copied 2026-07-18, code updated 2026-07-20 with
scripts 32-40. All files below are byte-identical to the source repo at
the date they were added; see `MANIFEST.md` for the per-file provenance
and verification note.

## Layout

The directory structure mirrors the source repo exactly, on purpose: the
analysis scripts locate their inputs/outputs relative to their own file
path (`common.py: REPO_ROOT = Path(__file__).resolve().parents[3]`), so
keeping this layout means the scripts run unmodified.

```
posegate-trajectory-triage-source-of-truth/
├── posegate/scripts/cad_paper_analysis/   # the analysis pipeline (scripts 25-40 + test suite + helpers)
├── pipeline/training/outputs/yau_recapture_screen/cad_states/  # data for scripts 25-31 only (15 MB)
│   ├── same_trajectory_source_of_truth/    # LOCKED 52-trajectory cohort — primary paper claims
│   └── expanded_source_of_truth/           # 73-trajectory exploratory cohort — model-complexity audit only
└── doc/YAU_competition_2026/yau_aj_trajectory_triage.tex  # the report itself, for reference
```

**Scripts 32-40 and `test_synthetic_frame_relative_theta.py` are code-only additions (2026-07-20): their supporting data is not included in this snapshot.** See "Scripts 32-40" below.

## What is and isn't reproducible from this snapshot

The raw MD trajectories (100 ns/20 ns OpenMM DCD runs) that everything
ultimately derives from are **not** included — they live outside any git
repo, at `common.py`'s `BASE_100NS` / `BASE_20NS`
(`/media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md` and
`simulation_20ns_md` on the original machine) and are multi-GB. This
snapshot starts one level up, from the already-computed per-frame/per-complex
CSVs in `cad_states/` — that is the practical "source of truth" the paper's
numbers are drawn from.

**Fully reproducible from this snapshot alone:**

```bash
cd posegate/scripts/cad_paper_analysis

# Locked 52-trajectory cohort: rebuilds same_trajectory_source_of_truth.csv
# and all supporting audit CSVs from the cad_states/ inputs already present.
python3 25_build_same_trajectory_source_of_truth.py

# Regenerates the paper's headline figures + key_claim_metrics.csv from the
# table 25 just wrote.
python3 26_plot_locked_claim_results.py

# Regenerates the headline registry, legacy-number reconciliation, and the
# dated prospective-lock audit from the same locked table.
python3 27_audit_headlines_validation_checkpoints.py

# Expanded 73-trajectory cohort (exploratory model-complexity audit cited in
# the paper's Results section): reuses the already-computed
# full_100ns_corrected_measurements.csv instead of re-deriving it from raw
# DCD trajectories.
python3 run_expanded_source_of_truth_pipeline.py --reuse-measurements \
    --result-dir ../../../pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth
```

Each script reads only from `cad_states/` (via `common.OUT_DIR`) and writes
back into it, so rerunning should reproduce the checked-in outputs exactly
(diff the regenerated CSVs against the copies already in this snapshot to
confirm).

**Not reproducible from this snapshot:**

- `28_generate_full_100ns_measurements.py` — this is the one script that
  reads the raw DCD trajectories via MDAnalysis. It is included for
  completeness/documentation but will fail without `BASE_100NS` populated.
  Its output (`full_100ns_corrected_measurements.csv`) is already present
  under `cad_states/`, which is why the expanded-cohort command above uses
  `--reuse-measurements` to skip it.

## Scripts 32-40 (code-only, added 2026-07-20)

Later analysis work, added here as code for reference/provenance but
**not runnable in this snapshot** -- their input CSVs (e.g.
`local_pocket_control_measurements.csv`, `prpd_candidate_measurements.csv`,
the cross-ligand 20 ns measurement outputs) are not included:

- `32_contact_weighted_pocket_geometry_deviation.py` -- prototype contact-weighted pocket-geometry-deviation metric (superseded/exploratory, not used in the paper's headline results).
- `33_pocket_relative_pose_decomposition.py` -- decomposes corrected RMSD into translation (T), rotation (Theta), and internal deformation (I) components.
- `34_prpd_frozen_lopo_joint_model.py` -- frozen leave-one-protein-out logistic model on `[T, Theta]`, with the exact spec (training-fold-only standardization, `C=1.0`, balanced class weights, no tuning) frozen in advance.
- `35_local_pocket_alignment_control.py` -- local-pocket-alignment control: re-derives T/Theta in a frozen local-pocket reference frame (instead of whole-protein) to test whether the rotational signal is a whole-protein-frame artifact. Includes the `frame_relative_theta` function exercised by the test suite below.
- `36_local_pocket_control_evaluation.py` -- evaluates the local-pocket control against the whole-protein result (AUROC, agreement statistics).
- `37_pocket_frame_diagnostic_audit.py` / `38_pocket_frame_agreement_analysis.py` -- full-cohort frame-divergence audit and Bland-Altman-style agreement analysis between the whole-protein and local-pocket frames.
- `39_plot_prpd_pocket_control_figure.py` -- generates the paper's 3-panel PRPD/pocket-control main-text figure.
- `40_cross_ligand_20ns_coordinate_check.py` -- tests the coordinate-ablation finding and a short-gap (2,7]ns-to-(14,20]ns forecast on two other ligands (Beauvericin, Verapamil) using recovered 20 ns trajectories, entirely independent of the primary milbemycin cohort.
- `test_synthetic_frame_relative_theta.py` -- synthetic unit-test suite (six cases, no trajectory files needed) validating `frame_relative_theta` against known, programmatically-applied rigid-body motions. This one **is** runnable standalone from this snapshot:
  ```bash
  cd posegate/scripts/cad_paper_analysis
  python3 test_synthetic_frame_relative_theta.py
  ```

Full narrative writeup of what these established (the pose-decomposition mechanism, the local-pocket-alignment control result, and the synthetic-test findings) is in `yau_aj.tex`'s "Pose decomposition identifies rigid-body reorientation as the dominant signal" subsection and Supplementary Material, in the source repo (not included in this snapshot -- see the top-level note).

## What "source of truth" means here

There is no external experimental (binding/activity) dataset. The label is
self-generated from the project's own MD trajectories: pose is "retained"
if median protein-relative, PBC-corrected ligand RMSD over the final
(70,100] ns window is < 3 Å, with complete coverage required. Full
provenance and the two-cohort distinction (locked vs. expanded) are
documented in `pipeline/training/outputs/yau_recapture_screen/cad_states/same_trajectory_source_of_truth/README.md`
and `.../expanded_source_of_truth/master/README.md` in this snapshot.

See also, in the source repo: `doc/YAU_competition_2026/yau_aj_trajectory_triage_DATASETS.md`.
