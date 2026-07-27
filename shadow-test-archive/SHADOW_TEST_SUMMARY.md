> **Legacy v1 results.** Report these separately from v2. At the launch level, two completed trajectories were false continues; one contributed both a 5 ns and 20 ns checkpoint call. See [`v2/MIGRATION_V1_TO_V2.md`](v2/MIGRATION_V1_TO_V2.md).

# Prospective shadow test — summary (as of 2026-07-25)

## What this is

The manuscript (`doc/YAU_competition_2026/yau_aj.tex`) makes a *retrospective*
claim: early-window ligand-pose geometry (5 ns or 20 ns into a 100 ns MD run)
predicts whether the pose is still "retained" at the end (70-100 ns). That
claim was fit and validated entirely on trajectories that had already finished
— so the model has, in a sense, already seen every outcome it was scored
against (via cross-validation, not directly, but still on completed data).

The shadow test is the check on whether that retrospective rule actually
works as a **live forecast**: freeze the exact fitted rule, then apply it to
brand-new 100 ns reruns *while they are still running*, writing down the
prediction before the trajectory reaches the outcome window. This is the only
part of the project where a prediction is committed before the answer exists.

## How it works (mechanics)

1. **Freeze once.** A `StandardScaler` + `LogisticRegression` is fit exactly
   once per checkpoint (5 ns, 20 ns) on the locked 52-row retrospective
   cohort, then written to a read-only `.joblib` file with a SHA-256 manifest.
   Every later shadow prediction reuses that exact object — never refit.
2. **Measure only the past.** For a new trajectory, the monitor reads the DCD
   only up through the checkpoint frame and refuses to run
   (`require_pre70=True`) if the file already exposes ≥70 ns — it is
   structurally unable to see or be influenced by the outcome window.
3. **Score and lock the record.** The frozen model's `predict_proba` gives
   `uncalibrated_retention_score`; `>= 0.5` → "continue", else "stop". The
   result is written with `os.O_EXCL` + `chmod 0o444` the instant it's
   computed, so it cannot be edited after the fact once the real outcome
   later becomes known.
4. **Validate later.** Once the trajectory reaches 100 ns, a separate script
   computes the true (70,100] ns outcome and compares it to the frozen
   prediction.

Full technical detail (feature math, exact hyperparameters, worked score
calculations): `SCORE_METHODOLOGY.md` in this directory.

## Results so far

| Trajectory | Checkpoint | Score | Decision | Actual outcome | Correct? |
|---|---|---|---|---|---|
| MDR1_CRYNH:pocket1 | 20 ns | 0.687 | continue | 3.575 Å, non-retained | No — clear miss |
| MDR1_CRYNH:pocket12 | 5 ns | 0.717 | continue | 3.047 Å, non-retained | No — boundary miss (0.047 Å past the 3.0 Å cutoff) |
| MDR1_CRYNH:pocket12 | 20 ns | 0.720 | continue | 3.047 Å, non-retained | No — same boundary miss |
| CDR2_CANAL:pocket14 | 5 ns | 0.761 | continue | 6.081 Å, non-retained (resolved 2026-07-26) | No — false continue |
| CDR2_CANAL:pocket14 | 20 ns | 0.685 | continue | 6.081 Å, non-retained (resolved 2026-07-26) | No — false continue |
| CIMG_06197:pocket15 | 20 ns | 0.314 | **stop** | not yet — trajectory still running to 100 ns | pending |

Three trajectories are now fully validated (MDR1_CRYNH:pocket1,
MDR1_CRYNH:pocket12, and CDR2_CANAL:pocket14, resolved 2026-07-26 at late
median RMSD 6.081 Å). All three were `continue` calls followed by
non-retention, so all are false continues at the launch level; one trajectory
(pocket12) contributed both a 5 ns and a 20 ns checkpoint record and must not
be counted twice as an independent prospective case. Misses ranged from clear
(3.575 Å, 6.081 Å) to a narrow boundary miss (0.047 Å past the 3.0 Å cutoff).
A fourth trajectory (CIMG_06197:pocket15) has a 20 ns checkpoint captured and
locked but no outcome yet — notably the **first `stop` decision** this v1
20 ns model has produced in the project (all prior calls were `continue`).
The strict label remains non-retained for the 3.047 Å case; threshold-sensitivity
context does not relabel it as correct.

This is a small-n live demonstration, not a validated error rate. It exists so
the paper can say the checkpoint rule was tested prospectively at least once,
not only cross-validated retrospectively.

## Why the forecasting rule is derived from the locked 52-row cohort, not the 73-row expanded cohort

Short answer: **the 73-row cohort is real and larger, but "locked" means
frozen *before* being looked at again — and the 73-row cohort has already
been looked at (used to compute the expanded-cohort audit numbers already in
the draft). Fitting the forecasting rule on it now would not be a live
forecast anymore; it would be curve-fitting with the benefit of hindsight.**

The longer reasoning (from `doc/YAU_competition_2026/cohort_expansion_80_status.md`
and `yau_aj_trajectory_triage_DATASETS.md`, both already in the repo):

1. **Freezing is what makes the reported numbers trustworthy.** If the
   "locked" cohort could just be enlarged whenever new trajectories finish,
   the reported AUROC/CI numbers would be quietly re-tuned by whichever data
   happens to exist at the moment — that's the classic "undisclosed multiple
   looks / optional stopping" statistical hazard. 52 was never claimed to be
   a ceiling, just a deliberate freeze point.
2. **Most of the growth isn't new information for the thing that matters.**
   Of the 21 trajectories added going from 52 → 73, only 1 protein
   (`CIMG_06197`, 5 rows) is genuinely new; the rest are additional
   pockets/replicas on the *same* 15 proteins already in the locked cohort.
   The validation scheme is leave-one-**protein**-out (LOPO), so protein
   count — not trajectory count — is what buys real statistical power. Class
   balance also didn't move much (13/39 = 25.0% retained at n=52 vs. 17/56 =
   23.3% at n=73), so the extra rows mostly add redundant within-protein
   signal, not new generalization evidence.
3. **The 73-row cohort has a lighter audit trail.** The locked 52 went
   through a heavier QC/redundancy chain (outcome-redundancy audit,
   SHA-256-pinned source manifest, legacy-number reconciliation — scripts
   25-27). The expanded cohort's pipeline (script 28-31) is simpler and
   hasn't been through that same provenance hardening yet. Promoting it to
   primary status as-is would quietly weaken the provenance story even if the
   headline numbers looked fine.
4. **It's already been "used."** The 73-row cohort is currently cited in the
   manuscript for a secondary check — the model-complexity audit (comparing
   corrected RMSD vs. nested logistic/spline/random-forest/gradient-boosting
   models). Once a cohort has been looked at and used to compute *any*
   reported number, it can't turn around and also be the cohort a forecasting
   rule gets "prospectively" validated against — that would be circular.

This isn't a permanent decision to ignore the bigger dataset — see
`cohort_expansion_80_status.md`'s decision log: once the primary cohort
finishes growing to n≈80, the plan (already written down, not yet executed)
is to run that n≈80 set through the *same* heavier audit chain the n=52 lock
used, and present it as an explicit **second locked cohort ("locked v2")**
that supersedes n=52 for the primary claim — reported side by side with the
original n=52 result, not silently overwriting it. The 73-row set as it
exists today is an intermediate/exploratory checkpoint on the way there, not
the final answer.

## Code index (packed for inspection)

All of the following are new/untracked (not yet committed) and are bundled
into `/Volumes/SharedFolder/git/valleyfevermutation/yau_competition/prospective_shadow/shadow_test_code_bundle.tar.gz`
for easy offline review:

**Monitors (score a running trajectory without touching its future frames)**
- `yau_competition/scripts/cad_paper_analysis/prospective_shadow_monitor.py` — 20 ns lock
- `yau_competition/scripts/cad_paper_analysis/50_prospective_shadow_monitor_5ns.py` — 5 ns lock

**Outcome validation (run only after a trajectory reaches 100 ns)**
- `yau_competition/scripts/cad_paper_analysis/51_validate_shadow_outcome.py`

**Orchestration (pause a live MD run cleanly at the checkpoint frame, keep the GPU/monitor decoupled)**
- `yau_competition/scripts/cad_paper_analysis/pause_shadow_run_at_20ns.sh`
- `yau_competition/scripts/cad_paper_analysis/pause_shadow_run_at_checkpoint.sh`
- `yau_competition/scripts/cad_paper_analysis/run_shadow_monitor_20ns_with_watchdog.sh`
- `yau_competition/scripts/cad_paper_analysis/run_shadow_monitor_with_watchdog.sh`

**Cohort-building (for context — what "locked 52" vs. "expanded 73" actually means)**
- `yau_competition/scripts/cad_paper_analysis/25_build_same_trajectory_source_of_truth.py` — builds the locked 52-row cohort
- `yau_competition/scripts/cad_paper_analysis/31_build_expanded_source_of_truth.py` — builds the 73-row expanded cohort
- `yau_competition/scripts/cad_paper_analysis/run_expanded_source_of_truth_pipeline.py` — expanded-cohort pipeline entry point

**State/records (this directory)**
- `frozen_model/*.manifest.json` — hash-audited frozen scaler+model specs, one per checkpoint lock
- `captures/*.json` — the immutable predictions themselves
- `outcome_validation_*.json` — validated outcomes for trajectories that reached 100 ns
- `shadow_predictions.csv`, `shadow_predictions_5ns.csv` — flat summary tables of all captures per lock
- `README.md`, `SCORE_METHODOLOGY.md`, `SUMMARY_MDR1_CRYNH_pocket12.md` — prior documentation
