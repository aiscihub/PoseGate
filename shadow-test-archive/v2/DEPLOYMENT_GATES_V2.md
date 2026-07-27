# Shadow-to-action deployment gates

The shadow test and an enforced stopping workflow answer different questions.
This file prevents a successful software demonstration from being mistaken for
permission to terminate production trajectories.

## Gate 0: development lock

Required before any v2 launch:

- strict 52-trajectory remeasurement completed with 10/150 exact frames;
- one direct X5 threshold selected under the declared development false-stop
  cap;
- policy and runtime hashes frozen;
- complete prospective enrollment frozen;
- no v1 outcome used to choose the v2 threshold.

Failure of any item blocks prospective capture.

## Gate 1: engineering integrity

Required for every counted launch:

- immutable pre-production arm record;
- production-start event bound to the same topology and DCD inode;
- capture by 6 ns with stable double read and no outcome-frame access;
- complete 100 ns endpoint with stable double read;
- capture-to-validation feature parity;
- no unresolved provenance violation.

A launch failing this gate remains in cohort flow as a technical disposition.

## Gate 2: completed prospective shadow pilot

Required before discussing operational validity:

- every enrolled primary launch is validated or has a sealed disposition;
- required launch and protein breadth are met;
- the stop branch is exercised by at least one validated stop candidate;
- false-stop count/rate, stop precision, failure capture, and counterfactual
  compute saving are reported with raw denominators;
- continuous X5-versus-late diagnostics are reported;
- dispositions after partial outcome availability are disclosed;
- no threshold, endpoint, or enrollment change occurred.

Passing this gate supports a **prospective pilot** claim, not guaranteed risk
control and not realized compute saving.

## Gate 3: decision to create an action pilot

This decision must be written before any enforced stop. It should specify:

- an acceptable false-stop target and why it is scientifically tolerable;
- minimum retained outcomes and protein breadth;
- the exact frozen threshold inherited from the completed shadow cohort;
- a randomized or predeclared audit-continuation fraction among stop candidates;
- how stopped and audited trajectories enter analysis;
- a rule for pausing the action pilot if false stops exceed the declared bound;
- no reuse of the action cohort to retune the same threshold.

The threshold must not be lowered simply because the shadow cohort produced too
few stop calls. That would create a new policy requiring a new prospective
cohort.

## Gate 4: realized compute-saving claim

Actual savings can be claimed only from an enforced action pilot. Because most
stopped trajectories lack their 70–100 ns counterfactual outcome, the audit
continuation design is essential for estimating false-stop risk. Report:

- wall-clock/GPU-hour and simulated-nanosecond savings;
- audit fraction and selection mechanism;
- observed false stops among audited stop candidates;
- uncertainty induced by unaudited stopped runs;
- any impact on the downstream scientific screen.

## Addendum, 2026-07-25: declared single-case scope reduction

The v2 threshold (`stop_threshold_A = 2.7634173197224974`) was frozen
2026-07-25. Reaching Gate 0's declared ≥30-launch/≥10-protein enrollment
requires launching ~28 new `replica_2` runs (only 2 of 80 known protein/pocket
slots have zero production data anywhere), at an estimated ~18h/launch on the
one available GPU — a ~3-week continuous-compute commitment. That commitment
was deliberately deferred, undated, not started, and no cohort was generated
for the single-case demonstration below.

In its place, `CDR2_CANAL__Milbemycin__pocket14__replica_1` was designated a
**declared single-case (n=1) scope reduction**: a single-trajectory,
pre-outcome implementation pilot (equivalently, a one-case prospective shadow
demonstration). This does **not** relax Gate 1. pocket14's production DCD
started 2026-07-24, before the v2 policy existed (frozen 2026-07-25), so
`production_dcd_must_be_absent_when_armed` is already violated and cannot be
satisfied without discarding that DCD — which was deliberately not done, to
avoid destroying already-completed GPU work. Consequently this case can never
become a validated prospective v2 launch under Gate 1, at any cohort size.

This launch also already had an immutable v1 checkpoint record (5 ns, scored
2026-07-24, decision `continue`) captured before its outcome was known. That
gives it two distinct, non-conflated notions of "prospective":

- **prospective with respect to the late outcome**: yes — the v2 pilot record
  below was sealed while only 20 ns (100 frames) of the trajectory existed,
  well before the (70,100] ns outcome window;
- **prospectively enrolled under v2 before simulation started**: no — the
  trajectory, and its v1 early result, existed before the v2 design was
  finalized.

Rather than weakening or bypassing the normal arm/capture path, a separate,
explicitly-labeled **pilot-import path**
(`yau_competition/scripts/cad_paper_analysis/57_shadow_pilot_import_v2.py`)
was added. It (1) loads the existing immutable v1 5 ns capture and verifies it
used exactly the frozen 25-frame prefix with no late-outcome access, (2)
verifies the topology sha256 against that capture, (3) re-derives the v1
capture script's own coordinate digest for frames 0..24 against the *current*
trajectory and requires an exact match — proving those frames are unchanged
since v1 captured them, (4) recomputes the strict v2 X5 feature independently
from frames 15..24 (not copied from v1, which used an 11-frame window and
different periodic-image handling), and (5) applies the frozen threshold and
seals a new, distinctly-labeled record before sealing if the (70,100] ns
outcome window is already exposed (in which case it would instead be labeled
`frozen_rule_replay_on_prospective_prefix`, not a clean pre-outcome pilot —
not the case here).

Sealed result (2026-07-25T19:49:06Z):
`yau_competition/prospective_shadow/v2/pilot_import/CDR2_CANAL__Milbemycin__pocket14__replica_1__pilot_import_v2.json`
— `capture_mode: existing_prefix_preoutcome_pilot`, `outcome_prospective: true`,
`launch_prospectively_enrolled_under_v2: false`, strict X5 = 1.292505092340566 Å,
signed margin = -1.4709122273819313 Å, action = `do_not_stop`. Outcome
(late (70,100] ns median RMSD, final category) is pending pocket14's resumption
to 100 ns.

This record must never be pooled with, or presented as, a Gate-2 multi-launch
pilot result, and does not support any claim of false-stop rate, stop
precision, sensitivity/specificity/AUROC, or threshold generalization — see
the record's own `excluded_claims` and `n1_evaluation_note` fields. It
demonstrates only that the frozen measurement code and rule run correctly
against one real, previously-unobserved-outcome launch.

pocket14 itself remains open under v1 (see `MIGRATION_V1_TO_V2.md`): its 5 ns
and 20 ns v1 checkpoint captures already exist and it is expected to be
resumed to 100 ns and validated with the legacy `51_validate_shadow_outcome.py`.
Once it reaches 100 ns, this pilot-import record should also be finalized with
the late outcome (a companion validation step, not yet built).

## Addendum, 2026-07-27: second single-case pilot, CIMG_06197:pocket15 — no v1 precedent, first v2 `stop_candidate`

`CIMG_06197:pocket15` is the next (and, after re-checking the actual remote
trajectory/log state rather than the stale
`mutation_pipeline/scripts/incomplete_100ns_rerun_manifest.csv`, the *only*
remaining) genuinely incomplete launch from the manifest's original 8-pocket
list. Unlike pocket14, it was never run under the v1 monitor, so there is no
existing v1 5 ns capture to cross-verify a pilot-import against. A new script,
`yau_competition/scripts/cad_paper_analysis/60_shadow_fresh_prefix_pilot_v2.py`,
implements the same declared single-case scope reduction without that
cross-check step: it measures the strict v2 X5 feature directly off whatever
prefix already exists (via `shadow_v2_core.find_production_inputs`/
`measure_prefix`), applies the frozen threshold, and seals a
`record_type: v2_fresh_prefix_pilot` / `capture_mode: fresh_prefix_preoutcome_pilot`
record with the same `excluded_claims`/`n1_evaluation_note` discipline as
`57_shadow_pilot_import_v2.py`.

**A real process error was caught before any GPU time was lost.** The script
was first run against the *original* truncated trajectory at
`simulation_100ns_md/CIMG_06197/simulation_explicit/pocket15/replica_1`
(paused at 39.8 ns since 2026-07-09), sealing a record referencing that exact
trajectory (X5 = 1.404 Å, `do_not_stop`). Immediately afterward,
`run_rerun_with_watchdog.sh CIMG_06197 pocket15 Milbemycin 100` was launched
expecting it to resume that same checkpoint — but that script's own header
says "clean 0-100 ns rerun under `simulation_100ns_md_reruns/`": it only
stages the prepared input PDB/SDF into a new directory tree and always starts
a fresh trajectory from step 0, never touching the original checkpoint/DCD.
This was caught while the new run was still in minimisation (zero DCD/
checkpoint frames written, no compute lost). The original sealed record was
left untouched (immutable by design, and it correctly described the original
trajectory at seal time) with a plain non-immutable sidecar,
`..._fresh_prefix_pilot_v2.SUPERSEDED_NOTE.json`, documenting that its target
trajectory was abandoned and it must never be finalized against the new one.
The new independent trajectory under `simulation_100ns_md_reruns/CIMG_06197/
simulation_explicit/pocket15/` was then let run (matching the user's explicit
choice to generate fresh rerun output there, consistent with how the other
manifest pockets — AFR1:pocket9, MDR1_CRYNH:pocket12, and the majority of the
manifest's remaining entries, which turned out to already be complete — were
actually handled).

**Filename collision bug fixed the same session:** `60_shadow_fresh_prefix_pilot_v2.py`
originally named its output only by protein/pocket/replica, which collided
with the orphaned record above once run a second time against the correct
(reruns-tree) trajectory. Fixed by suffixing the filename with the first 8
hex characters of `trajectory_sha256` (`srctraj-<hash>`), so two genuinely
different trajectories for the same pocket can never collide or silently
overwrite each other's sealed record.

**Both checkpoints are now sealed against the correct, single trajectory**
(`simulation_100ns_md_reruns/CIMG_06197/simulation_explicit/pocket15/replica_1`,
`trajectory_sha256` starting `e9c08dec`):

- 20 ns v1 (legacy, run alongside v2 at the user's request since it's "more
  useful" than v2 alone): corrected RMSD (18,20] = 4.303 Å, centroid
  displacement = 3.644 Å, decision **`stop`** — the first `stop` this v1 20 ns
  model has produced in the project (all prior calls were `continue`). Shadow
  only; the run correctly continued regardless
  (`captures/CIMG_06197__Milbemycin__pocket15__replica_1__same_trajectory_checkpoint_20ns_v1.json`).
- 5 ns v2 (sealed at 32.8 ns available, well pre-70ns embargo): strict X5 =
  4.169 Å vs. frozen threshold 2.763 Å, signed margin = **+1.405 Å**,
  **`stop_candidate`** — the first `stop_candidate` the frozen v2 rule has
  produced in the project (pocket14 was `do_not_stop`)
  (`pilot_import/CIMG_06197__Milbemycin__pocket15__replica_1__srctraj-e9c08dec__fresh_prefix_pilot_v2.json`).

v1 and v2 agree directionally on this trajectory: both flag early divergence
from the docked pose. As with pocket14, this remains n=1 and must not be
pooled into any rate claim. **Still open:** no finalizer exists yet for
`v2_fresh_prefix_pilot`-type records once the trajectory reaches 100 ns —
`58_shadow_pilot_import_outcome_v2.py` expects an `imported_v1_capture` field
this record type doesn't have, so it needs to be extended or a small
analogous finalizer written before this case can be closed out.
