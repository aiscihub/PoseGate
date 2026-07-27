# Prospective shadow test v2

This directory replaces the live **decision protocol**, not the retrospective
scientific result. The core observation is deliberately simple: measure one
trajectory's early protein-relative ligand-pose deviation and ask whether it is
large enough to justify an early-stop candidate.

V2 has **no runtime logistic regression and no retention probability**. It
reports the measured 5 ns RMSD, one frozen threshold in ångströms, the signed
threshold margin, and one of two actions:

- `stop_candidate`: the early RMSD exceeds the frozen cutoff;
- `do_not_stop`: the cutoff was not exceeded. This is not a prediction that the
  trajectory will retain its pose.

Every shadow trajectory continues to 100 ns. The action is counterfactual until
the outcome is known.

## What is fixed

The implementation enforces all of the following:

- one primary checkpoint: **5 ns**;
- exactly ten saved frames for `(3,5]` ns: zero-based indices `15..24`, nominal
  times `3.2, 3.4, ..., 5.0` ns;
- whole-ligand periodic imaging, followed by protein-CA Kabsch alignment;
- one direct RMSD threshold selected on a strict remeasurement of the locked
  development cohort under a predeclared false-stop cap;
- complete prospective enrollment before production;
- an immutable arm record created only while the enrolled production DCD is
  absent;
- capture no later than 6 ns and before any `(70,100]` outcome frame exists;
- exactly 150 outcome frames: zero-based indices `350..499`, nominal times
  `70.2, 70.4, ..., 100.0` ns;
- stable double reads of both checkpoint and outcome coordinates;
- one unique launch ID, one decision, and one outcome per enrolled launch;
- immutable, SHA-256-sealed policy, enrollment, arm, event, capture,
  disposition, and validation records;
- runtime code hashes and software versions checked before capture and
  validation.

## Required workflow

Run these steps in order from a clean, committed implementation in the actual
project repository and the same Python/MDAnalysis environment used for live
capture.

```bash
# 1. Re-measure all 52 locked development trajectories with exact frame indices.
python yau_competition/scripts/cad_paper_analysis/54_shadow_reimplementation_audit_v2.py

# 2. Select a direct threshold. Declare this cap before running the command.
#    Zero is the clearest safety-oriented development rule; one may be used only
#    when explicitly justified as the accepted development tradeoff.
python yau_competition/scripts/cad_paper_analysis/55_select_shadow_threshold_v2.py \
  --max-false-stops 0

# 3. Commit the complete v2 runtime implementation. Then freeze the policy.
python yau_competition/scripts/cad_paper_analysis/52_shadow_monitor_v2.py \
  freeze-policy

# 4. Replace every placeholder in enrollment/enrollment_template_v2.csv, list the entire
#    prospective cohort before any feature is measured, and freeze enrollment.
python yau_competition/scripts/cad_paper_analysis/52_shadow_monitor_v2.py \
  freeze-enrollment \
  --csv /path/to/enrollment_v2.csv \
  --cohort-label yau-shadow-v2-pilot \
  --minimum-primary-launches 30 \
  --minimum-distinct-proteins 10
```

Start each monitor **before** starting its enrolled MD launch. The monitor command
blocks while waiting, so run it in a separate terminal, scheduler sidecar, or
supervised background process. Do not launch production until the immutable arm
JSON and sidecar exist.

```bash
PYTHON_BIN=/path/to/frozen/environment/bin/python \
  yau_competition/scripts/cad_paper_analysis/run_shadow_monitor_v2.sh LAUNCH_ID

# In a second shell/scheduler step, after confirming the arm record exists:
# start the enrolled MD launch that writes the declared DCD path.
```

After the same append-only trajectory reaches at least 500 saved frames:

```bash
python yau_competition/scripts/cad_paper_analysis/53_validate_shadow_outcome_v2.py \
  --launch-id LAUNCH_ID

python yau_competition/scripts/cad_paper_analysis/56_summarize_shadow_v2.py
```

When a frozen launch cannot produce an outcome for a documented technical
reason, seal that fact rather than deleting or replacing the launch:

```bash
python yau_competition/scripts/cad_paper_analysis/52_shadow_monitor_v2.py \
  record-disposition \
  --launch-id LAUNCH_ID \
  --code failed_after_capture_before_outcome \
  --reason "GPU node failed and the append-only DCD could not be resumed"
```

The disposition records whether partial or complete outcome data were already
available. A completed 100 ns trajectory cannot be excluded as
`failed_after_capture_before_outcome`; it must be validated. Technical
non-outcomes remain in cohort flow and never authorize post-freeze replacement.

The delivered bundle contains no pre-frozen v2 policy. The strict 52-trajectory
remeasurement and threshold selection must run against the actual raw development
trajectories before a scientifically valid cutoff can be sealed.

## Primary interpretation

| Captured action | Late retained | V2 category | Operational meaning |
|---|---:|---|---|
| `stop_candidate` | no | `successful_stop` | compute could have been saved |
| `stop_candidate` | yes | `false_stop` | scientifically costly premature stop |
| `do_not_stop` | yes | `safe_continue` | retained trajectory was preserved |
| `do_not_stop` | no | `false_continue` | compute was spent, but no retained run was lost |

Headline operating metrics are false-stop risk, stop precision, failure capture,
and counterfactual compute saving. The summary also reports descriptive
Spearman association between X5 and late RMSD and raw-X5 AUROC for late
non-retention. Those test the observation itself independently of the chosen
threshold.

A completed cohort with zero validated `stop_candidate` calls does not test the
stop branch. The summary labels that state explicitly rather than calling it a
validated operating policy.

## Staged use

V2 is a **shadow** protocol. It can prospectively evaluate the decision rule,
but it cannot demonstrate realized compute savings because every run continues
to 100 ns. Enforced stopping belongs in a later, separately frozen action pilot.
That pilot should preserve an audit sample of stop candidates to 100 ns or use a
randomized continuation design; otherwise false stops become unobservable.

## V1 status

All existing 5 ns and 20 ns logistic captures remain immutable **legacy v1**
records. They are not converted, relabeled, or pooled with v2. The two v1
checkpoint calls on one trajectory are correlated evaluations of one outcome,
not two independent prospective cases.

See `PROTOCOL_V2.md`, `ANALYSIS_PLAN_V2.md`, and `DEPLOYMENT_GATES_V2.md`.
