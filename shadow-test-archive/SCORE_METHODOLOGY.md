> **Legacy v1 methodology.** The uncalibrated logistic score is retired for new work. See [`v2/PROTOCOL_V2.md`](v2/PROTOCOL_V2.md) for the direct-RMSD policy.

# How the shadow-test retention score is calculated

Reference document for `uncalibrated_retention_score` as it appears in every
`captures/*.json` prediction record. Covers both checkpoint locks currently in
use: `same_trajectory_checkpoint_5ns_v1` and `same_trajectory_checkpoint_20ns_v1`.

Source: `yau_competition/scripts/cad_paper_analysis/prospective_shadow_monitor.py`
(20 ns lock) and `50_prospective_shadow_monitor_5ns.py` (5 ns lock).

## 1. Geometric feature extraction (`measure_prefix`)

Input: the trajectory's own topology PDB and DCD, read only through the
checkpoint frame (never beyond — this is what makes the capture prospective).

For each frame up to the checkpoint:

1. Select protein CA atoms and the ligand (`resname UNK`, heavy atoms only).
2. Center the frame on the CA centroid: `protein_center = mean(ca_positions)`.
3. Express the ligand position relative to that centroid, minimum-image
   wrapped against the periodic box (`MDAnalysis.lib.distances.minimize_vectors`)
   — this removes protein whole-molecule translation and PBC image jumps.
4. Kabsch-align the current frame's centered CA coordinates onto frame 0's
   centered CA coordinates (`MDAnalysis.analysis.align.rotation_matrix`),
   and apply that same rotation to the ligand-relative vector. This removes
   protein whole-molecule rotation, so what remains is ligand motion
   *relative to the protein frame*, not lab-frame drift.
5. `delta = ligand_aligned(t) - ligand_relative(frame_0)`
6. Per-frame **pose RMSD** = `sqrt(mean(sum(delta**2, axis=atoms)))` — root-mean-
   square displacement of the ligand's coordinate-corrected pose vs. its
   starting pose.
7. Per-frame **centroid displacement** = `norm(mean(delta, axis=atoms))` — the
   translation-only component of that same delta (ligand centroid drift,
   ignoring internal ligand pose changes).

These two per-frame series are then averaged over the checkpoint's **feature
window** (not the full prefix):

| Checkpoint lock | Feature window | Features used |
|---|---|---|
| `same_trajectory_checkpoint_5ns_v1` | (3, 5] ns | `same_traj_pose_rmsd_5ns_A` only |
| `same_trajectory_checkpoint_20ns_v1` | (18, 20] ns | `same_traj_pose_rmsd_18_20_A`, `same_traj_centroid_displacement_18_20_A` |

Windowing (rather than using the single last frame) matches how the same
features were computed for the 52-trajectory locked retrospective cohort —
using anything else would make the shadow score incomparable to the model it
was fit on.

A SHA-256 digest of every frame's box + CA + ligand coordinates used in the
prefix is stored (`prefix_coordinates_sha256`) so the exact input to a given
score is independently reproducible and auditable later.

## 2. The frozen model

One `StandardScaler` + `LogisticRegression` pair per checkpoint lock, fit
**once**, on all 52 rows of the locked retrospective cohort
(`same_trajectory_source_of_truth.csv`, filtered to
`eligible_locked_primary_pose_model`), against the true (70,100]ns outcome
label `late_pose_retained_locked` (13 retained / 39 non-retained across 15
proteins). Frozen to disk as a read-only `.joblib` (`chmod 0o444`) the first
time either monitor script runs; every later shadow prediction on any new
trajectory reuses that exact fitted object — the scaler and model are never
refit per trajectory, and a manifest hash-checks the model file and its
locked training table on every load to guarantee that.

```
StandardScaler(with_mean=True, with_std=True)
LogisticRegression(
    penalty="l2", C=1.0, fit_intercept=True,
    class_weight="balanced", solver="lbfgs",
    max_iter=100, tol=1e-4, random_state=7,
)
```

- 5 ns lock: single input feature (RMSD only), matching the locked 5 ns row
  of the manuscript's own checkpoint-curve sweep
  (`checkpoint_model_results.csv`).
- 20 ns lock: two input features (RMSD + centroid displacement), matching
  the manuscript's primary two-feature policy model
  (`primary_model_validation_specification.csv`).
- `class_weight="balanced"` reweights the loss for the 13/39 class
  imbalance in the training cohort.

## 3. Scoring a new trajectory

```python
x = [[measured_feature(s)_for_this_trajectory]]
scaled = frozen_scaler.transform(x)          # uses the TRAINING mean/scale, not this trajectory's
score = frozen_model.predict_proba(scaled)[0, 1]   # P(class = "retained"), per the fitted sigmoid
decision = "continue" if score >= 0.5 else "stop"
```

`predict_proba(...)[0, 1]` is the logistic model's estimated probability that
this trajectory belongs to the "retained" class, evaluated at this
trajectory's own (RMSD[, centroid displacement]) value(s), using the
decision boundary learned from the 52 locked examples. It is called
**uncalibrated** because there is no held-out calibration set — it is a raw
sigmoid output of a model fit directly on the full development cohort, not a
validated probability (see `manifest["limitations"]` in the frozen model's
own manifest file).

Decision rule is a fixed 0.5 threshold on that score — not tuned per
trajectory, not tuned after seeing any outcome.

## 4. Worked examples (values already captured)

| Trajectory | Checkpoint | Feature(s) | Score | Decision |
|---|---|---|---|---|
| MDR1_CRYNH:pocket1 | 20 ns | RMSD, centroid disp. | 0.687 | continue |
| MDR1_CRYNH:pocket12 | 5 ns | RMSD (3,5]=1.569 Å | 0.717 | continue |
| MDR1_CRYNH:pocket12 | 20 ns | RMSD (18,20]=2.052 Å, centroid disp.=1.847 Å | 0.720 | continue |
| CDR2_CANAL:pocket14 | 5 ns | RMSD (3,5] | 0.761 | continue |
| CDR2_CANAL:pocket14 | 20 ns | RMSD (18,20]=2.258 Å, centroid disp.=2.065 Å | 0.685 | continue |

In each case, `score` is exactly the fitted model's sigmoid output at that
trajectory's own measured feature value(s) — there is no other transformation
between "measured geometry" and "reported score" besides standardization
(subtract training mean, divide by training std) and the logistic function.

## 5. Why this is a genuinely prospective score, not a fitted lookup

- The scaler/model are fit **only** on the 52 locked rows, frozen before any
  shadow trajectory's checkpoint data exists, and hash-verified unchanged on
  every reuse.
- `measure_prefix` refuses to run (`require_pre70=True`) if the trajectory
  already exposes ≥70 ns of frames — the score is structurally incapable of
  seeing (or being fit on) the very outcome window it is trying to predict.
- The capture JSON is written with `os.O_EXCL` + `chmod 0o444`
  (`atomic_json(..., immutable=True)`) at the moment of scoring, so the
  prediction cannot be edited or overwritten after the fact once the
  trajectory later reaches its true (70,100]ns outcome.

## 6. Where this is not yet validated

The score answers "how far past the retained/non-retained decision boundary
does this trajectory's early-window geometry fall, under the frozen
same-trajectory model" — it says nothing about calibration accuracy or a
guaranteed false-stop rate. See `SUMMARY_MDR1_CRYNH_pocket12.md` for a
worked case where a "continue" call (score 0.720) turned out to be a
boundary miss against the true outcome (3.047 Å vs. the 3.0 Å locked
cutoff) — the score being confidently above 0.5 does not guarantee the
downstream 3 Å binary outcome lands on the same side.
