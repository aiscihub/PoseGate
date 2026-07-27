> **Legacy v1.** Do not use these logistic 5 ns/20 ns monitors for new captures. The corrected one-checkpoint, direct-threshold protocol is in [`v2/README.md`](v2/README.md). Existing records remain immutable historical evidence.

# Prospective shadow captures

This directory stores timestamped, one-shot decisions issued from a new 100 ns
trajectory's own `(18,20] ns` prefix before its `(70,100] ns` outcome exists.
Every trajectory remains a shadow run and must continue to 100 ns regardless
of the recorded `stop` or `continue` decision.

The workflow is implemented by
`yau_competition/scripts/cad_paper_analysis/prospective_shadow_monitor.py`.

## Frozen specification

- Checkpoint lock: `same_trajectory_checkpoint_20ns_v1`
- Feature window: `(18,20] ns`
- Features: corrected protein-relative pose RMSD and centroid displacement
- Model: all-development-data StandardScaler plus class-weighted L2 logistic regression
- Threshold: `0.5`
- Score: uncalibrated retention score, not an absolute probability
- Outcome embargo: the prediction must be written before 70 ns data exist

## Records

- `frozen_model/`: read-only model and manifest with source hashes
- `events/`: immutable production-start detection records
- `captures/`: immutable per-candidate prediction JSON files
- `shadow_predictions.csv`: regenerable summary index
- `monitor_logs/`: watchdog logs
- `state/`: process lock and PID files

Captures are same-trajectory prospective monitoring evidence. They are not
independent-replica validation, external validation, binding evidence, or a
guarantee of biological activity.
