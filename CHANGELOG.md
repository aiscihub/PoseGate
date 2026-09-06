# Changelog

## 0.1.0.dev0

- Define a strict, versioned scientific configuration contract.
- Add portable 5 ns and 20 ns frozen standardized-logistic policies.
- Extract protein-relative PBC correction and ligand pose geometry.
- Add inspect, measure, shadow, watch, and validate command foundations.
- Add immutable, hash-audited shadow records.
- Add synthetic geometry, policy parity, and optional real-trajectory tests.
- Record schema 1.1: per-frame prefix series, structured run identity, named
  quality-control checks, and a measured ligand-diameter domain profile.
- Add the one-sided `threshold_rule` policy type alongside
  `standardized_logistic`. It reports a signed stop margin in angstrom and no
  score, and returns `DEFER` rather than `POSE_RETAINED` when the stop
  condition is not met.
- Add an optional `applicability` configuration block; a policy refuses to
  measure a ligand outside its declared domain instead of extrapolating.
- Add `posegate report`, which renders one sealed record as a self-contained
  HTML page with no external requests.
- Report the audit category (`successful_stop`, `false_stop`, `safe_continue`,
  `false_continue`, and the two deferral cases) from `validate`, and re-apply
  the sealed policy to the sealed measurement to detect an edited forecast.
- Report the feature-window frame count on revalidation, so a window boundary
  that moved between two readings of the same coordinates is visible.
- Fix the real-trajectory parity test, which compared PoseGate's recommendation
  vocabulary against the legacy one and so could only ever pass on a continue.

Bundled v1 policy configurations are unchanged and hash exactly as before.
