# Migration from shadow v1 to v2

## Do not edit or delete v1

Keep the existing joblib files, manifests, captures, and outcomes unchanged.
They document what was actually run. Mark them as legacy in prose and never
combine checkpoint rows with the new launch-level v2 cohort.

The completed v1 results should be described as **two completed launches, both
false continues**. One launch had both 5 ns and 20 ns records; those are
correlated checkpoint evaluations of one outcome, not two independent cases.
The pending launch remains pending. No v1 stop call has yet been validated.

## Why v2 is a new policy

The old live scripts could select 11 feature frames and 151 late frames because
they multiplied one-based frame number by noisy DCD metadata and then used
floating comparisons. The nominal 3.0 and 70.0 ns frames could therefore enter
open-left windows. V1 also used atom-wise ligand imaging, exposed an
uncalibrated logistic score as a retention forecast, did not freeze enrollment,
and recorded but did not enforce the complete extraction-policy hash.

V2 changes all of these items. It cannot inherit v1's prospective status. The
strict remeasurement audit quantifies differences between corrected values and
the old locked table before a new threshold is selected.

## Cutover checklist

1. Copy the v2 scripts and protocol directory into the repository.
2. Run `yau_competition/tests/test_shadow_v2.py`.
3. Run the strict 52-trajectory remeasurement audit.
4. Inspect old-versus-strict feature and outcome differences, especially label
   changes near 3 Å.
5. Declare the development false-stop cap before threshold selection.
6. Select the direct threshold without reading any shadow outcomes.
7. Commit the complete runtime implementation.
8. Freeze the policy and complete enrollment.
9. Start each monitor before its production DCD exists and verify the immutable
   arm record before launching MD.
10. Never change threshold or enrollment after the first capture.
11. Validate every complete endpoint; do not dispose a completed trajectory as
    missing.
12. Seal genuine technical non-outcomes with `record-disposition`; never replace
    them after inspecting a feature or outcome.
13. Report v1 and v2 separately.
14. Do not enforce stops until the shadow cohort passes its predeclared gate and
    a separately frozen action pilot is designed.
