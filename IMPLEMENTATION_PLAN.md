# Implementation status

## Released contract

Version 0.2.0 implements whole-ligand imaging, exact nominal windows, the
one-sided X5 rule, and same-trajectory monitoring and shadow triage. It adds
current results, an 81-trajectory manifest, provenance mapping, checksums,
and compact saved-frame fixtures.

Legacy v1 per-atom/logistic configurations remain hash-identical test fixtures,
not installed policies. Obsolete research/capture archives remain in Git history.
The current pilot configuration is a 5 ns checkpoint with a 20 ns outcome,
not the legacy 20 ns checkpoint.

## Implemented

- Configuration, selection, cadence, periodic-box, and prefix validation.
- Whole-ligand correction, protein alignment, corrected RMSD, and controls.
- Strict one-sided threshold; unflagged trajectories remain unresolved.
- `inspect`, `measure`, `shadow`, `watch`, `validate`, and HTML `report`.
- Immutable records, revalidation, synthetic and saved-frame regression tests.
- Current manuscript data, mapping, numerical checks, and release checksums.

## Outside this release

- Original prospective enrollment/arming/deadline infrastructure.
- OpenMM directory discovery, GPU/job monitoring, and a registry-wide console.
- Automatic pause/resume, scheduling, or permanent termination.
- Distribution of complete raw trajectories, which are too large to bundle.

These are separate development tasks, not claims of this release. See `RELEASE.md`.
