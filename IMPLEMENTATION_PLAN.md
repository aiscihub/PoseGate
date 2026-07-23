# PoseGate-MD implementation plan

## Product boundary

PoseGate-MD 0.1 is a local command-line sidecar that measures a running
trajectory prefix, applies a frozen same-trajectory policy, and writes one
immutable shadow record before late frames exist. It does not modify or
terminate the simulation.

The numbered manuscript scripts remain unchanged as the research provenance
layer. PoseGate-MD is the reusable implementation for new trajectories.

## Phases

### 1. Freeze the scientific contract

- Freeze the 5 ns `(3,5]` and 20 ns `(18,20]` checkpoint specifications.
- Preserve the validated per-atom minimum-image representation relative to
  the protein Cα alignment centroid.
- Keep corrected RMSD primary and centroid displacement supporting.
- Treat frame-relative reorientation and internal deformation as diagnostics.
- Freeze the `(70,100]` ns outcome and 3 Å retention threshold.
- Version every policy and record its training/model provenance.

Changing to bond-connectivity molecule reconstruction would define a new
measurement and policy version; it must not silently alter the bundled v1
policies.

### 2. Extract the measurement library

- Implement strict topology, selection, cadence, checkpoint, and PBC checks.
- Extract protein alignment and ligand coordinate correction.
- Compute corrected RMSD, centroid displacement, frame-relative
  reorientation, and internal deformation.
- Match synthetic invariances and the existing scripts exactly.

### 3. Deliver shadow mode

- Provide `inspect`, `measure`, `shadow`, `watch`, and `validate`.
- Use portable frozen policy coefficients rather than version-sensitive
  model deserialization.
- Seal configuration, topology, coordinate-prefix, and optional OpenMM
  checkpoint hashes.
- Make records immutable and idempotent by run and policy.
- Refuse a genuine shadow capture after the late-outcome embargo begins.
- Produce JSON first; add concise HTML reporting after the record schema is
  stable.

### 4. Add OpenMM sidecar adapters

- Discover topology, DCD, log, and checkpoint files from a run directory.
- Monitor reporter-safe checkpoint boundaries.
- Report run state and GPU-independent measurement progress.
- Never modify the live simulation in shadow mode.

### 5. Add reversible scheduling

- Pause and resume only at aligned DCD/checkpoint boundaries.
- Reprioritize queued continuations and use otherwise idle GPU periods.
- Keep permanent early termination disabled until prospective false-stop and
  resource-savings evidence is sufficient.

## Current development status

- Strict configuration schema: implemented.
- Portable frozen 5 ns and 20 ns policies: implemented.
- Corrected geometry library: implemented.
- Immutable shadow records: implemented.
- `inspect`, `measure`, `shadow`, `watch`, and `validate`: implemented as
  direct-file commands.
- Synthetic/configuration/immutability tests: implemented.
- Real `MDR1_CRYNH:pocket12` 20 ns golden parity: passing exactly for RMSD,
  centroid displacement, coordinate-prefix hash, and score.
- Wheel packaging with bundled policy configurations: verified.

## Remaining release gates

- Add 3–5 redistributable small golden fixtures rather than relying on a
  private full trajectory.
- Add OpenMM run-directory discovery and status reporting.
- Add HTML report generation.
- Add continuous integration across supported Python versions.
- Confirm package metadata and release documentation remain consistent with
  the repository's MIT license.
- Add contributor documentation, security guidance for untrusted files, and
  a tagged archival release.
