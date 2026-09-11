# PoseGate-MD

**Same-trajectory monitoring and shadow triage.** PoseGate-MD measures
protein-relative ligand-pose change, applies a frozen one-sided rule to an
observed trajectory prefix, and records a shadow decision for later evaluation.
It does not stop simulations. The early-clean and separate-launch results do
not establish advance forecasting or replica-independent prediction.

![Measurement and evaluation](paper/current/figures/primary_idea.png)

## Manuscript release

The manuscript-aligned release is **v0.2.0**. Cite its full commit ID, not a
moving branch: `git rev-parse 'v0.2.0^{commit}'`. The annotated tag also records
that ID. A local tag is not a published hosted release.

- [Current table extracts and CSVs](paper/current/tables/) and [figure index](paper/current/figure_index.json).
- [MANUSCRIPT_MAPPING.md](MANUSCRIPT_MAPPING.md): every table and figure, scripts, inputs, outputs, and audit scope.
- [Final 81-trajectory manifest](paper/current/cohorts/final_expanded_81.json): 16 proteins, including the 52 primary trajectories.
- [Checksums](SHA256SUMS), [source provenance](paper/current/manifest.json), and [release notes](RELEASE.md).

The manuscript has not been submitted. Full manuscript and supplementary PDFs
and TeX documents are not bundled or published. Their source hashes are
recorded in the provenance manifest; table extracts and figure assets are included.

## Installation and verification

Python 3.10 or newer, from the release checkout:

```bash
python -m pip install -e '.[test,paper]'
posegate --version
python tools/release_checksums.py
python tools/verify_paper_results.py
python -m pytest -q
```

The verifier audits distributed measurements without raw trajectories. Tests
include synthetic invariances, two saved-frame fixtures, threshold equality,
configuration hashes, and immutable records. No full raw trajectory is required.

## Current measurement and Table 5 rule

Use [configs/posegate_5ns_v2.yaml](configs/posegate_5ns_v2.yaml). For each frame
independently, a common lattice shift selects one periodic image of the whole
ligand using its centroid relative to the protein C-alpha centroid. After
protein centering and alignment to the first saved frame, the same rotation is
applied to the ligand. RMSD compares corresponding ligand heavy atoms. This
requires whole saved components; it does not unwrap the protein or reconstruct
a broken ligand. PBC-naive and ligand-self-aligned RMSD are diagnostic controls,
not policy inputs.

At the nominal 0.2 ns output cadence, frame 0 is 0.2 ns. X5 averages indices
15 through 24, corresponding to (3,5] ns; L100 is the median of indices 350
through 499, corresponding to (70,100] ns. Retention means L100 < 3 angstrom;
equality is non-retention. The Table 5 rule is exactly:

```text
X5 > 2.7634173197224974 angstrom: stop candidate (shadow only)
X5 <= 2.7634173197224974 angstrom: do not stop; retention unresolved
```

No logistic regression, probability, or numerical tolerance is added to that
strict comparison. Unflagged trajectories return `DEFER`, not a retention
prediction. Prospective validation flagged nine trajectories: eight later
non-retained and one retained. Table 5 counts 855 ns of potential remaining
time from the nominal 5 ns checkpoint. Actual saved time is zero: all completed.

[configs/posegate_pilot20ns_v2.yaml](configs/posegate_pilot20ns_v2.yaml) uses the
same 5 ns rule with the pilot's (14,20] ns outcome, **not a 20 ns checkpoint**.
Declare a different nominal cadence in a separate configuration before applying
the method to differently sampled files; historical cross-ligand data differ.

## Commands

```bash
posegate inspect --topology system.pdb --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v2.yaml
posegate measure --topology system.pdb --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v2.yaml --output measurements.json
posegate shadow --run-id pocket12_replica1 \
  --topology system.pdb --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v2.yaml --registry shadow_registry/
posegate watch --run-id pocket12_replica1 \
  --topology system.pdb --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v2.yaml --registry shadow_registry/
posegate validate --record shadow_registry/<record-id>.json \
  --topology system.pdb --trajectory trajectory.dcd --output validation.json
posegate report --record shadow_registry/<record-id>.json \
  --validation validation.json --output report.html
```

`measure` can audit completed trajectories. `shadow` refuses a new record once
the late-outcome embargo begins. This CLI implements the measurement and rule,
not the full original prospective enrollment, arming, and capture-deadline
infrastructure. New records are not original prospective observations. The JSON
field `forecast` and enum `FORECASTS_THIS_TRAJECTORY_ONLY` are retained for record
compatibility, not as evidence of advance warning. No action is applied.

## Raw and derived data

Complete raw DCD trajectories for all cohorts are **not distributed because of
their size**. No confidentiality or licensing restriction is asserted here.
Absolute paths in manifests identify sources, not download URLs. Coordinate
digests cover selected atoms, boxes, and indices; they are not full-DCD hashes.

Distributed inputs sufficient for statistical auditing include:

- Primary/expanded comparisons, checkpoints, early-clean subsets, and timing:
  the 81-row table, membership manifest, and per-frame scalar traces.
- Prospective validation/pilot and triage: early/late measurements, eligibility
  flags, frozen threshold, and outcome definitions.
- Separate-launch and cross-ligand analyses: paired or cohort measurement rows.
- Local-frame/protein-side controls: local measurements, C-alpha traces,
  summary rows, and bootstrap outputs; archived radius controls are separate.
- Contact retention: archived predictor/target inputs and result rows.
- Boundary sensitivity: continuous early/late values and original statistical
  implementations; prospective threshold selection remains unchanged.

These data do not reconstruct full trajectories or independently verify every
atom selection. Small coordinate fixtures support regression tests; selected
PDB snapshots support Figure 4. The mapping identifies each audit level.

## Historical material

Obsolete research snapshots, capture archives, and artwork have been removed
from this release tree. They remain in Git history at commit
`df4f6c206c8de38206aa09469db34e21c7eb5eb8`. The current threshold's frozen
policy and development audit are retained under `paper/current/policy/`.
Only two small v1 configuration fixtures remain under `tests/fixtures/legacy/`
to test record compatibility; they are not installed as current policies.

## License

PoseGate-MD is released under the MIT License.
