# Frozen protocol design: same-trajectory 5 ns RMSD screen v2

## 1. Scientific question

For a trajectory that is already running, does its own first 5 ns show enough
protein-relative ligand-pose deviation to mark it as a candidate for early
termination?

This is a one-sided resource-allocation screen. It is not a claim about binding,
biological activity, equilibrium stability, or the outcome of another launch.

## 2. Unit of analysis

One uniquely identified production **launch** contributes:

1. one pre-outcome 5 ns measurement;
2. one frozen action; and
3. one exact 70–100 ns outcome.

A launch ID is required because protein–pocket–replica labels can be reused by a
rerun. Two checkpoint predictions from one trajectory are not independent;
v2 therefore has one primary checkpoint.

## 3. Measurement

The first saved DCD frame is the reference frame. Under the nominal 0.2 ns
reporting cadence it is saved frame 1 at 0.2 ns.

For every required frame:

1. select protein Cα atoms and ligand heavy atoms;
2. determine one periodic shift from the ligand centroid relative to the
   protein-Cα centroid and apply that shift to every ligand atom;
3. reject a raw frame whose ligand diameter exceeds 45% of the shortest box
   edge, indicating that the ligand was written torn across periodic images;
4. center and Kabsch-align protein Cα atoms to the first saved frame;
5. apply the protein-derived rotation to the ligand;
6. compute ligand heavy-atom RMSD against the first saved frame.

The feature is the mean of exactly saved frames 16–25, corresponding to nominal
3.2–5.0 ns and therefore `(3,5]` ns. Integer saved-frame arithmetic is
authoritative. Floating DCD metadata such as `0.200000000596 ns` is checked for
cadence but never used to decide boundary inclusion.

## 4. Decision rule

```text
X5 = mean corrected pose RMSD over exact frames 16–25
stop_candidate if X5 > tau
otherwise do_not_stop
```

`tau` is selected from a strict remeasurement of the locked 52-row development
cohort. The selected threshold is the lowest cutoff that remains within a
predeclared empirical false-stop cap, maximizing nominal compute saving subject
to that cap. With a zero-false-stop cap, this is the retained-envelope boundary:
no retained development trajectory lies above the threshold.

No standardization, sigmoid, classifier, or probability is evaluated at
runtime. The development cap is not a prospective guarantee.

## 5. Prospective capture rules

The full primary cohort and launch order are frozen before production. For each
launch, an immutable arm event must be written while the declared production
DCD is absent. A stale event, capture, validation, or disposition under the same
launch ID blocks arming.

Capture is rejected when:

- no valid pre-production arm event exists;
- fewer than 25 saved frames exist;
- more than 30 saved frames exist, meaning the frozen 6 ns deadline was missed;
- any frame in `(70,100]` already exists;
- code, policy, enrollment, software, topology, path, or file identity fails a
  provenance check;
- two consecutive reads of the exact 25-frame prefix disagree.

The MD job is not paused or killed. Every shadow launch continues to 100 ns.

## 6. Outcome

A completed trajectory must contain saved frame 500. The endpoint is the median
corrected pose RMSD over exactly frames 351–500, nominal 70.2–100.0 ns. The
trajectory is retained only when this median is strictly below 3.0 Å.

The validator reads the exact early and late coordinates twice, permits only
append-only growth, and requires the selected-coordinate hashes, early feature,
late median, and label to agree. It recomputes the early feature and requires
agreement with the capture. A ±0.1 Å boundary flag is descriptive only; it never
changes the strict label.

## 7. Primary metrics

The protocol is asymmetric:

- **false-stop rate among retained trajectories** is the primary safety metric;
- **stop precision** asks how often a proposed stop later fails the endpoint;
- **failure capture** asks how many non-retained trajectories are identified;
- **counterfactual compute saving** assigns 95 nominal ns to each validated
  stop candidate.

A false continue costs compute. A false stop can discard a scientifically useful
retained trajectory. Generic accuracy obscures this difference and is secondary.

Threshold-independent diagnostics are also reported:

- Spearman association between X5 and late median RMSD;
- raw-X5 AUROC for late non-retention.

These diagnostics test whether the observation itself forecasts later behavior,
not whether one selected cutoff is optimal.

## 8. Cohort control and missingness

All primary launches, their order, and minimum launch/protein breadth are frozen
before the first v2 capture. Technical failures remain visible through immutable
disposition records and cannot be replaced after freeze.

A disposition records the DCD frame count and whether partial or complete
outcome data were available. A complete 100 ns endpoint cannot be excluded as
`failed_after_capture_before_outcome`; it must be validated. Dispositions after
partial outcome availability are flagged in the final cohort flow. When a
damaged or unreadable DCD prevents that determination, availability is recorded
as unknown and counted separately rather than treated as absent.

The final primary summary is produced only after every launch is validated or
has a sealed technical disposition. A completed cohort with no validated stop
candidate leaves the key stop branch untested and is labeled accordingly.

## 9. Cohort size and claim level

A practical pilot may freeze 30 launches across at least 10 proteins, but it
remains a pilot. To obtain a one-sided 95% binomial upper bound below 10% after
zero false stops requires 29 retained outcomes because
`1 - 0.05^(1/n) < 0.10`. At roughly 25% retention prevalence this is about 116
total trajectories. This calculation assumes independent retained outcomes;
repeated pockets within a protein can be correlated, so protein breadth and the
per-protein table remain essential.

## 10. Shadow versus action deployment

The shadow stage estimates what would have happened while keeping every outcome
observable. It does not realize compute savings. Enforced stopping requires a
new policy version and a separately frozen action pilot with randomized or
predeclared audit continuation of some stop candidates. Without such an audit,
false stops among terminated runs cannot be measured.

## 11. Versioning

V1 and v2 are different policies because v2 changes frame-boundary handling,
whole-ligand imaging, threshold selection, terminology, provenance, and cohort
control. V1 records remain historical. Any later change to a window, threshold,
action rule, extractor, or endpoint creates a new version and a new prospective
cohort.
