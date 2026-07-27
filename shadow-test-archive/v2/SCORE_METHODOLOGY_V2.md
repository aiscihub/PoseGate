# How the shadow-test v2 action is calculated

V2 does not produce a retention probability. It computes one geometric
measurement in ångströms and compares it with one frozen cutoff.

## 1. Exact early observation

For an enrolled launch, the monitor reads only the first 25 saved DCD frames.
With the frozen 0.2 ns reporting cadence, saved frame 0 is nominally 0.2 ns and
saved frame 24 is nominally 5.0 ns.

The feature window `(3,5]` is selected by integer index, never by floating DCD
time metadata:

- zero-based indices `15..24`;
- nominal times `3.2, 3.4, ..., 5.0` ns;
- exactly 10 frames.

For every selected frame, the implementation:

1. selects protein Cα atoms and ligand heavy atoms;
2. checks that the raw ligand has not been torn across periodic images;
3. determines one periodic shift from the ligand centroid relative to the
   protein Cα centroid and applies that same shift to all ligand atoms;
4. centers and Kabsch-aligns the current protein Cα coordinates to saved frame
   0; and
5. applies the protein-derived rotation to the ligand before measuring ligand
   displacement from its saved-frame-0 pose.

For ligand heavy atoms `a = 1,...,N`, the per-frame observable is

```text
RMSD(t) = sqrt((1/N) * sum_a ||r_aligned,a(t) - r_reference,a||^2)
```

The live feature is

```text
X5 = mean of RMSD(t) over saved-frame indices 15..24.
```

The coordinate prefix is read twice. The capture proceeds only when the exact
coordinate hash and X5 agree across both reads. The DCD may append new frames,
but already written prefix coordinates may not change and the file inode may
not be replaced.

## 2. Development-only threshold selection

Before prospective enrollment, all 52 locked development trajectories must be
remeasured with the identical strict extractor:

- 10 early frames at indices `15..24`;
- 150 late frames at indices `350..499`;
- whole-ligand periodic imaging; and
- the unchanged late endpoint, median corrected RMSD `< 3.0 Å`.

A maximum empirical development false-stop count is declared before running the
selection script. The selector enumerates the observed X5 values and chooses
the lowest threshold that does not exceed that cap. Because a stop is proposed
when `X5 > threshold`, this is the most compute-saving development operating
point that satisfies the declared cap.

No logistic regression, standardization, sigmoid, calibration, v1 shadow
outcome, or v2 shadow outcome is used in selection.

The development cap is not a prospective guarantee. It merely defines the
operating point that the separately frozen v2 cohort will test.

## 3. Live action rule

For the frozen threshold `τ`:

```text
signed_stop_margin_A = X5 - τ

stop_candidate  if X5 > τ
do_not_stop     if X5 <= τ
```

Equality is deliberately conservative: it is `do_not_stop`.

`do_not_stop` means only that the one-sided early-stop condition was not met.
It is not a forecast or probability of late retention. A later non-retained
outcome after `do_not_stop` is a **false continue**, which wastes compute but
does not prematurely discard a retained trajectory.

Every v2 capture records the raw X5 value, the frozen threshold, the signed
margin, exact frame indices, coordinate hash, policy and enrollment hashes,
pre-production arm event, production-start event, topology hash, DCD
path/device/inode, runtime implementation hashes, and software versions.

## 4. Exact late outcome

Every shadow launch continues to at least 500 saved frames. The late endpoint
uses integer indices:

- zero-based indices `350..499`;
- nominal times `70.2, 70.4, ..., 100.0` ns;
- exactly 150 frames.

```text
late_median = median corrected pose RMSD over indices 350..499
late_retained = late_median < 3.0 Å
```

The strict comparison is never changed for a boundary case. The validator also
reports signed distance from 3.0 Å and whether the outcome is within 0.1 Å, but
those fields are descriptive and do not relabel the outcome.

The exact early and late coordinates are read twice during validation. The
validator requires the same append-only DCD identity, complete 100 ns coverage,
stable selected-coordinate hashes, and agreement between the original captured
X5 and its later recomputation.

## 5. Outcome categories

| Captured action | Late endpoint | Category | Meaning |
|---|---|---|---|
| `stop_candidate` | non-retained | `successful_stop` | compute could have been saved |
| `stop_candidate` | retained | `false_stop` | a useful retained run would have been lost |
| `do_not_stop` | retained | `safe_continue` | retained run preserved |
| `do_not_stop` | non-retained | `false_continue` | remaining compute would be spent |

False-stop count/rate, stop precision, non-retained capture, and
counterfactual nominal saving are primary. Overall binary accuracy is
secondary. Spearman association between X5 and late median RMSD and raw-X5
AUROC for late non-retention are reported separately to evaluate the core
observation without depending on the selected threshold.

## 6. Relationship to legacy v1

The old 5 ns logistic model's score-0.5 boundary is mathematically equivalent
to approximately `2.518358545938333 Å` under the old extractor. V2 can derive
that number for historical comparison, but it does not deploy it automatically.
The old extractor may have admitted 11 early frames and 151 outcome frames due
to floating-boundary handling. A v2 cutoff may be frozen only after the strict
52-trajectory remeasurement and direct threshold selection are complete.

V1 captures remain immutable historical records and are never converted into
v2 observations.
