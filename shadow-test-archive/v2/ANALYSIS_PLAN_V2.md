# Prospective analysis plan v2

## Frozen question

Can a direct 5 ns corrected-RMSD threshold identify a subset of ongoing
trajectories that later fail the exact 70–100 ns pose-retention endpoint while
keeping false stops acceptably rare?

This is a one-sided early-stop question. `do_not_stop` is not a forecast of
retention.

## Freeze order

Before the first v2 production DCD exists:

1. remeasure the locked 52-row development cohort with the strict v2 extractor;
2. declare the maximum empirical development false-stop count;
3. select and freeze one RMSD threshold;
4. commit and hash the complete runtime implementation;
5. list every prospective launch, its order, primary status, minimum total
   launch count, and minimum distinct-protein count;
6. freeze enrollment; and
7. arm each monitor before its enrolled production DCD exists.

No threshold, window, endpoint, action rule, enrollment row, or primary metric
changes after the first capture. A change creates a new policy version and a new
prospective cohort.

## Analysis unit and endpoint

The unit is one enrolled launch. It contributes at most one primary 5 ns action
and one strict late outcome.

- Feature: mean corrected ligand-pose RMSD over exact zero-based saved-frame
  indices `15..24` (nominal 3.2–5.0 ns).
- Action: `stop_candidate` only when `X5 > tau`; equality is `do_not_stop`.
- Outcome: median corrected RMSD over exact indices `350..499` (nominal
  70.2–100.0 ns); retained only when the median is strictly below 3.0 Å.

The 3.0 Å endpoint is never changed for a borderline case. Signed distance to
3.0 Å and a descriptive ±0.1 Å boundary flag are reported separately.

## Primary outcome categories

| Action | Late outcome | Category |
|---|---|---|
| `stop_candidate` | non-retained | `successful_stop` |
| `stop_candidate` | retained | `false_stop` |
| `do_not_stop` | retained | `safe_continue` |
| `do_not_stop` | non-retained | `false_continue` |

A false stop is the primary scientific safety error. A false continue is a
compute-only error.

## Primary operating estimands

All are calculated only on launches with a sealed capture and complete outcome:

1. false-stop count and rate among retained trajectories;
2. precision of stop candidates;
3. capture rate of non-retained trajectories; and
4. counterfactual nominal compute saving, assigning 95 ns to each validated
   stop candidate.

Raw numerators and denominators are reported first. Launch-level Wilson
intervals are descriptive because launches within a protein may be correlated.
The per-protein outcome table is always reported beside them. Overall binary
accuracy is secondary.

## Threshold-independent scientific diagnostics

To test the core observation rather than only one operating point, the summary
also reports on completed launches:

- pooled Spearman rank correlation between X5 and late median RMSD;
- pooled raw-X5 AUROC for late non-retention.

These are descriptive launch-level diagnostics. A formal report may add a
predeclared protein-clustered bootstrap, but no post-hoc choice between pooled
and clustered results may alter the frozen action rule.

## Cohort flow and missingness

Every enrolled primary launch appears in the final flow as exactly one of:

- validated outcome;
- immutable technical disposition; or
- unresolved at the reporting date.

A technical disposition does not erase primary membership and does not permit a
replacement after freeze. It is excluded from performance denominators because
no evaluable late outcome exists, but its code, reason, DCD frame count, and
whether partial or complete outcome data were available remain visible.

A trajectory with a complete 100 ns endpoint cannot be classified as
`failed_after_capture_before_outcome`; it must be validated. Dispositions after
partial outcome-window availability are counted explicitly as a missingness
risk. When a damaged or unreadable DCD prevents the availability check, that
state is counted as unknown rather than treated as evidence that no outcome data
existed. Final claims should include a conservative sensitivity discussion
whenever either kind of disposition occurs.

## Timing and peeking

Interim summaries may be generated for engineering checks, but they do not
change the rule, sample, endpoint, or claim threshold. The final primary summary
is produced only after every planned launch has either a validated outcome or a
sealed technical disposition. V1 and v2 outcomes are never pooled.

## Stop-branch requirement

The main operational benefit is the stop branch. Therefore:

- zero validated `stop_candidate` calls means stop precision and prospective
  false-stop behavior at the operating point are unestimated;
- a completed cohort in that state is reported as
  `completed_flow_but_stop_branch_unexercised`, not as validation of the
  stopping rule;
- the threshold is not lowered after seeing this result. A different threshold
  requires a new policy and cohort.

## Claim levels

- **Engineering demonstration:** at least one valid capture and later outcome;
  establishes that mechanics work, not predictive validity.
- **Prospective pilot:** the pre-enrolled cohort is complete, the stop branch is
  exercised, all primary estimands and continuous diagnostics are reported, and
  protein breadth/cohort flow are shown. No guaranteed error rate is claimed.
- **Risk-control claim:** requires a separately justified target, sufficient
  retained outcomes across broad protein groups, and a predeclared pass
  criterion. A simple zero-error binomial calculation is insufficient when
  outcomes are concentrated within a few proteins.
- **Realized compute-saving claim:** requires a later action pilot. The shadow
  stage alone supports only counterfactual savings.
