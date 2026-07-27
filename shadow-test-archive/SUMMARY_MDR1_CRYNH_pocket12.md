# Legacy v1 shadow summary — MDR1_CRYNH:pocket12

Status: **complete**. Both legacy checkpoint records were captured before any
`(70,100]` ns data existed; the trajectory later finished to 100 ns and its
outcome was validated.

## Predictions versus outcome

| Checkpoint | Feature(s) | Score | Legacy action | Captured |
|---|---|---:|---|---|
| 5 ns | corrected RMSD `(3,5]` = 1.569 Å | 0.717 | continue | 2026-07-22T19:10:42Z |
| 20 ns | RMSD `(18,20]` = 2.052 Å; centroid displacement = 1.847 Å | 0.720 | continue | 2026-07-23T16:05:33Z |

**Strict outcome:** median corrected RMSD over `(70,100]` ns = **3.047 Å**,
therefore non-retained under the frozen `< 3.0 Å` definition. Both checkpoint
records are false continues under the strict rule.

## Boundary context without relabeling

The outcome lies 0.047 Å above the threshold. That makes it boundary-adjacent,
but it does not make the strict result correct. The manuscript's retrospective
threshold-sensitivity analysis shows that nearby late-outcome cutoffs preserve
the qualitative ranking signal; it does **not** provide a measurement-noise
estimate and does not authorize changing this prospective label after the fact.

The correct reporting is therefore:

- strict category: false continue;
- raw late median: 3.047 Å;
- signed distance from the 3.0 Å endpoint: -0.047 Å on the retained-margin
  convention;
- boundary context: descriptive only.

## Combined legacy v1 record to date

| Launch | Checkpoint records | Legacy action(s) | Actual outcome | Launch-level category |
|---|---:|---|---|---|
| MDR1_CRYNH:pocket1 | 1 (20 ns) | continue | 3.575 Å, non-retained | false continue |
| MDR1_CRYNH:pocket12 | 2 (5 and 20 ns) | continue, continue | 3.047 Å, non-retained | false continue |

There are **two completed launches**, and both are false continues. There are
three checkpoint records only because one launch was scored twice. Those three
records are not independent prospective cases and should not be summarized as
an independent “0/3” validation result. No legacy stop call has yet been tested,
so the proposed compute-saving branch remains unexercised.

Source artifacts remain unmodified and hash-audited:
`captures/MDR1_CRYNH__Milbemycin__pocket12__replica_1__same_trajectory_checkpoint_{5,20}ns_v1.json`
and
`outcome_validation_MDR1_CRYNH__Milbemycin__pocket12__replica_1.json`.
