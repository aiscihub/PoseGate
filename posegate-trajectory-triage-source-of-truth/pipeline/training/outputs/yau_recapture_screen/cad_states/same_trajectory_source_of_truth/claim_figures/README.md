# Locked claim figures

Generated from `same_trajectory_source_of_truth.csv` by
`26_plot_locked_claim_results.py`.

All primary comparisons use the identical locked cohort: 52 trajectories,
15 proteins, 13 late retained and
39 late non-retained. Seven truncated trajectories
from the legacy outcome are excluded.

## Main files

- `figure_key_locked_claim.png/.pdf`: within-run persistence, coordinate controls,
  same-versus-separate-launch comparison, and retrospective policy frontier.
- `figure_checkpoint_and_oof.png/.pdf`: common-cohort checkpoint curve and outer-LOPO
  probability distribution.
- `key_claim_metrics.csv`: exact AUROCs and protein-clustered intervals.
- `paired_auroc_differences.csv`: paired AUROC differences and clustered intervals.
- `checkpoint_locked_common_cohort.csv`: checkpoint results on the same 52 rows.
- `retrospective_policy_frontier.csv`: every plotted policy threshold.

Intervals resample whole proteins with replacement. They treat each continuous
score or saved OOF probability as fixed and do not refit the complete modeling
procedure within each bootstrap draw. The policy frontier is retrospective and
must not be described as prospective false-stop control. Because the
class-weighted logistic output was not calibrated, its saved `predict_proba`
value is interpreted as a model score rather than an absolute probability.
