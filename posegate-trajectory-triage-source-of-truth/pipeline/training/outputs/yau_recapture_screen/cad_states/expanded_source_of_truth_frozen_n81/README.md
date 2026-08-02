# Expanded cohort, frozen n=81

Canonical version: `expanded_same_trajectory_frozen_n81`.

This is the exact 81-trajectory cohort behind the published expanded-cohort
numbers in `yau_aj.tex`: Table~\ref{tab:checkpoint} ("Expanded $n=81$
AUROC"), Table~\ref{tab:post-lock-additions} ("Combined expanded" row), and
Figure 2 Panel A/B. It supersedes `expanded_source_of_truth_v2_incl_reruns`
(n=79) and `expanded_source_of_truth` (n=76) for any analysis that needs to
match currently-published paper numbers.

## How this cohort was assembled

Starting point: `expanded_source_of_truth_v2_incl_reruns` (n=79 eligible,
84 discovered). Two changes were applied:

**1. Four prospective shadow-test trajectories, freshly completed, folded in
as ordinary primary-cohort rows once their ground truth became available.**
Their raw DCD trajectories in the original `simulation_100ns_md` tree were
truncated (writer/log mismatch — the OpenMM simulation log reported
completion to ~100.1 ns but the DCD reporter had only flushed a fraction of
that): MDR1_TRIRC|pocket4 (402 frames / 80.4 ns) and PDH1_CANGA|pocket13
(166 frames / 33.2 ns) were missing their outcome window entirely;
CDR2_CANAL|pocket8 (474 frames / 94.8 ns) and CIMG_06197|pocket15 (199
frames / 39.8 ns) were close but incomplete. Complete 500-frame / 100 ns
reruns exist for all four under `simulation_100ns_md_reruns/`, which is
where the shadow-test monitoring pipeline (`yau_competition/prospective_shadow/`)
had already been reading from — its `outcome_validation_*.json` files for
these four cases report the identical `same_traj_pose_rmsd_18_20_A` values
recovered here (cross-check delta 0.0 Å in every case). Measurements were
regenerated directly from those reruns using the unmodified
`28_generate_full_100ns_measurements.py::measure()` routine (no new logic
written), then folded into `31_build_expanded_source_of_truth.py`:

| complex_id | rerun trajectory source | n_frames | outcome |
|---|---|---|---|
| MDR1_TRIRC\|Milbemycin\|pocket4\|replica_1 | `simulation_100ns_md_reruns/MDR1_TRIRC/.../pocket4/replica_1/` | 500 (100.0 ns) | non-retained |
| PDH1_CANGA\|Milbemycin\|pocket13\|replica_1 | `simulation_100ns_md_reruns/PDH1_CANGA/.../pocket13/replica_1/` | 500 (100.0 ns) | non-retained |
| CDR2_CANAL\|Milbemycin\|pocket8\|replica_1 | `simulation_100ns_md_reruns/CDR2_CANAL/.../pocket8/replica_1/` | 500 (100.0 ns) | non-retained |
| CIMG_06197\|Milbemycin\|pocket15\|replica_1 | `simulation_100ns_md_reruns/CIMG_06197/.../pocket15/replica_1/` | 500 (100.0 ns) | non-retained |

All four were prospectively predicted "stop" by the shadow monitor at both
the 5 ns and 20 ns checkpoints, and all four turned out non-retained — the
four correct "stop" calls in the shadow-test record.

**2. Two negative-control pockets explicitly excluded from the headline
cohort**, despite having complete, eligible 100 ns data (see
`master/negative_controls_excluded.csv`):

| complex_id | reason |
|---|---|
| CDR2_CANAL\|Milbemycin\|pocket14\|replica_1 | negative-control pocket, not part of the locked headline cohort |
| MDR1_CRYNH\|Milbemycin\|pocket12\|replica_1 | negative-control pocket, not part of the locked headline cohort |

Two other trajectories that were also flagged as negative-control candidates
— MDR1_CRYNH|pocket1 and PDR5_YEAST|pocket14 — were checked and are **not**
both excluded: MDR1_CRYNH|pocket1 is confirmed present in the locked n=81 (it
was already complete in the original, non-rerun tree, unrelated to shadow
monitoring); PDR5_YEAST|pocket14 is confirmed absent (its data was still
incomplete in the snapshot the freeze used, and it is not part of n=81 — its
shadow-test "continue" call, later shown wrong, is not yet reflected in this
retrospective cohort).

Net effect: 79 − 2 (excluded) + 4 (shadow completions) = **81**.

## Independent verification

`doc/YAU_competition_2026/figures/figure2_frozen_expanded_3_5ns_raw_data.csv`
is the exact source table behind Figure 2 and the published 5 ns AUROC
values (0.856 frozen-52 / 0.821 expanded-81). Its 81-row
`expanded_audit_n81` cohort's `complex_id` set is **identical** to this
frozen cohort's eligible set. For the 77 rows this cohort shares with the
pre-existing v2_incl_reruns snapshot (i.e. excluding the 4 freshly-completed
shadow rows above), `early_corrected_rmsd_mean_3_5_A` and `late_pose_retained`
match that file exactly — 0.0 Å max delta, 0 label mismatches — confirming
the underlying measurement pipeline is consistent with whatever produced the
paper's canonical figure data.

Running `59_checkpoint_curve_and_robustness.py` on this locked master
reproduces **all five** published Table 3 expanded-cohort checkpoints,
point estimate and 95% protein-clustered CI, to the paper's own rounding:

| Checkpoint | Computed AUROC (95% CI) | Published Table 3 |
|---|---|---|
| 5 ns  | 0.8208 [0.7477, 0.8842] | 0.821 [0.748, 0.884] |
| 10 ns | 0.8015 [0.6902, 0.8836] | 0.801 [0.690, 0.884] |
| 15 ns | 0.8143 [0.7070, 0.8867] | 0.814 [0.707, 0.887] |
| 20 ns | 0.8180 [0.7312, 0.8884] | 0.818 [0.731, 0.888] |
| 30 ns | 0.8612 [0.7682, 0.9316] | 0.861 [0.768, 0.932] |

See `build_manifest.json` for the exact machine-readable side-by-side and
file hashes.

## Files

- `master/expanded_trajectory_source_of_truth.csv`: full 84-row discovered
  table; use `eligible_primary_cohort == True` for the locked 81-row
  headline cohort (81 rows satisfy this).
- `master/negative_controls_excluded.csv`: the 2 rows deliberately excluded
  from the primary cohort despite complete data, and why.
- `full_100ns_corrected_measurements.csv`: the raw-trajectory measurement
  table this master was built from (84 rows; 4 replaced with fresh
  reruns-tree measurements relative to v2_incl_reruns).
- `checkpoint_and_robustness/checkpoint_curve.csv`,
  `robustness_subsets.csv`: outputs of `59_checkpoint_curve_and_robustness.py`
  on the locked master, reproducing Table 3's expanded-cohort column.
- `build_manifest.json`: row counts, verification numbers, and SHA-256
  hashes of the master/measurement tables.

## Provenance / how to regenerate

1. Base measurements: `expanded_source_of_truth_v2_incl_reruns/full_100ns_corrected_measurements.csv`.
2. Replace the 4 rows listed above with fresh `28_generate_full_100ns_measurements.py::measure()`
   output computed against `simulation_100ns_md_reruns/<protein>/simulation_explicit/<pocket>/replica_1/`.
3. Run `31_build_expanded_source_of_truth.py` on the merged 84-row measurement table.
4. Set `eligible_primary_cohort = False`, `late_pose_retained = <NA>`,
   `primary_exclusion_reason = "excluded_negative_control_pocket_not_in_locked_headline_cohort"`
   for `CDR2_CANAL|Milbemycin|pocket14|replica_1` and
   `MDR1_CRYNH|Milbemycin|pocket12|replica_1`.
5. Run `59_checkpoint_curve_and_robustness.py` on the result.

This dataset is retrospective. It is the row-level numerical source of
truth for the paper's currently-published expanded-cohort headline numbers,
not evidence of prospective validation or experimental binding.
