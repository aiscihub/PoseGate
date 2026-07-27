# Code fixes implemented in v2

| V1 problem | V2 correction |
|---|---|
| Logistic score obscured a one-feature threshold | Runtime classifier removed; direct RMSD cutoff only |
| `continue` sounded like predicted retention | One-sided actions: `stop_candidate` and `do_not_stop` |
| 5 ns and 20 ns calls on one run inflated apparent n | One primary 5 ns decision per launch |
| Floating DCD `dt` admitted the nominal 3.0 ns frame | Integer indices select exactly frames 16–25 |
| Floating boundary admitted the nominal 70.0 ns frame | Integer indices select exactly frames 351–500 |
| Outcome coverage accepted runs short of 100 ns | Saved frame 500 is mandatory |
| Atom-wise PBC imaging could split a ligand | One centroid-derived periodic shift is applied to every ligand atom |
| A differently written DCD could contain an already torn ligand | Raw ligand diameter guard rejects full-box atom-wise tears |
| Script hash was recorded but not enforced | Policy binds and verifies every runtime implementation file |
| Frozen model depended on sklearn/joblib versions | No runtime sklearn object or unpickling |
| Reruns could reuse the same complex label | Required unique launch ID and frozen enrollment |
| New cases could be selected opportunistically | Entire primary cohort and launch order are frozen before production |
| “Prospective” monitoring was not proven to start before the DCD | Immutable arm event is created only while the production DCD is absent |
| A stale production event could be reused | Arming rejects any pre-existing launch event, capture, validation, or disposition |
| The topology can exist before the DCD, causing the pre-armed monitor to exit | A one-file-only startup state is treated as not ready; the monitor keeps waiting |
| Tiny or under-diverse cohorts could be frozen implicitly | Enrollment declares minimum launch and protein counts |
| Failed launches could disappear or be replaced | Immutable technical dispositions remain in cohort flow |
| Technical exclusions could be recorded after outcome data existed without disclosure | Disposition records snapshot frame availability and flag partial/complete outcome availability |
| An unreadable DCD could make outcome availability look falsely absent | Unknown availability is preserved and counted as a separate missingness risk |
| Capture could occur long after 5 ns | Frozen 6 ns deadline plus outcome embargo |
| Live DCD could change during read | Prefix is read twice and coordinate hashes must agree |
| Final DCD could still be flushing during validation | Exact early/outcome coordinates are read twice and must agree |
| A crash could expose a partially written final JSON | Records are written to a temporary inode and atomically linked/replaced into place |
| Concurrent index rebuilds could lose a capture row | Global file lock plus atomic index replacement |
| Hash integrity alone did not verify policy semantics | Validator and summarizer recompute frame, deadline, embargo, threshold, and margin invariants |
| Validator claimed stop calls were uncheckable | Every shadow run continues, so both branches are validated |
| Overall accuracy treated both errors alike | False-stop safety and stop utility are primary |
| Near-3 Å outcomes invited post-hoc relabeling | Strict label is unchanged; signed margin and 0.1 Å boundary flag are separate |
| The observation itself was hidden by threshold metrics | Summary also reports descriptive Spearman association and raw-X5 AUROC |
| A completed cohort with no stop calls could look “validated” | Summary explicitly flags `completed_flow_but_stop_branch_unexercised` |
| Existing v1 threshold silently assumed valid | Strict 52-row remeasurement precedes direct threshold selection |
