# Shadow-test v2 verification report

This report documents checks that can be completed from the uploaded bundle.
It does not claim that the development remeasurement or a prospective cohort
has been run.

## Scope

The verification covers the v2 implementation, synthetic trajectory behavior,
artifact immutability, provenance chaining, direct threshold selection, and a
synthetic arm-to-capture-to-outcome workflow.

## Automated tests

The test suite exercises:

- exact integer frame selection despite noisy DCD time metadata;
- removal of global rigid motion and detection of relative ligand rotation;
- whole-ligand periodic imaging and rejection of a torn ligand;
- one-sided action/category semantics;
- rejection of non-finite or negative timing values;
- semantic capture-contract verification beyond file hashes;
- crash-safe, write-once JSON publication;
- false-stop-constrained direct threshold selection;
- policy binding to the strict audit and selected threshold;
- robust CSV boolean parsing;
- pre-production arming and stale-event rejection;
- the normal topology-before-DCD startup state;
- write-once arm and production-start events;
- technical-disposition timing and cohort-flow preservation;
- exact 10-frame early and 150-frame outcome extraction from a synthetic DCD;
- refusal to dispose a complete endpoint as a pre-outcome failure;
- threshold-independent Spearman and AUROC diagnostics; and
- synthetic end-to-end arm, capture, append-only growth, validation, and
  launch-level summary.

Commands executed in this review environment:

```bash
python -m py_compile \
  yau_competition/scripts/cad_paper_analysis/shadow_v2_core.py \
  yau_competition/scripts/cad_paper_analysis/52_shadow_monitor_v2.py \
  yau_competition/scripts/cad_paper_analysis/53_validate_shadow_outcome_v2.py \
  yau_competition/scripts/cad_paper_analysis/54_shadow_reimplementation_audit_v2.py \
  yau_competition/scripts/cad_paper_analysis/55_select_shadow_threshold_v2.py \
  yau_competition/scripts/cad_paper_analysis/56_summarize_shadow_v2.py

sh -n yau_competition/scripts/cad_paper_analysis/run_shadow_monitor_v2.sh
pytest -q yau_competition/tests/test_shadow_v2.py
```

Result: all Python files compiled, shell syntax passed, and **21 tests passed**.
All five command-line entry points loaded with `--help`. The historical
comparison command re-derived the legacy boundary as
`2.518358545938333 Å`.

## Checks intentionally not performed here

The raw 52 development trajectories are not present in the uploaded bundle.
Therefore the following deployment-gate steps cannot be honestly completed in
this environment:

1. strict remeasurement of all 52 development trajectories;
2. inspection of old-versus-strict X5 and late-outcome differences;
3. selection of the final v2 threshold under the declared false-stop cap;
4. freezing a deployable v2 policy; or
5. running a new prospective cohort.

No threshold or policy artifact is fabricated. The legacy 2.518 Å-equivalent
boundary remains comparison-only.

## Acceptance criteria before live use

Live use is blocked unless all of these pass in the actual project repository:

- all 52 trajectories remeasure successfully;
- every early window has 10 frames and every outcome window has 150 frames;
- any changed late labels are reviewed and retained as v2 differences rather
  than silently reconciled;
- the false-stop cap is declared before threshold selection;
- all runtime files are committed and clean at policy freeze;
- the complete enrollment is frozen before any enrolled production DCD exists;
- each monitor is armed and its immutable arm record confirmed before MD starts;
- every launch is later validated or receives a visible technical disposition;
- v1 and v2 results remain separate.
