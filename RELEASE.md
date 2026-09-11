# Manuscript release v0.2.0

This release freezes the manuscript evidence bundle, whole-ligand measurement, and
one-sided Table 5 rule. Its annotated Git tag points to a content-addressed
commit. Obtain the immutable full ID with `git rev-parse 'v0.2.0^{commit}'`;
the tag annotation also records it. Do not force-move the tag.

`SHA256SUMS` covers release files, frozen inputs, results, configurations, and
fixtures. It excludes itself; the Git commit binds the checksum list.
`paper/current/manifest.json` also records original paths and file hashes.
Neither scheme claims to hash undistributed full DCD files.

The full manuscript and supplement are unsubmitted and must not be published
with this release. Only the requested table extracts, numerical data, figure
assets, and provenance mapping are included; full-document hashes contain no text.

## Scientific contract

- Whole-ligand periodic shift and whole-protein C-alpha reference alignment.
- First saved production frame, fixed correspondence, exact nominal windows.
- X5 over (3,5] ns; strict one-sided threshold 2.7634173197224974 angstrom.
- Late non-retention means L100 >= 3 angstrom, not irreversible departure.
- Monitoring and shadow triage; no applied stop or demonstrated advance warning.
- Legacy configuration bytes remain in compatibility tests; obsolete models and
  captures are recoverable from Git history, not bundled with this paper release.

Tested numerical environment: Python 3.11, NumPy 2.4.6, MDAnalysis 2.10.0,
pandas 2.3.3, SciPy 1.17.1, PyYAML 6.0.3, pytest 8.4.2. Original measurement
metadata preserve their producing versions. Saved-frame regression tolerance
is 2e-5 in the reported coordinate/angle unit, not a policy threshold tolerance.

## Audit boundary

The portable verifier recomputes point estimates, checks tables against frozen
results, and reruns selected clustered intervals. It does not rerun every
bootstrap or raw-coordinate analysis. Whole proteins are resampled; original
RNG choices remain distinct. The primary corrected interval uses Python Random;
current coordinate controls use NumPy's generator.

The exact 20 ns local corrected AUROC is 0.8125, displayed as 0.813 by half-up
rounding. The original validation summary preserves accounting from 6 ns
(846 ns); Table 5 counts from 5 ns (855 ns). Both are counterfactual, with zero
actual saved time. The verifier uses the current table's convention.

Byte-preserved `paper/current/source_analysis/` scripts are provenance, not
repackaged CLI entry points. Their original main functions may require research
paths and external raw files. Use `tools/verify_paper_results.py` for portable
auditing. Authored diagrams and protocol tables have no statistical generator.

## Maintainer procedure

1. Export with `python tools/sync_paper_results.py --source-root /path/to/valleyfevermutation`.
2. Run `python tools/export_manuscript_tables.py`, then the table verifier and tests; review the mapping and hashes.
3. Run `python tools/release_checksums.py --write`, then verify without `--write`.
4. Commit reviewed files and create an annotated `v0.2.0` tag recording the commit.
5. Publish the commit/tag only with authorization. A local tag is not a hosted release.

Do not overwrite this snapshot for later analyses; create a new release.
