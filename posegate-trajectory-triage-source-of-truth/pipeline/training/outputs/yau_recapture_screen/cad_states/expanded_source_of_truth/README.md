# Expanded same-trajectory source of truth

This directory is the canonical numerical source for the expanded trajectory
analysis. Version: `expanded_same_trajectory_v1`.

## Canonical cohort

- 81 discovered raw DCD trajectories across 16 proteins;
- 73 primary eligible trajectories across all 16 proteins;
- 17 late-retained and 56 late-non-retained outcomes;
- 8 excluded trajectories retained in the master table with explicit reasons.

The canonical row-level table is
`master/expanded_trajectory_source_of_truth.csv`. Use
`eligible_primary_cohort == True` for primary same-trajectory analyses. Do not
construct the cohort by dropping missing values or by reading legacy tables.

## Directory layout

- `full_100ns_corrected_measurements.csv`: direct measurements from raw DCDs;
- `master/expanded_trajectory_source_of_truth.csv`: canonical row-level table;
- `master/cohort_flow.csv`: explicit cohort counts;
- `master/primary_cohort_exclusions.csv`: excluded systems and reasons;
- `master/build_manifest.json`: hashes, versions, and build metadata;
- `primary_analysis/`: expanded corrected-RMSD and LOPO results;
- `physics_guided_ai_audit/`: nested protein-held-out model comparison.

## Full reproduction on the Linux data host

From the repository root, using the remote `openmm-env`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/home/zhenli/miniconda3/envs/openmm-env/bin/python \
yau_competition/scripts/cad_paper_analysis/run_expanded_source_of_truth_pipeline.py \
  --raw-root /media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md \
  --workers 4
```

This rescans every raw trajectory, rebuilds the master dataset, regenerates the
primary analysis, and repeats the nested AI audit. It does not modify raw DCDs.

To rebuild downstream artifacts from an already validated measurement table:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
/home/zhenli/miniconda3/envs/openmm-env/bin/python \
yau_competition/scripts/cad_paper_analysis/run_expanded_source_of_truth_pipeline.py \
  --reuse-measurements
```

The dataset is retrospective and supports a same-trajectory geometric
persistence claim. It is not prospective validation or experimental evidence
of binding.
