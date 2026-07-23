# PoseGate-MD

PoseGate-MD is a local, open-source command-line tool for physically corrected,
same-trajectory geometric forecasting beside a running molecular-dynamics job.
It measures only the coordinate prefix available at a frozen checkpoint,
records an immutable shadow forecast, and later checks that forecast against a
predefined late trajectory window.

Version 0.1 is deliberately a measurement and shadow-recording product. It
does not terminate simulations, claim binding affinity, predict a separately
launched replica, or report a per-trajectory confidence interval.

## Installation

From a source checkout:

```bash
python -m pip install -e ".[test]"
posegate --version
```

## Scientific scope

The validated primary observable is protein-relative ligand heavy-atom RMSD.
PoseGate-MD:

1. requires periodic box vectors;
2. represents ligand atoms by minimum-image vectors from the protein C-alpha
   centroid;
3. aligns protein C-alpha coordinates to frame 0;
4. applies the protein-derived rotation to the ligand;
5. measures corrected ligand RMSD and centroid displacement; and
6. optionally reports frame-relative reorientation and internal deformation.

The bundled policies reproduce the frozen 5 ns and 20 ns shadow specifications
used by the accompanying manuscript. Their scores are uncalibrated ranking
scores, not probabilities of binding or drug activity.

## Initial commands

```bash
# Verify inputs, selections, timing, and checkpoint coverage.
posegate inspect \
  --topology system.pdb \
  --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v1.yaml

# Measure the configured prefix without creating a forecast.
posegate measure \
  --topology system.pdb \
  --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v1.yaml \
  --output measurements.json

# Create one immutable shadow record before late frames exist.
posegate shadow \
  --run-id pocket12_replica1 \
  --topology system.pdb \
  --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v1.yaml \
  --registry shadow_registry/

# Watch a growing trajectory until the configured checkpoint is ready.
posegate watch \
  --run-id pocket12_replica1 \
  --topology system.pdb \
  --trajectory trajectory.dcd \
  --config configs/posegate_5ns_v1.yaml \
  --registry shadow_registry/

# Compare the sealed forecast with a completed trajectory.
posegate validate \
  --record shadow_registry/<record-id>.json \
  --topology system.pdb \
  --trajectory trajectory.dcd
```

## Output language

Forecasts use `POSE_RETAINED`, `POSE_NONRETAINED`, or `DEFER`. Reversible
scheduling recommendations use `CONTINUE`, `CANDIDATE_FOR_PAUSE`, or
`NO_ACTION_SHADOW_MODE`.

Every record states:

> This forecast applies only to continuation of the trajectory already
> observed. It has not been validated for predicting a separately launched
> simulation.

## Development boundary

The numbered manuscript-analysis scripts remain the research provenance
layer. PoseGate-MD extracts their reusable coordinate operations and must pass
synthetic tests plus golden parity checks before release. OpenMM directory
adapters, HTML reports, and reversible GPU scheduling follow the core
measurement and shadow workflow. Automatic permanent termination is outside
version 0.1.

## Reproducibility archive

The repository also contains
`posegate-trajectory-triage-source-of-truth/`, the preserved scripts, tables,
and manuscript snapshot underlying the initial scientific analysis. That
archive remains separate from the reusable package so historical results can
be reproduced without silently changing the original numbered scripts.

## License

PoseGate-MD is released under the MIT License.
