# PoseGate Core Development Guide

## Purpose

This document explains the core technical workflow behind **PoseGate**, from the OpenMM molecular-dynamics simulation to the final same-trajectory ligand-pose measurement.

It is written for developers, reviewers, and future contributors who need to understand:

- how the molecular-dynamics trajectory is produced;
- what information is stored in the PDB and DCD files;
- how MDAnalysis reads and processes the trajectory;
- why coordinate correction is necessary;
- how the corrected ligand RMSD is calculated;
- what the measurement means physically;
- what the measurement does **not** mean;
- which implementation details still require careful testing.

The central observable discussed here is the value currently represented by names such as:

```text
same_traj_pose_rmsd_18_20_A
```

or, in the checkpoint analysis script:

```text
rmsd_corrected_20ns
```

This quantity is a **physical distance in ångströms**, not a probability and not inherently a machine-learning score.

---

# 1. End-to-End System Overview

The complete workflow has four major stages:

```text
Prepared protein–ligand system
        ↓
OpenMM molecular-dynamics simulation
        ↓
PDB topology + DCD trajectory
        ↓
MDAnalysis trajectory processing
        ↓
Protein-relative, PBC-corrected ligand measurements
        ↓
Checkpoint feature calculation
        ↓
Same-trajectory stop/continue evaluation
```

The responsibilities are different:

| Component | Main responsibility |
|---|---|
| OpenMM | Simulates the motion of atoms over time |
| PDB file | Describes atom identities, residues, and ordering |
| DCD file | Stores a sequence of atomic-coordinate snapshots |
| MDAnalysis | Reads the files and calculates measurements |
| PoseGate analysis | Converts trajectory measurements into same-run forecasting features |

MDAnalysis does **not** perform the molecular-dynamics simulation. It analyzes coordinates that OpenMM has already generated.

---

# 2. OpenMM Molecular-Dynamics Simulation

## 2.1 What OpenMM is doing

A molecular-dynamics simulation treats atoms as particles moving under forces.

At a high level, every integration step performs a process like:

```text
Current atomic coordinates
        ↓
Calculate forces from the force field
        ↓
Update velocities
        ↓
Update coordinates
        ↓
Advance simulation time
```

The forces include contributions from:

- chemical bonds;
- bond angles;
- torsional rotations;
- electrostatic interactions;
- van der Waals interactions;
- water molecules;
- dissolved ions;
- pressure and temperature control.

The simulation does not directly calculate whether a ligand is “bound” or “stable.” It only propagates atomic positions and velocities according to the selected physical model.

---

## 2.2 System contents

The current systems contain:

- an ABC-transporter protein;
- a milbemycin ligand in a predicted pocket;
- explicit water;
- sodium and chloride ions;
- a periodic simulation box.

The present scientific scope is therefore:

> Ligand-pose behavior in explicit-solvent simulations of predicted ABC-transporter pockets.

It is not a direct model of:

- experimental binding affinity;
- transporter inhibition;
- antifungal activity;
- resistance reversal;
- full membrane-embedded transporter behavior.

---

## 2.3 Preparation stages

Before production MD, the system goes through several preparation stages.

### A. Protein–ligand starting structure

The starting complex contains:

- a prepared protein structure;
- a selected pocket;
- a docked ligand pose.

The docked pose is a hypothesis. It is not yet proof that the ligand will remain in that position.

### B. Force-field parameter assignment

The protein and ligand are assigned molecular-mechanics parameters.

These parameters describe:

- equilibrium bond lengths;
- equilibrium angles;
- torsional preferences;
- atomic partial charges;
- van der Waals properties.

### C. Solvation and ions

Water molecules and ions are added around the complex.

The water box is periodic, meaning the simulation box repeats in every direction. An atom or molecule leaving one side of the box can appear through the opposite side.

### D. Energy minimization

Energy minimization removes severe steric clashes and unusually high-energy contacts.

Minimization is not a time-resolved molecular trajectory. It is a geometry-relaxation step.

### E. Equilibration

Equilibration allows the system to adjust before the production run.

During equilibration:

- water reorganizes;
- ions move;
- temperature stabilizes;
- pressure and box dimensions adjust;
- local contacts relax.

### F. Production MD

Production MD creates the trajectory used for the scientific analysis.

The current nominal protocol uses:

- 100 ns production time;
- 4 fs integration time step;
- 25,000,000 production steps;
- coordinates saved every 50,000 steps;
- one saved frame every 0.2 ns;
- approximately 500 frames for a complete 100 ns trajectory.

Calculation:

```text
50,000 steps × 4 fs/step
= 200,000 fs
= 200 ps
= 0.2 ns per saved frame
```

A complete 100 ns trajectory therefore has approximately:

```text
100 ns ÷ 0.2 ns/frame = 500 frames
```

The simulation evaluates millions of integration steps, but the DCD only stores selected snapshots.

---

# 3. What a Trajectory Is

A molecular-dynamics trajectory is an ordered sequence of atomic-coordinate snapshots:

\[
X_0, X_1, X_2, \ldots, X_n
\]

Each frame contains the positions of all atoms at one saved simulation time.

Conceptually:

```text
Frame 1      Frame 2      Frame 3                Frame 500
0.2 ns       0.4 ns       0.6 ns      ...        100.0 ns
```

The trajectory is similar to a movie:

- the atoms are the actors;
- each trajectory frame is one image;
- the full sequence shows how the system changes over time.

Because only one frame is saved every 0.2 ns, very fast events that occur and reverse between saved frames may not appear in the DCD.

---

# 4. PDB and DCD Files

The analysis loads two files:

```python
u = mda.Universe(topology, trajectory)
```

Typically:

- `topology` is a PDB file;
- `trajectory` is a DCD file.

## 4.1 PDB file

The PDB supplies structural identity information, including:

- atom names;
- residue names;
- residue numbers;
- protein and ligand membership;
- atom ordering;
- initial structural metadata.

A useful analogy is:

> The PDB is the cast list.

It tells MDAnalysis which coordinate belongs to which atom.

## 4.2 DCD file

The DCD contains the changing atomic coordinates across saved frames.

It may also include periodic-box dimensions for each frame.

A useful analogy is:

> The DCD is the movie.

The DCD depends on the atom order defined by the topology. If the PDB and DCD atom orders do not match, the analysis is invalid.

## 4.3 Important reference-frame detail

The current code uses the **first saved DCD frame** as the reference frame.

It does not automatically use:

- the docked structure before minimization;
- the original PDB coordinates;
- exact simulation time zero.

The relevant implementation is:

```python
if frame_i == 0:
    ref_ca_centered = ca_pos - c_t
    ref_lig_rel_pbc = lig_rel_pbc.copy()
```

Therefore, the current measurement compares each later frame against the first coordinate frame stored in the DCD.

With a 0.2 ns reporting interval, that first frame is normally near 0.2 ns.

This behavior must be stated clearly in code documentation and scientific writing.

---

# 5. How MDAnalysis Reads the Trajectory

## 5.1 Creating a Universe

MDAnalysis combines the topology and trajectory into a `Universe`:

```python
u = mda.Universe(topology, trajectory)
```

The Universe provides:

- atom identities from the topology;
- coordinates from the trajectory;
- frame iteration;
- atom selections;
- periodic-box information.

## 5.2 Atom selections

The current code selects protein Cα atoms:

```python
ca = u.select_atoms("protein and name CA")
```

This produces one backbone reference atom per protein residue.

The ligand heavy atoms are selected using:

```python
lig = u.select_atoms(
    f"resname {LIGAND_RESNAME} and not name H*"
)
```

The ligand measurement therefore excludes hydrogen atoms.

## 5.3 Iterating through frames

The code loops through the trajectory:

```python
for frame_i, ts in enumerate(u.trajectory):
```

At each iteration:

```python
ca_pos = ca.positions
lig_pos = lig.positions
```

contain the coordinates for the current frame.

The atom selections stay the same, but their coordinates update as MDAnalysis advances through the DCD.

---

# 6. Why Raw Coordinates Cannot Be Compared Directly

The raw trajectory coordinates are physically valid, but they are not automatically suitable for measuring ligand motion relative to the protein.

Three effects must be considered.

## 6.1 Whole-system translation

The protein and ligand can drift together through the water box.

Example:

```text
Earlier frame:
[protein + ligand]

Later frame:
                         [protein + ligand]
```

The ligand may not have moved relative to the pocket, even though its raw Cartesian coordinates changed substantially.

## 6.2 Whole-system rotation

The complete protein–ligand complex can rotate in the simulation box.

If the ligand remains fixed in the pocket while the whole system rotates, raw coordinate comparison would incorrectly report ligand movement.

## 6.3 Periodic-boundary wrapping

The simulation box repeats periodically.

A ligand near one side of the box can be written on the opposite side in a later frame:

```text
Frame A:
|                           ligand |

Frame B:
| ligand                            |
```

The physical motion can be small even though the raw coordinate jump appears large.

These effects are why the analysis must create a **protein-relative and periodic-boundary-consistent coordinate system**.

---

# 7. The Scientific Coordinate System

The simulation box provides a laboratory-style coordinate system:

\[
x,\ y,\ z
\]

but the scientific question is not:

> Where is the ligand inside the simulation box?

The scientific question is:

> Where is the ligand relative to the protein pocket?

Therefore, the protein defines the reference coordinate system.

The analysis removes:

- protein translation;
- protein rotation;
- artificial periodic-image jumps.

It preserves:

- ligand translation relative to the protein;
- ligand rotation relative to the protein;
- ligand internal conformational change.

---

# 8. Corrected Ligand-Pose Calculation

The current implementation is in the `pose_traces` function.

For every frame, it calculates several related measurements:

- corrected ligand RMSD;
- corrected ligand-centroid displacement;
- PBC-naive RMSD control;
- ligand-self-aligned RMSD control.

---

## 8.1 Step 1 — calculate the protein Cα centroid

```python
c_t = ca_pos.mean(axis=0)
```

This is the geometric center of the protein Cα atoms:

\[
\mathbf c_t =
\frac{1}{N_{\mathrm{CA}}}
\sum_i \mathbf r_i(t)
\]

The centroid is used to remove whole-protein translation.

---

## 8.2 Step 2 — express coordinates relative to the protein center

```python
lig_rel_raw = lig_pos - c_t
ca_centered = ca_pos - c_t
```

This moves the protein Cα centroid to the origin.

After this step:

- shared translation of the protein and ligand is removed;
- ligand position is expressed relative to the protein center.

---

## 8.3 Step 3 — apply minimum-image PBC correction

```python
lig_rel_pbc = minimize_vectors(lig_rel_raw, box=box)
```

This converts each ligand-to-protein-center vector to its nearest periodic representation.

The purpose is to prevent a ligand crossing a periodic box boundary from appearing to teleport across the box.

A precise implementation description is:

> Apply the minimum-image convention to ligand coordinates expressed relative to the protein Cα centroid.

This is slightly more specific than the general phrase “PBC reconstruction.”

---

## 8.4 Step 4 — determine the protein rotation

```python
R, _ = rotation_matrix(ca_centered, ref_ca_centered)
```

The Kabsch procedure finds the rotation that best superposes the current protein Cα structure onto the reference protein Cα structure.

This removes whole-protein rigid-body rotation.

The protein is used to define the frame because the scientific question concerns ligand movement relative to the protein.

---

## 8.5 Step 5 — apply the protein-derived rotation to the ligand

```python
aligned_pbc = lig_rel_pbc @ R.T
```

The ligand is not aligned to itself.

Instead, the rotation calculated from the protein is applied unchanged to the ligand.

This is essential because ligand self-alignment would remove ligand translation and orientation changes relative to the pocket.

After this step:

- protein translation is removed;
- protein rotation is removed;
- periodic wrapping is corrected;
- true ligand motion relative to the protein remains.

---

## 8.6 Step 6 — calculate corrected ligand RMSD

```python
rmsd = np.sqrt(
    np.mean(
        np.sum(
            (aligned_pbc - ref_lig_rel_pbc) ** 2,
            axis=1
        )
    )
)
```

Mathematically:

\[
d_{\mathrm{RMSD}}(t)
=
\sqrt{
\frac{1}{N_L}
\sum_{a=1}^{N_L}
\left\|
\widetilde{\mathbf r}_a(t)
-
\widetilde{\mathbf r}_a(0)
\right\|^2
}
\]

where:

- \(N_L\) is the number of ligand heavy atoms;
- \(\widetilde{\mathbf r}_a(t)\) is the corrected protein-relative coordinate of ligand atom \(a\);
- frame \(0\) is the first saved DCD frame.

The output unit is ångströms.

---

# 9. What Corrected Ligand RMSD Measures

Corrected ligand RMSD captures the combined effect of:

1. ligand translation relative to the protein;
2. ligand rotation or reorientation;
3. ligand internal conformational change.

It is therefore a **pose-deviation measurement**.

It does not separately identify which type of motion caused the increase.

Examples:

| Physical change | Corrected RMSD response |
|---|---|
| Protein and ligand translate together | Remains low |
| Protein and ligand rotate together | Remains low |
| Ligand translates within or out of pocket | Increases |
| Ligand rotates inside pocket | Increases |
| Ligand changes internal conformation | Increases |

---

# 10. Corrected Ligand-Centroid Displacement

The code also calculates:

```python
disp_corrected.append(
    np.linalg.norm(
        aligned_pbc.mean(axis=0)
        - ref_lig_rel_pbc.mean(axis=0)
    )
)
```

Mathematically:

\[
d_{\mathrm{disp}}(t)
=
\left\|
\mathbf C_L(t)
-
\mathbf C_L(0)
\right\|
\]

where \(\mathbf C_L\) is the unweighted geometric centroid of the selected ligand heavy atoms.

This mainly measures translation of the ligand as a whole.

## RMSD versus centroid displacement

| Measurement | Main interpretation |
|---|---|
| Corrected centroid displacement | Whole-ligand translation |
| Corrected ligand RMSD | Translation + rotation + internal deformation |

A ligand can have:

- low centroid displacement but high RMSD if it rotates in place;
- high displacement and high RMSD if it moves away;
- low values for both if it remains near its original pose.

Because both use an unweighted average of selected heavy atoms, the document should use the term **geometric centroid**, not mass-weighted center of mass.

---

# 11. Checkpoint-Window Feature Calculation

The per-frame corrected RMSD trace is converted into a checkpoint feature by averaging over the final two nanoseconds before a checkpoint.

For the 20 ns checkpoint:

```python
window_mean(
    pose["time_ns"],
    pose["rmsd_corrected"],
    20.0
)
```

The mask is:

```python
(times > 18.0) & (times <= 20.0)
```

Therefore:

\[
\overline d_{\mathrm{RMSD}}^{18-20}
=
\frac{1}{N_W}
\sum_{t \in (18,20]}
d_{\mathrm{RMSD}}(t)
\]

With a 0.2 ns frame interval, this normally includes approximately:

```text
18.2, 18.4, 18.6, ..., 20.0 ns
```

The average is less sensitive to one unusual frame than a single-frame measurement.

The field name:

```text
same_traj_pose_rmsd_18_20_A
```

should be interpreted as:

> Mean protein-relative, PBC-corrected ligand heavy-atom RMSD from the first saved trajectory frame, averaged over 18–20 ns, in ångströms.

---

# 12. The Measurement Is Not the Model Score

Several related values must not be confused.

## 12.1 Physical observable

```text
Corrected ligand RMSD in Å
```

This is a direct geometric measurement.

## 12.2 AUROC ranking score

For AUROC, the analysis may use:

```python
score = -corrected_rmsd
```

The sign is reversed because lower RMSD indicates greater likelihood of late pose retention.

This sign change is only for ranking direction. It does not change the physical measurement.

## 12.3 Logistic model output

A logistic model may transform one or more standardized measurements into a fitted score.

That output is:

- model-dependent;
- not itself a physical distance;
- not automatically a calibrated probability.

The software should distinguish these fields clearly:

```text
corrected_rmsd_A
centroid_displacement_A
model_score
decision
```

Avoid naming all of them simply `score`.

---

# 13. Same-Trajectory Forecasting Scope

The early measurement is used to forecast the late behavior of the **same continuously running trajectory**.

Conceptually:

```text
0 ns                  20 ns                         70–100 ns
|----------------------|-------------------------------|
     early prefix              forecast gap              late outcome
```

The question is:

> Based on the trajectory observed so far, is this same run likely to retain its ligand pose later?

It is not asking:

> Will a separately launched simulation from the same starting structure behave the same way?

The current results indicate that the signal is much stronger for the same running trajectory than for a separately launched realization.

The interpretation is therefore:

> trajectory-conditioned pose persistence

and not:

- universal molecular stability;
- replica-level reproducibility;
- binding free energy;
- a formal committor;
- biological efficacy.

---

# 14. Control Measurements

The code calculates deliberately weaker alternatives to show why the coordinate definition matters.

## 14.1 Protein-aligned but PBC-naive RMSD

```python
aligned_raw = lig_rel_raw @ R.T
```

This removes protein translation and rotation but skips minimum-image correction.

A periodic box jump can therefore appear as a large ligand movement.

Purpose:

> Isolate the value contributed by PBC correction.

## 14.2 Ligand-self-aligned RMSD

```python
lig_centered = lig_pos - lig_centroid
R_lig, _ = rotation_matrix(lig_centered, ref_lig_centered)
aligned_lig = lig_centered @ R_lig.T
```

This centers the ligand on its own centroid and rotates it to best match its reference shape.

It removes:

- ligand translation relative to the pocket;
- ligand rigid-body rotation relative to the pocket.

It mainly preserves internal ligand deformation.

A rigid ligand could leave the pocket and still produce low ligand-self-aligned RMSD.

Therefore, ligand-self-aligned RMSD answers:

> Did the ligand change its internal shape?

It does not answer:

> Did the ligand change pose relative to the protein?

---

# 15. Why the Core Measurement Is Physically Reasonable

The measurement matches the scientific question because:

- the protein defines the reference frame;
- whole-system drift is removed;
- whole-system rotation is removed;
- periodic wrapping is corrected;
- ligand motion relative to the protein remains;
- ligand atomic-pose change is quantified directly.

The implementation is scientifically appropriate for measuring:

> Protein-relative ligand-pose deviation within one MD trajectory.

The measurement should not be described as direct binding stability.

---

# 16. Important Implementation Risks and Required Checks

The core approach is sound, but production software should explicitly test the following issues.

## 16.1 Molecule wholeness under PBC

The current code applies `minimize_vectors` to ligand-atom vectors independently.

For a compact ligand this will often work correctly, but if the ligand itself is split across the periodic boundary, different atoms could theoretically be mapped into inconsistent images.

A more robust production workflow is:

1. make the protein whole;
2. make the ligand whole using molecular connectivity;
3. place the whole ligand in the periodic image nearest the protein or pocket;
4. center and align the protein;
5. apply the same transform to the ligand.

## 16.2 Whole-protein versus local-pocket alignment

The current coordinate frame uses all protein Cα atoms:

```text
protein and name CA
```

This defines ligand motion relative to the whole protein.

For a large, flexible transporter, domain motion may influence the result.

A useful sensitivity analysis should compare:

- all-protein Cα alignment;
- local-pocket backbone alignment;
- stable-core alignment.

The main method is not wrong, but its reference frame must be explicit.

## 16.3 Reference-frame identity

The current reference is the first saved DCD frame.

The code and documentation should not casually call this exact simulation time zero unless verified.

Recommended metadata:

```text
reference_frame_index
reference_frame_time_ns
reference_source
```

## 16.4 Time calculation

The current code calculates time using:

```python
t = (frame_i + 1) * frame_dt_ns
```

This assumes the first saved frame occurs one reporting interval after time zero.

A more robust implementation should compare this value with the timestamp stored in the trajectory:

```python
ts.time
```

The selected approach should be covered by a unit test.

## 16.5 Atom-order consistency

The PDB and DCD must contain the same atoms in the same order.

The software should verify:

- number of atoms;
- ligand selection size;
- protein Cα count;
- expected residue name;
- no missing coordinate arrays.

## 16.6 Unit consistency

The analysis and output should explicitly record units:

- time in ns;
- coordinates in Å;
- RMSD in Å;
- centroid displacement in Å.

---

# 17. Required Synthetic Unit Tests

These tests should be part of the core repository.

## Test 1 — shared translation invariance

Translate the protein and ligand by the same vector.

Expected result:

```text
corrected RMSD ≈ 0 Å
```

## Test 2 — shared rotation invariance

Rotate the whole protein–ligand complex together.

Expected result:

```text
corrected RMSD ≈ 0 Å
```

## Test 3 — ligand-only translation

Translate a rigid ligand 5 Å relative to the protein.

Expected result:

```text
centroid displacement ≈ 5 Å
corrected RMSD ≈ 5 Å
```

## Test 4 — ligand-only rotation

Rotate the ligand within the pocket without translating its centroid.

Expected result:

```text
centroid displacement remains low
corrected RMSD increases
```

## Test 5 — periodic-image equivalence

Write the same physical ligand position using a different periodic image.

Expected result:

```text
PBC-corrected RMSD remains nearly unchanged
PBC-naive RMSD becomes artificially large
```

## Test 6 — ligand-self-alignment failure case

Move a rigid ligand away from the pocket without changing its internal shape.

Expected result:

```text
protein-relative corrected RMSD increases
ligand-self-aligned RMSD remains low
```

## Test 7 — reference-frame timestamp

Verify that:

- the recorded frame index;
- `ts.time`;
- computed checkpoint windows;

all refer to the intended simulation times.

---

# 18. Recommended Software Architecture

A clean PoseGate implementation should separate the following responsibilities.

```text
posegate/
├── io/
│   ├── topology.py
│   ├── trajectory.py
│   └── metadata.py
├── preprocess/
│   ├── pbc.py
│   ├── alignment.py
│   └── selections.py
├── features/
│   ├── pose_rmsd.py
│   ├── displacement.py
│   └── pocket_retention.py
├── checkpoints/
│   ├── windows.py
│   └── validation.py
├── decision/
│   ├── thresholds.py
│   └── model.py
├── reports/
│   ├── console.py
│   ├── csv.py
│   └── plots.py
└── tests/
    ├── test_translation_invariance.py
    ├── test_rotation_invariance.py
    ├── test_pbc_equivalence.py
    ├── test_ligand_motion.py
    └── test_time_windows.py
```

The physical measurement layer should remain independent from the decision-model layer.

---

# 19. Recommended Core API

A useful internal API could look like:

```python
result = analyze_trajectory(
    topology="complex.pdb",
    trajectory="production.dcd",
    ligand_selection="resname MIL and not name H*",
    protein_selection="protein and name CA",
    reference_frame=0,
    checkpoint_ns=20.0,
    window_ns=2.0,
)
```

Suggested result object:

```python
{
    "reference_frame_index": 0,
    "reference_time_ns": 0.2,
    "checkpoint_ns": 20.0,
    "window_start_ns": 18.0,
    "window_end_ns": 20.0,
    "corrected_rmsd_A": 2.41,
    "centroid_displacement_A": 1.76,
    "pocket_retention": 0.72,
    "model_score": 0.81,
    "decision": "CONTINUE",
    "scope": "same_trajectory_only",
}
```

---

# 20. Recommended User-Facing Output

The software should lead with the physical observable.

```text
PoseGate Analysis
-----------------
Reference frame:              first saved DCD frame
Reference time:               0.2 ns
Checkpoint:                   20.0 ns
Averaging window:             18.0–20.0 ns

Corrected ligand RMSD:        2.41 Å
Centroid displacement:        1.76 Å
Pocket retention:             0.72

Recommendation: CONTINUE

Interpretation:
The ligand remains relatively close to its original
protein-relative pose in this same running trajectory.

Scope:
This result forecasts continuation of this trajectory.
It does not predict binding affinity or a separate replica.
```

---

# 21. Naming Conventions

Use names that preserve meaning and units.

Recommended:

```text
pose_rmsd_corrected_A
pose_rmsd_mean_18_20_A
ligand_centroid_displacement_A
pocket_retention_fraction
same_trajectory_model_score
stop_continue_decision
```

Avoid ambiguous names such as:

```text
score
stability
binding_score
probability
```

unless their exact meaning is defined.

---

# 22. Core Interpretation

The central measurement should be described as:

> The mean heavy-atom RMS deviation of the ligand from its first saved protein-relative pose, after minimum-image periodic correction and protein Cα alignment, averaged over the selected checkpoint window.

The central same-trajectory claim should be described as:

> A correctly processed early ligand-pose measurement can forecast whether that same running trajectory later retains its pose.

The mechanism should be described as:

> Trajectory-conditioned pose persistence.

The project should not claim that the measurement proves:

- binding;
- affinity;
- inhibition;
- thermodynamic stability;
- replica-level reproducibility;
- a universal reaction coordinate.

---

# 23. Development Checklist

Before declaring the core implementation stable:

- [ ] Confirm PDB and DCD atom ordering.
- [ ] Confirm ligand heavy-atom selection.
- [ ] Confirm protein Cα selection.
- [ ] Verify periodic box data exist for every analyzed frame.
- [ ] Make the ligand whole before image placement.
- [ ] Audit whether the protein also needs unwrapping.
- [ ] Confirm first saved frame time.
- [ ] Compare computed time with `ts.time`.
- [ ] Add translation-invariance test.
- [ ] Add rotation-invariance test.
- [ ] Add periodic-image-equivalence test.
- [ ] Add ligand-only translation test.
- [ ] Add ligand-only rotation test.
- [ ] Add ligand-self-alignment failure test.
- [ ] Compare whole-protein and pocket-local alignment.
- [ ] Record all output units.
- [ ] Separate physical measurements from model scores.
- [ ] Record software version and Git commit.
- [ ] Save run-level metadata and analysis configuration.
- [ ] Preserve same-trajectory scope in all reports.

---

# 24. Final Technical Summary

OpenMM generates a time-dependent molecular system by integrating atomic motion under a molecular-mechanics force field.

The DCD stores selected coordinate snapshots from that simulation, while the PDB provides atom identities and ordering.

MDAnalysis combines those files, selects the protein Cα atoms and ligand heavy atoms, and processes every saved frame.

Raw coordinates cannot be compared directly because the complete system can translate, rotate, and wrap across periodic boundaries.

The analysis therefore:

1. expresses ligand coordinates relative to the protein Cα centroid;
2. applies minimum-image periodic correction;
3. aligns the current protein Cα coordinates to the first saved trajectory frame;
4. applies the protein-derived rotation to the ligand;
5. calculates ligand heavy-atom RMSD and centroid displacement;
6. averages those values over a checkpoint window.

The resulting corrected ligand RMSD is a physically meaningful measure of ligand-pose deviation relative to the protein in the same trajectory.

It is not inherently a machine-learning score, binding affinity, or universal stability measurement.

Its value comes from answering a narrow operational question:

> Has the ligand already developed a pose pattern that is informative about the later continuation of this same running MD trajectory?
