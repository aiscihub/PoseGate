#!/usr/bin/env python3
"""Synthetic unit tests for frame_relative_theta (script 35).

Pure-geometry tests with no trajectory/topology files -- constructs
synthetic point clouds with KNOWN, programmatically-applied rigid-body
motions and checks that the actual production function
(`35_local_pocket_alignment_control.frame_relative_theta`) recovers the
expected T/Theta in each case, within numerical tolerance. This directly
exercises the real code path (not a reimplementation), so it is a
correctness check on the actual pipeline, not just an argument about it.

Six cases:
  A. Domain motion only (pocket+ligand rotate together, ligand fixed
     relative to pocket): whole-frame Theta should register ~domain
     angle; pocket-frame Theta should be ~0.
  B. Ligand-only rotation (protein entirely fixed): both whole- and
     pocket-frame Theta should recover the ligand's true rotation angle.
  C. Shared translation of the entire system: T and Theta should be ~0 in
     both frames (translation is fully removed by centroid-relative
     alignment).
  D. Shared/global rotation of the entire system (protein+ligand rotate
     together, ligand fixed relative to protein): T and Theta should be
     ~0 in both frames. This is the direct regression test for the bug
     fixed in script 35 -- script 33's original Theta (a pure ligand-self
     Kabsch fit in raw lab coordinates, no protein reference) would
     INCORRECTLY report Theta ~= the global rotation angle here, since it
     never removed the complex's own tumbling. Both the naive and the
     fixed calculation are run so the contrast is explicit, not asserted.
  E. Periodic-image equivalence: shifting the ligand's raw stored
     coordinates by an integer multiple of a box vector (simulating that
     the MD writer stored a different periodic image) must not change T
     or Theta once minimum-image correction is applied.
  F. Small-angle self-consistency: a small known rotation should recover
     with correspondingly small residual error, confirming the fit is not
     saturating or degenerate at small angles (relevant since real Theta
     values are mostly in the ~5-25 degree range).

Run: python test_synthetic_frame_relative_theta.py
Exits non-zero if any assertion fails.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import minimize_vectors

sys.path.insert(0, str(Path(__file__).resolve().parent))
import importlib

control_mod = importlib.import_module("35_local_pocket_alignment_control")
frame_relative_theta = control_mod.frame_relative_theta

RNG = np.random.default_rng(20260718)
TOL_DEG = 0.5  # numerical tolerance for exact synthetic cases (noiseless Kabsch fits)


def random_rotation_matrix(angle_deg: float, axis: np.ndarray | None = None) -> np.ndarray:
    if axis is None:
        axis = RNG.normal(size=3)
    axis = axis / np.linalg.norm(axis)
    theta = np.radians(angle_deg)
    kx, ky, kz = axis
    k = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    identity = np.eye(3)
    return identity + np.sin(theta) * k + (1 - np.cos(theta)) * (k @ k)


def make_system():
    non_pocket = RNG.normal(scale=15.0, size=(170, 3)) + np.array([0.0, 0.0, 0.0])
    pocket = RNG.normal(scale=4.0, size=(30, 3)) + np.array([20.0, 0.0, 0.0])
    ligand = RNG.normal(scale=2.0, size=(15, 3)) + np.array([22.0, 0.0, 0.0])
    return non_pocket, pocket, ligand


def theta_from_whole_and_pocket(non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t, box=None):
    whole0 = np.vstack([non_pocket0, pocket0])
    wholet = np.vstack([non_pocket_t, pocket_t])

    whole_center0 = whole0.mean(axis=0)
    whole_centert = wholet.mean(axis=0)
    whole0_centered = whole0 - whole_center0
    wholet_centered = wholet - whole_center0  # centered on frame-0 style below is redone per-call

    pocket_center0 = pocket0.mean(axis=0)
    pocket_centert = pocket_t.mean(axis=0)

    if box is not None:
        ligand_rel0_whole = minimize_vectors(ligand0 - whole_center0, box=box)
        ligand_relt_whole = minimize_vectors(ligand_t - whole_centert, box=box)
        ligand_rel0_pocket = minimize_vectors(ligand0 - pocket_center0, box=box)
        ligand_relt_pocket = minimize_vectors(ligand_t - pocket_centert, box=box)
    else:
        ligand_rel0_whole = ligand0 - whole_center0
        ligand_relt_whole = ligand_t - whole_centert
        ligand_rel0_pocket = ligand0 - pocket_center0
        ligand_relt_pocket = ligand_t - pocket_centert

    t_whole, theta_whole, _ = frame_relative_theta(
        wholet - whole_centert, whole0 - whole_center0, ligand_relt_whole, ligand_rel0_whole
    )
    t_pocket, theta_pocket, _ = frame_relative_theta(
        pocket_t - pocket_centert, pocket0 - pocket_center0, ligand_relt_pocket, ligand_rel0_pocket
    )
    return t_whole, np.degrees(theta_whole), t_pocket, np.degrees(theta_pocket)


def naive_lab_frame_theta(ligand0, ligand_t):
    """Reproduces script 33's ORIGINAL (buggy) Theta: pure ligand-self
    Kabsch in raw lab coordinates, no protein/pocket reference at all."""
    c0 = ligand0 - ligand0.mean(axis=0)
    ct = ligand_t - ligand_t.mean(axis=0)
    q_lig, _ = rotation_matrix(ct, c0)
    return float(np.degrees(np.arccos(np.clip((np.trace(q_lig) - 1.0) / 2.0, -1.0, 1.0))))


def check(name: str, actual: float, expected: float, tol: float = TOL_DEG) -> bool:
    ok = abs(actual - expected) <= tol
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {name}: actual={actual:.4f}  expected~{expected:.4f}  (tol={tol})")
    return ok


def main() -> None:
    results = []

    # --- Case A: domain motion only, ligand fixed relative to pocket ---
    print("Case A: domain+ligand rotate together by 15 deg; rest of protein fixed.")
    non_pocket0, pocket0, ligand0 = make_system()
    domain_angle = 15.0
    axis = np.array([0.0, 0.0, 1.0])
    r_domain = random_rotation_matrix(domain_angle, axis)
    pocket_center0 = pocket0.mean(axis=0)
    pocket_t = (pocket0 - pocket_center0) @ r_domain.T + pocket_center0
    ligand_t = (ligand0 - pocket_center0) @ r_domain.T + pocket_center0  # ligand rigidly fixed to pocket
    non_pocket_t = non_pocket0.copy()  # rest of protein does not move

    t_whole, theta_whole, t_pocket, theta_pocket = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t
    )
    results.append(check("A: whole-frame Theta registers domain rotation", theta_whole, domain_angle, tol=2.0))
    results.append(check("A: pocket-frame Theta ~ 0 (ligand fixed relative to pocket)", theta_pocket, 0.0))

    # --- Case B: ligand-only rotation, protein entirely fixed ---
    print("\nCase B: only ligand rotates by 12 deg; protein (whole and pocket) fixed.")
    non_pocket0, pocket0, ligand0 = make_system()
    lig_angle = 12.0
    r_lig = random_rotation_matrix(lig_angle)
    ligand_center0 = ligand0.mean(axis=0)
    ligand_t = (ligand0 - ligand_center0) @ r_lig.T + ligand_center0
    non_pocket_t = non_pocket0.copy()
    pocket_t = pocket0.copy()

    t_whole, theta_whole, t_pocket, theta_pocket = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t
    )
    results.append(check("B: whole-frame Theta recovers ligand rotation", theta_whole, lig_angle))
    results.append(check("B: pocket-frame Theta recovers ligand rotation", theta_pocket, lig_angle))

    # --- Case C: shared translation of everything ---
    print("\nCase C: entire system (protein+ligand) translates by a fixed vector.")
    non_pocket0, pocket0, ligand0 = make_system()
    shift = np.array([37.0, -12.0, 5.0])
    non_pocket_t = non_pocket0 + shift
    pocket_t = pocket0 + shift
    ligand_t = ligand0 + shift

    t_whole, theta_whole, t_pocket, theta_pocket = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t
    )
    results.append(check("C: whole-frame T ~ 0 under shared translation", t_whole, 0.0, tol=1e-6))
    results.append(check("C: whole-frame Theta ~ 0 under shared translation", theta_whole, 0.0))
    results.append(check("C: pocket-frame T ~ 0 under shared translation", t_pocket, 0.0, tol=1e-6))
    results.append(check("C: pocket-frame Theta ~ 0 under shared translation", theta_pocket, 0.0))

    # --- Case D: shared/global rotation of everything (the bug-fix regression test) ---
    print("\nCase D: entire complex (protein+ligand) tumbles together by 25 deg; ligand does not move relative to protein.")
    non_pocket0, pocket0, ligand0 = make_system()
    global_angle = 25.0
    r_global = random_rotation_matrix(global_angle)
    origin = np.vstack([non_pocket0, pocket0, ligand0]).mean(axis=0)
    non_pocket_t = (non_pocket0 - origin) @ r_global.T + origin
    pocket_t = (pocket0 - origin) @ r_global.T + origin
    ligand_t = (ligand0 - origin) @ r_global.T + origin

    t_whole, theta_whole, t_pocket, theta_pocket = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t
    )
    naive_theta = naive_lab_frame_theta(ligand0, ligand_t)
    results.append(check("D: whole-frame Theta ~ 0 under global tumbling (fixed method)", theta_whole, 0.0))
    results.append(check("D: pocket-frame Theta ~ 0 under global tumbling (fixed method)", theta_pocket, 0.0))
    print(
        f"  [INFO] D: naive lab-frame Theta (script 33's ORIGINAL method) = {naive_theta:.2f} deg"
        f"  -- incorrectly reports ~the global tumbling angle ({global_angle} deg); "
        f"this is the bug the frame-relative fix in script 35 corrects."
    )
    results.append(check("D: naive method WOULD have misreported ~global angle (confirms bug existed)", naive_theta, global_angle, tol=2.0))

    # --- Case E: periodic-image equivalence ---
    print("\nCase E: ligand's raw stored coordinates shifted by one box vector (simulated periodic image).")
    non_pocket0, pocket0, ligand0 = make_system()
    lig_angle = 8.0
    r_lig = random_rotation_matrix(lig_angle)
    ligand_center0 = ligand0.mean(axis=0)
    ligand_t = (ligand0 - ligand_center0) @ r_lig.T + ligand_center0
    non_pocket_t = non_pocket0.copy()
    pocket_t = pocket0.copy()
    box = np.array([80.0, 80.0, 80.0, 90.0, 90.0, 90.0])
    box_vector = np.array([80.0, 0.0, 0.0])

    t_whole_a, theta_whole_a, t_pocket_a, theta_pocket_a = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t, box=box
    )
    t_whole_b, theta_whole_b, t_pocket_b, theta_pocket_b = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t + box_vector, box=box
    )
    results.append(check("E: whole-frame T unchanged by periodic-image shift", t_whole_b, t_whole_a, tol=1e-4))
    results.append(check("E: whole-frame Theta unchanged by periodic-image shift", theta_whole_b, theta_whole_a, tol=1e-3))
    results.append(check("E: pocket-frame T unchanged by periodic-image shift", t_pocket_b, t_pocket_a, tol=1e-4))
    results.append(check("E: pocket-frame Theta unchanged by periodic-image shift", theta_pocket_b, theta_pocket_a, tol=1e-3))

    # --- Case F: small-angle recovery ---
    print("\nCase F: small known ligand rotation (2 deg), protein fixed.")
    non_pocket0, pocket0, ligand0 = make_system()
    small_angle = 2.0
    r_small = random_rotation_matrix(small_angle)
    ligand_center0 = ligand0.mean(axis=0)
    ligand_t = (ligand0 - ligand_center0) @ r_small.T + ligand_center0
    non_pocket_t = non_pocket0.copy()
    pocket_t = pocket0.copy()

    t_whole, theta_whole, t_pocket, theta_pocket = theta_from_whole_and_pocket(
        non_pocket0, pocket0, ligand0, non_pocket_t, pocket_t, ligand_t
    )
    results.append(check("F: whole-frame Theta recovers small ligand rotation", theta_whole, small_angle, tol=0.1))
    results.append(check("F: pocket-frame Theta recovers small ligand rotation", theta_pocket, small_angle, tol=0.1))

    print(f"\n{sum(results)}/{len(results)} checks passed.")
    if not all(results):
        sys.exit(1)


if __name__ == "__main__":
    main()
