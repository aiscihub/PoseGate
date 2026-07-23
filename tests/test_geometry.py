import numpy as np
import pytest

from posegate.geometry import (
    frame_relative_geometry,
    protein_relative_coordinates,
)


RNG = np.random.default_rng(20260723)
BOX = np.array([100.0, 100.0, 100.0, 90.0, 90.0, 90.0])


def rotation_matrix(angle_degrees: float, axis: np.ndarray) -> np.ndarray:
    axis = np.asarray(axis, dtype=float)
    axis /= np.linalg.norm(axis)
    theta = np.radians(angle_degrees)
    x, y, z = axis
    cross = np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])
    return np.eye(3) + np.sin(theta) * cross + (1 - np.cos(theta)) * (cross @ cross)


def system() -> tuple[np.ndarray, np.ndarray]:
    protein = RNG.normal(scale=5.0, size=(80, 3)) + 50.0
    ligand = RNG.normal(scale=1.2, size=(12, 3)) + np.array([54.0, 50.0, 50.0])
    return protein, ligand


def measure(
    protein0: np.ndarray,
    ligand0: np.ndarray,
    protein_t: np.ndarray,
    ligand_t: np.ndarray,
):
    frame0, relative0 = protein_relative_coordinates(protein0, ligand0, BOX)
    frame_t, relative_t = protein_relative_coordinates(protein_t, ligand_t, BOX)
    return frame_relative_geometry(frame_t, frame0, relative_t, relative0)


def test_shared_translation_produces_zero_corrected_motion() -> None:
    protein, ligand = system()
    shift = np.array([3.0, -2.0, 1.0])
    result = measure(protein, ligand, protein + shift, ligand + shift)
    assert result.corrected_pose_rmsd_angstrom == pytest.approx(0.0, abs=1e-5)
    assert result.centroid_displacement_angstrom == pytest.approx(0.0, abs=1e-5)
    assert np.degrees(result.frame_relative_reorientation_radians) == pytest.approx(
        0.0, abs=1e-4
    )


def test_shared_rotation_produces_zero_corrected_motion() -> None:
    protein, ligand = system()
    rotation = rotation_matrix(25.0, np.array([0.3, 0.5, 0.8]))
    origin = np.vstack([protein, ligand]).mean(axis=0)
    protein_t = (protein - origin) @ rotation.T + origin
    ligand_t = (ligand - origin) @ rotation.T + origin
    result = measure(protein, ligand, protein_t, ligand_t)
    assert result.corrected_pose_rmsd_angstrom == pytest.approx(0.0, abs=2e-5)
    assert result.centroid_displacement_angstrom == pytest.approx(0.0, abs=2e-5)
    assert np.degrees(result.frame_relative_reorientation_radians) == pytest.approx(
        0.0, abs=1e-4
    )


def test_ligand_translation_is_recovered() -> None:
    protein, ligand = system()
    shift = np.array([1.0, -2.0, 0.5])
    result = measure(protein, ligand, protein, ligand + shift)
    expected = float(np.linalg.norm(shift))
    assert result.corrected_pose_rmsd_angstrom == pytest.approx(expected, abs=1e-5)
    assert result.centroid_displacement_angstrom == pytest.approx(expected, abs=1e-5)
    assert result.internal_deformation_angstrom == pytest.approx(0.0, abs=1e-5)


def test_ligand_only_rotation_recovers_angle() -> None:
    protein, ligand = system()
    rotation = rotation_matrix(12.0, np.array([0.2, 0.7, 0.1]))
    center = ligand.mean(axis=0)
    ligand_t = (ligand - center) @ rotation.T + center
    result = measure(protein, ligand, protein, ligand_t)
    assert result.centroid_displacement_angstrom == pytest.approx(0.0, abs=1e-5)
    assert np.degrees(result.frame_relative_reorientation_radians) == pytest.approx(
        12.0, abs=1e-4
    )
    assert result.internal_deformation_angstrom == pytest.approx(0.0, abs=1e-5)


def test_periodic_image_shift_does_not_change_measurement() -> None:
    protein, ligand = system()
    shift = np.array([0.5, -0.25, 0.75])
    baseline = measure(protein, ligand, protein, ligand + shift)
    imaged = measure(
        protein,
        ligand,
        protein,
        ligand + shift + np.array([BOX[0], 0.0, 0.0]),
    )
    assert imaged.corrected_pose_rmsd_angstrom == pytest.approx(
        baseline.corrected_pose_rmsd_angstrom, abs=1e-5
    )
    assert imaged.centroid_displacement_angstrom == pytest.approx(
        baseline.centroid_displacement_angstrom, abs=1e-5
    )
