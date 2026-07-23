"""Reusable coordinate operations extracted from the manuscript analysis."""

from __future__ import annotations

from dataclasses import dataclass

from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import minimize_vectors
import numpy as np

from .exceptions import InputValidationError


@dataclass(frozen=True)
class FrameGeometry:
    """Ligand motion relative to a protein-derived reference frame."""

    corrected_pose_rmsd_angstrom: float
    centroid_displacement_angstrom: float
    frame_relative_reorientation_radians: float
    internal_deformation_angstrom: float


def _xyz(value: np.ndarray, name: str) -> np.ndarray:
    array = np.asarray(value)
    if not np.issubdtype(array.dtype, np.floating):
        array = array.astype(float)
    if array.ndim != 2 or array.shape[1] != 3 or array.shape[0] == 0:
        raise InputValidationError(f"{name} must have shape (n_atoms, 3)")
    if not np.isfinite(array).all():
        raise InputValidationError(f"{name} contains non-finite coordinates")
    return array


def validate_box(box: np.ndarray) -> np.ndarray:
    """Return a finite MDAnalysis box or fail instead of silently disabling PBC."""
    array = np.asarray(box)
    if not np.issubdtype(array.dtype, np.floating):
        array = array.astype(float)
    if array.shape != (6,):
        raise InputValidationError(
            "periodic box must contain [lx, ly, lz, alpha, beta, gamma]"
        )
    if not np.isfinite(array).all() or np.any(array[:3] <= 0):
        raise InputValidationError("periodic box lengths must be finite and positive")
    if np.any(array[3:] <= 0) or np.any(array[3:] >= 180):
        raise InputValidationError(
            "periodic box angles must lie between 0 and 180 degrees"
        )
    return array


def protein_relative_coordinates(
    alignment_positions: np.ndarray,
    ligand_positions: np.ndarray,
    box: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Center the protein frame and map ligand atoms to its nearest images.

    This reproduces the validated paper implementation: each ligand atom is
    represented by its minimum-image vector from the protein C-alpha centroid.
    It intentionally does not claim bond-connectivity molecule reconstruction.
    """
    frame = _xyz(alignment_positions, "alignment_positions")
    ligand = _xyz(ligand_positions, "ligand_positions")
    periodic_box = validate_box(box)
    center = frame.mean(axis=0)
    frame_centered = frame - center
    ligand_relative = minimize_vectors(ligand - center, box=periodic_box)
    return frame_centered, np.asarray(ligand_relative)


def frame_relative_geometry(
    frame_centered: np.ndarray,
    reference_frame_centered: np.ndarray,
    ligand_relative: np.ndarray,
    reference_ligand_relative: np.ndarray,
) -> FrameGeometry:
    """Measure ligand pose after alignment by the protein reference frame."""
    frame = _xyz(frame_centered, "frame_centered")
    reference_frame = _xyz(reference_frame_centered, "reference_frame_centered")
    ligand = _xyz(ligand_relative, "ligand_relative")
    reference_ligand = _xyz(reference_ligand_relative, "reference_ligand_relative")
    if frame.shape != reference_frame.shape:
        raise InputValidationError("current and reference alignment atom counts differ")
    if ligand.shape != reference_ligand.shape:
        raise InputValidationError("current and reference ligand atom counts differ")

    protein_rotation, _ = rotation_matrix(frame, reference_frame)
    ligand_aligned = ligand @ protein_rotation.T
    delta = ligand_aligned - reference_ligand
    corrected_rmsd = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
    centroid_displacement = float(np.linalg.norm(delta.mean(axis=0)))

    ligand_centered = ligand_aligned - ligand_aligned.mean(axis=0)
    reference_ligand_centered = reference_ligand - reference_ligand.mean(axis=0)
    ligand_rotation, _ = rotation_matrix(ligand_centered, reference_ligand_centered)
    theta = float(
        np.arccos(np.clip((float(np.trace(ligand_rotation)) - 1.0) / 2.0, -1.0, 1.0))
    )
    ligand_self_aligned = ligand_centered @ ligand_rotation.T
    internal_deformation = float(
        np.sqrt(
            np.mean(
                np.sum(
                    (ligand_self_aligned - reference_ligand_centered) ** 2,
                    axis=1,
                )
            )
        )
    )
    return FrameGeometry(
        corrected_pose_rmsd_angstrom=corrected_rmsd,
        centroid_displacement_angstrom=centroid_displacement,
        frame_relative_reorientation_radians=theta,
        internal_deformation_angstrom=internal_deformation,
    )


def zero_geometry() -> FrameGeometry:
    return FrameGeometry(
        corrected_pose_rmsd_angstrom=0.0,
        centroid_displacement_angstrom=0.0,
        frame_relative_reorientation_radians=0.0,
        internal_deformation_angstrom=0.0,
    )
