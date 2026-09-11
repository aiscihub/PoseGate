"""Reusable coordinate operations extracted from the manuscript analysis."""

from __future__ import annotations

from dataclasses import dataclass

from MDAnalysis.analysis.align import rotation_matrix
from MDAnalysis.lib.distances import minimize_vectors
import numpy as np

from .exceptions import InputValidationError
from .config import SUPPORTED_PBC_METHOD, WHOLE_LIGAND_PBC_METHOD


@dataclass(frozen=True)
class FrameGeometry:
    """Ligand motion relative to a protein-derived reference frame."""

    corrected_pose_rmsd_angstrom: float
    centroid_displacement_angstrom: float
    frame_relative_reorientation_radians: float
    internal_deformation_angstrom: float
    pbc_naive_pose_rmsd_angstrom: float | None = None


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
    *,
    method: str = WHOLE_LIGAND_PBC_METHOD,
) -> tuple[np.ndarray, np.ndarray]:
    """Center the protein and select one periodic image of the whole ligand.

    The explicitly named v1 method is retained only for legacy records.
    Neither method reconstructs molecular components or unwraps the protein.
    """
    frame = _xyz(alignment_positions, "alignment_positions")
    ligand = _xyz(ligand_positions, "ligand_positions")
    periodic_box = validate_box(box)
    if method == WHOLE_LIGAND_PBC_METHOD:
        frame, ligand = frame.astype(float), ligand.astype(float)
    center = frame.mean(axis=0)
    frame_centered = frame - center
    if method == WHOLE_LIGAND_PBC_METHOD:
        ligand_relative, _ = whole_ligand_relative_to_protein(ligand, center, periodic_box)
    elif method == SUPPORTED_PBC_METHOD:
        ligand_relative = minimize_vectors(ligand - center, box=periodic_box)
    else:
        raise InputValidationError(f"unsupported coordinate method: {method}")
    return frame_centered, np.asarray(ligand_relative)


def whole_ligand_relative_to_protein(
    ligand_positions: np.ndarray, protein_center: np.ndarray, box: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Return common-shift coordinates and the atom-wise shift disagreement.

    The diagnostic does not alter coordinates or discard a trajectory.
    """
    ligand = _xyz(ligand_positions, "ligand_positions").astype(float)
    center = np.asarray(protein_center, dtype=float)
    if center.shape != (3,) or not np.isfinite(center).all():
        raise InputValidationError("protein_center must be a finite three-vector")
    periodic_box = validate_box(box)
    relative = ligand - center
    centroid = relative.mean(axis=0)
    shift = minimize_vectors(centroid[None, :], box=periodic_box)[0] - centroid
    atomwise_shifts = minimize_vectors(relative, box=periodic_box) - relative
    disagreement = float(np.linalg.norm(atomwise_shifts - shift, axis=1).max())
    return relative + shift, disagreement


def validate_raw_ligand_wholeness(
    ligand_positions: np.ndarray, box: np.ndarray, *, max_diameter_box_fraction: float = 0.45,
) -> float:
    """Check raw ligand extent before imaging, without attempting reconstruction."""
    ligand = _xyz(ligand_positions, "ligand_positions").astype(float)
    periodic_box = validate_box(box)
    if not 0 < max_diameter_box_fraction < 0.5:
        raise InputValidationError("ligand diameter fraction must lie between 0 and 0.5")
    differences = ligand[:, None, :] - ligand[None, :, :]
    diameter = float(np.sqrt(np.sum(differences * differences, axis=2).max()))
    if diameter > max_diameter_box_fraction * periodic_box[:3].min():
        raise InputValidationError("raw ligand exceeds the declared diameter limit; whole coordinates required")
    return diameter


def frame_relative_geometry(
    frame_centered: np.ndarray,
    reference_frame_centered: np.ndarray,
    ligand_relative: np.ndarray,
    reference_ligand_relative: np.ndarray,
    *,
    raw_ligand_relative: np.ndarray | None = None,
    reference_raw_ligand_relative: np.ndarray | None = None,
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
    naive_rmsd = None
    if (raw_ligand_relative is None) != (reference_raw_ligand_relative is None):
        raise InputValidationError("both raw ligand coordinate arrays are required")
    if raw_ligand_relative is not None:
        raw = _xyz(raw_ligand_relative, "raw_ligand_relative")
        reference_raw = _xyz(reference_raw_ligand_relative, "reference_raw_ligand_relative")
        if raw.shape != ligand.shape or reference_raw.shape != ligand.shape:
            raise InputValidationError("raw and corrected ligand atom counts differ")
        naive_rmsd = float(np.sqrt(np.mean(np.sum((raw @ protein_rotation.T - reference_raw) ** 2, axis=1))))
    return FrameGeometry(
        corrected_pose_rmsd_angstrom=corrected_rmsd,
        centroid_displacement_angstrom=centroid_displacement,
        frame_relative_reorientation_radians=theta,
        internal_deformation_angstrom=internal_deformation,
        pbc_naive_pose_rmsd_angstrom=naive_rmsd,
    )


def zero_geometry() -> FrameGeometry:
    return FrameGeometry(
        corrected_pose_rmsd_angstrom=0.0,
        centroid_displacement_angstrom=0.0,
        frame_relative_reorientation_radians=0.0,
        internal_deformation_angstrom=0.0,
        pbc_naive_pose_rmsd_angstrom=0.0,
    )
