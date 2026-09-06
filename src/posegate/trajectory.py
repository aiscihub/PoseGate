"""Strict measurement of corrected geometry from trajectory prefixes."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
from pathlib import Path
from typing import Any
import warnings

import MDAnalysis as mda
import numpy as np

from .config import PoseGateConfig
from .exceptions import (
    CheckpointNotReady,
    InputValidationError,
    OutcomeEmbargoError,
)
from .geometry import (
    FrameGeometry,
    frame_relative_geometry,
    protein_relative_coordinates,
    zero_geometry,
)


TIME_TOLERANCE_NS = 1e-4


@dataclass(frozen=True)
class InputInspection:
    topology_path: str
    trajectory_path: str
    topology_atoms: int
    alignment_atoms: int
    ligand_atoms: int
    frame_count_available: int
    frame_interval_ns: float
    trajectory_end_ns: float
    checkpoint_ns: float
    required_checkpoint_frames: int
    checkpoint_ready: bool
    box_vectors_valid_at_first_frame: bool


@dataclass(frozen=True)
class PrefixSeries:
    """Per-frame geometry for every frame of the measured prefix.

    The scalar means remain the policy input. This series exists so a sealed
    record can be re-plotted without re-reading the trajectory, and it is never
    consulted by :mod:`posegate.policy`.
    """

    time_ns: tuple[float, ...]
    corrected_pose_rmsd_angstrom: tuple[float, ...]
    centroid_displacement_angstrom: tuple[float, ...]
    frame_relative_reorientation_degrees: tuple[float, ...]
    internal_deformation_angstrom: tuple[float, ...]
    in_feature_window: tuple[bool, ...]


@dataclass(frozen=True)
class DomainProfile:
    """Ligand extent relative to the periodic box across the measured prefix.

    A ligand approaching the box size makes minimum-image reconstruction
    ambiguous, so this is the applicability check for the coordinate method
    rather than a property of the pose.
    """

    max_ligand_diameter_angstrom: float
    min_box_length_angstrom: float
    ligand_diameter_box_fraction: float


@dataclass(frozen=True)
class PrefixMeasurement:
    topology_path: str
    trajectory_path: str
    alignment_selection: str
    ligand_selection: str
    alignment_atom_count: int
    ligand_atom_count: int
    trajectory_frames_available: int
    prefix_frames_used: int
    frame_interval_ns: float
    max_available_time_ns: float
    checkpoint_ns: float
    window_start_ns: float
    window_end_ns: float
    window_frame_count: int
    prefix_coordinates_sha256: str
    corrected_pose_rmsd_mean_angstrom: float
    corrected_centroid_displacement_mean_angstrom: float
    frame_relative_reorientation_mean_degrees: float
    internal_deformation_mean_angstrom: float
    domain_profile: DomainProfile
    series: PrefixSeries

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    def policy_measurements(self) -> dict[str, float]:
        return {
            "corrected_pose_rmsd_mean_angstrom": (
                self.corrected_pose_rmsd_mean_angstrom
            ),
            "corrected_centroid_displacement_mean_angstrom": (
                self.corrected_centroid_displacement_mean_angstrom
            ),
        }


@dataclass(frozen=True)
class LateOutcome:
    trajectory_frames: int
    trajectory_end_ns: float
    late_window_start_ns: float
    late_window_end_ns: float
    late_window_frame_count: int
    late_pose_rmsd_median_angstrom: float
    retained_rmsd_threshold_angstrom: float
    pose_retained: bool

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def _universe(topology: Path, trajectory: Path) -> mda.Universe:
    if not topology.is_file():
        raise InputValidationError(f"topology does not exist: {topology}")
    if not trajectory.is_file():
        raise InputValidationError(f"trajectory does not exist: {trajectory}")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return mda.Universe(str(topology), str(trajectory))
    except (OSError, ValueError) as exc:
        raise InputValidationError(
            f"could not read topology/trajectory: {exc}"
        ) from exc


def _frame_interval_ns(universe: mda.Universe) -> float:
    dt_ps = float(getattr(universe.trajectory, "dt", float("nan")))
    dt_ns = dt_ps / 1000.0
    if not np.isfinite(dt_ns) or dt_ns <= 0:
        raise InputValidationError(
            "trajectory frame interval is missing or invalid; no fallback is allowed"
        )
    return dt_ns


def _required_frames(checkpoint_ns: float, dt_ns: float) -> int:
    required = int(round(checkpoint_ns / dt_ns))
    if required <= 0 or abs(required * dt_ns - checkpoint_ns) > TIME_TOLERANCE_NS:
        raise InputValidationError(
            f"checkpoint {checkpoint_ns:g} ns is not representable at "
            f"the trajectory cadence {dt_ns:g} ns"
        )
    return required


def _selections(
    universe: mda.Universe, config: PoseGateConfig
) -> tuple[mda.core.groups.AtomGroup, mda.core.groups.AtomGroup]:
    protein = universe.select_atoms(config.selections.protein)
    alignment = universe.select_atoms(config.selections.alignment_atoms)
    ligand = universe.select_atoms(config.selections.ligand)
    if protein.n_atoms == 0:
        raise InputValidationError(
            f"protein selection is empty: {config.selections.protein}"
        )
    if alignment.n_atoms == 0:
        raise InputValidationError(
            f"alignment selection is empty: {config.selections.alignment_atoms}"
        )
    if ligand.n_atoms == 0:
        raise InputValidationError(
            f"ligand selection is empty: {config.selections.ligand}"
        )
    if not np.isin(alignment.indices, protein.indices).all():
        raise InputValidationError(
            "alignment selection contains atoms outside the protein selection"
        )
    if np.intersect1d(ligand.indices, protein.indices).size:
        raise InputValidationError("ligand and protein selections overlap")
    return alignment, ligand


def inspect_inputs(
    topology: str | Path, trajectory: str | Path, config: PoseGateConfig
) -> InputInspection:
    topology_path = Path(topology).resolve()
    trajectory_path = Path(trajectory).resolve()
    universe = _universe(topology_path, trajectory_path)
    alignment, ligand = _selections(universe, config)
    dt_ns = _frame_interval_ns(universe)
    available = len(universe.trajectory)
    required = _required_frames(config.checkpoint.time_ns, dt_ns)
    universe.trajectory[0]
    box = np.asarray(universe.trajectory.ts.dimensions)
    box_valid = bool(
        box.shape == (6,)
        and np.isfinite(box).all()
        and np.all(box[:3] > 0)
        and np.all(box[3:] > 0)
        and np.all(box[3:] < 180)
    )
    return InputInspection(
        topology_path=str(topology_path),
        trajectory_path=str(trajectory_path),
        topology_atoms=universe.atoms.n_atoms,
        alignment_atoms=alignment.n_atoms,
        ligand_atoms=ligand.n_atoms,
        frame_count_available=available,
        frame_interval_ns=dt_ns,
        trajectory_end_ns=available * dt_ns,
        checkpoint_ns=config.checkpoint.time_ns,
        required_checkpoint_frames=required,
        checkpoint_ready=available >= required,
        box_vectors_valid_at_first_frame=box_valid,
    )


def _measure_trace(
    universe: mda.Universe,
    config: PoseGateConfig,
    *,
    frame_limit: int | None,
    digest_prefix: bool,
    profile: bool = False,
) -> tuple[np.ndarray, list[FrameGeometry], str | None, DomainProfile | None]:
    alignment, ligand = _selections(universe, config)
    dt_ns = _frame_interval_ns(universe)
    reference_frame_centered: np.ndarray | None = None
    reference_ligand_relative: np.ndarray | None = None
    times: list[float] = []
    values: list[FrameGeometry] = []
    digest = hashlib.sha256() if digest_prefix else None
    max_diameter = 0.0
    min_box_length = float("inf")

    for frame_index, ts in enumerate(universe.trajectory):
        if frame_limit is not None and frame_index >= frame_limit:
            break
        time_ns = (frame_index + 1) * dt_ns
        box = np.asarray(ts.dimensions)
        frame_positions = np.asarray(alignment.positions)
        ligand_positions = np.asarray(ligand.positions)
        if digest is not None:
            digest.update(np.asarray(box, dtype="<f4").tobytes())
            digest.update(np.asarray(frame_positions, dtype="<f4").tobytes())
            digest.update(np.asarray(ligand_positions, dtype="<f4").tobytes())

        frame_centered, ligand_relative = protein_relative_coordinates(
            frame_positions, ligand_positions, box
        )
        if profile:
            spread = ligand_relative[:, None, :] - ligand_relative[None, :, :]
            max_diameter = max(
                max_diameter, float(np.sqrt((spread * spread).sum(axis=-1)).max())
            )
            min_box_length = min(min_box_length, float(np.asarray(box)[:3].min()))
        if frame_index == config.reference_frame:
            reference_frame_centered = frame_centered.copy()
            reference_ligand_relative = ligand_relative.copy()
            geometry = zero_geometry()
        else:
            if reference_frame_centered is None or reference_ligand_relative is None:
                raise InputValidationError("reference frame was not initialized")
            geometry = frame_relative_geometry(
                frame_centered,
                reference_frame_centered,
                ligand_relative,
                reference_ligand_relative,
            )
        times.append(time_ns)
        values.append(geometry)
    domain: DomainProfile | None = None
    if profile and times:
        domain = DomainProfile(
            max_ligand_diameter_angstrom=max_diameter,
            min_box_length_angstrom=min_box_length,
            ligand_diameter_box_fraction=max_diameter / min_box_length,
        )
    return (
        np.asarray(times, dtype=float),
        values,
        digest.hexdigest() if digest is not None else None,
        domain,
    )


def measure_prefix(
    topology: str | Path,
    trajectory: str | Path,
    config: PoseGateConfig,
    *,
    require_pre_outcome: bool = False,
) -> PrefixMeasurement:
    topology_path = Path(topology).resolve()
    trajectory_path = Path(trajectory).resolve()
    universe = _universe(topology_path, trajectory_path)
    alignment, ligand = _selections(universe, config)
    dt_ns = _frame_interval_ns(universe)
    available = len(universe.trajectory)
    available_end_ns = available * dt_ns
    if require_pre_outcome and available_end_ns >= config.outcome.window_start_ns:
        raise OutcomeEmbargoError(
            f"trajectory exposes {available_end_ns:.3f} ns; a genuine shadow "
            f"record requires less than {config.outcome.window_start_ns:g} ns"
        )
    required = _required_frames(config.checkpoint.time_ns, dt_ns)
    if available < required:
        raise CheckpointNotReady(
            f"checkpoint incomplete: frames={available}, required={required}, "
            f"end={available_end_ns:.3f} ns"
        )

    times, geometries, prefix_hash, domain = _measure_trace(
        universe, config, frame_limit=required, digest_prefix=True, profile=True
    )
    if len(times) != required:
        raise InputValidationError(
            f"prefix read was incomplete: used={len(times)}, required={required}"
        )
    window = (times > config.checkpoint.window_start_ns) & (
        times <= config.checkpoint.window_end_ns + TIME_TOLERANCE_NS
    )
    if not window.any():
        raise InputValidationError("configured checkpoint window contains no frames")
    if domain is None:
        raise InputValidationError("domain profile was not measured")
    limit = config.applicability
    if limit is not None and (
        domain.ligand_diameter_box_fraction > limit.max_ligand_diameter_box_fraction
    ):
        raise InputValidationError(
            "ligand extent is outside the declared domain of this policy: "
            f"diameter/box = {domain.ligand_diameter_box_fraction:.4f} exceeds "
            f"{limit.max_ligand_diameter_box_fraction:g}; minimum-image "
            "reconstruction is not unambiguous at this size"
        )

    rmsd = np.asarray([item.corrected_pose_rmsd_angstrom for item in geometries])
    centroid = np.asarray([item.centroid_displacement_angstrom for item in geometries])
    theta = np.asarray(
        [item.frame_relative_reorientation_radians for item in geometries]
    )
    internal = np.asarray([item.internal_deformation_angstrom for item in geometries])
    return PrefixMeasurement(
        topology_path=str(topology_path),
        trajectory_path=str(trajectory_path),
        alignment_selection=config.selections.alignment_atoms,
        ligand_selection=config.selections.ligand,
        alignment_atom_count=alignment.n_atoms,
        ligand_atom_count=ligand.n_atoms,
        trajectory_frames_available=available,
        prefix_frames_used=required,
        frame_interval_ns=dt_ns,
        max_available_time_ns=available_end_ns,
        checkpoint_ns=config.checkpoint.time_ns,
        window_start_ns=config.checkpoint.window_start_ns,
        window_end_ns=config.checkpoint.window_end_ns,
        window_frame_count=int(window.sum()),
        prefix_coordinates_sha256=str(prefix_hash),
        corrected_pose_rmsd_mean_angstrom=float(np.mean(rmsd[window])),
        corrected_centroid_displacement_mean_angstrom=float(np.mean(centroid[window])),
        frame_relative_reorientation_mean_degrees=float(
            np.degrees(np.mean(theta[window]))
        ),
        internal_deformation_mean_angstrom=float(np.mean(internal[window])),
        domain_profile=domain,
        series=PrefixSeries(
            time_ns=tuple(float(item) for item in times),
            corrected_pose_rmsd_angstrom=tuple(float(item) for item in rmsd),
            centroid_displacement_angstrom=tuple(float(item) for item in centroid),
            frame_relative_reorientation_degrees=tuple(
                float(item) for item in np.degrees(theta)
            ),
            internal_deformation_angstrom=tuple(float(item) for item in internal),
            in_feature_window=tuple(bool(item) for item in window),
        ),
    )


def measure_late_outcome(
    topology: str | Path,
    trajectory: str | Path,
    config: PoseGateConfig,
) -> LateOutcome:
    topology_path = Path(topology).resolve()
    trajectory_path = Path(trajectory).resolve()
    universe = _universe(topology_path, trajectory_path)
    dt_ns = _frame_interval_ns(universe)
    available = len(universe.trajectory)
    end_ns = available * dt_ns
    if end_ns < config.outcome.window_end_ns - TIME_TOLERANCE_NS:
        raise CheckpointNotReady(
            f"late outcome incomplete: end={end_ns:.3f} ns, "
            f"required={config.outcome.window_end_ns:g} ns"
        )
    times, geometries, _, _ = _measure_trace(
        universe, config, frame_limit=None, digest_prefix=False
    )
    late = (times > config.outcome.window_start_ns) & (
        times <= config.outcome.window_end_ns + TIME_TOLERANCE_NS
    )
    if not late.any():
        raise InputValidationError("configured late outcome window contains no frames")
    rmsd = np.asarray([item.corrected_pose_rmsd_angstrom for item in geometries])
    late_median = float(np.median(rmsd[late]))
    threshold = config.outcome.retained_rmsd_threshold_angstrom
    return LateOutcome(
        trajectory_frames=available,
        trajectory_end_ns=end_ns,
        late_window_start_ns=config.outcome.window_start_ns,
        late_window_end_ns=config.outcome.window_end_ns,
        late_window_frame_count=int(late.sum()),
        late_pose_rmsd_median_angstrom=late_median,
        retained_rmsd_threshold_angstrom=threshold,
        pose_retained=late_median < threshold,
    )
