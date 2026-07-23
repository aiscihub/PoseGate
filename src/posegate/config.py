"""Strict loading of the versioned PoseGate-MD scientific contract."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import yaml

from .exceptions import ConfigurationError


SUPPORTED_SCHEMA_VERSION = "1.0"
SUPPORTED_PBC_METHOD = "per_atom_minimum_image_to_alignment_centroid_v1"
SUPPORTED_POLICY_TYPE = "standardized_logistic"
SUPPORTED_FEATURES = {
    "corrected_pose_rmsd_mean_angstrom",
    "corrected_centroid_displacement_mean_angstrom",
}


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ConfigurationError(f"{path} must be a mapping")
    return value


def _require(mapping: Mapping[str, Any], key: str, path: str) -> Any:
    if key not in mapping:
        raise ConfigurationError(f"missing required field {path}.{key}")
    return mapping[key]


def _reject_unknown(mapping: Mapping[str, Any], allowed: set[str], path: str) -> None:
    unknown = sorted(set(mapping).difference(allowed))
    if unknown:
        raise ConfigurationError(f"unknown fields in {path}: {unknown}")


def _finite_float(value: Any, path: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ConfigurationError(f"{path} must be numeric") from exc
    if not (-float("inf") < number < float("inf")):
        raise ConfigurationError(f"{path} must be finite")
    return number


def _boolean(value: Any, path: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigurationError(f"{path} must be a YAML boolean")
    return value


def _float_tuple(value: Any, path: str) -> tuple[float, ...]:
    if not isinstance(value, (list, tuple)) or not value:
        raise ConfigurationError(f"{path} must be a non-empty list")
    return tuple(
        _finite_float(item, f"{path}[{index}]") for index, item in enumerate(value)
    )


@dataclass(frozen=True)
class Selections:
    protein: str
    alignment_atoms: str
    ligand: str


@dataclass(frozen=True)
class PBC:
    method: str
    require_box_vectors: bool


@dataclass(frozen=True)
class Checkpoint:
    time_ns: float
    window_start_ns: float
    window_end_ns: float


@dataclass(frozen=True)
class Measurements:
    corrected_rmsd: bool
    centroid_displacement: bool
    frame_relative_reorientation: bool
    internal_deformation: bool


@dataclass(frozen=True)
class Policy:
    mode: str
    type: str
    feature_names: tuple[str, ...]
    scaler_mean: tuple[float, ...]
    scaler_scale: tuple[float, ...]
    coefficients: tuple[float, ...]
    intercept: float
    decision_threshold: float
    score_calibrated: bool


@dataclass(frozen=True)
class Outcome:
    window_start_ns: float
    window_end_ns: float
    retained_rmsd_threshold_angstrom: float


@dataclass(frozen=True)
class PoseGateConfig:
    schema_version: str
    policy_id: str
    selections: Selections
    reference_frame: int
    pbc: PBC
    checkpoint: Checkpoint
    measurements: Measurements
    policy: Policy
    outcome: Outcome
    provenance: Mapping[str, Any]
    configuration_sha256: str

    def scientific_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value.pop("configuration_sha256")
        value["reference"] = {"frame": value.pop("reference_frame")}
        return value


def config_from_mapping(
    raw: Mapping[str, Any], *, configuration_sha256: str
) -> PoseGateConfig:
    _reject_unknown(
        raw,
        {
            "schema_version",
            "policy_id",
            "selections",
            "reference",
            "pbc",
            "checkpoint",
            "measurements",
            "policy",
            "outcome",
            "provenance",
        },
        "config",
    )
    schema_version = str(_require(raw, "schema_version", "config"))
    if schema_version != SUPPORTED_SCHEMA_VERSION:
        raise ConfigurationError(
            f"schema_version must be {SUPPORTED_SCHEMA_VERSION}, got {schema_version}"
        )
    policy_id = str(_require(raw, "policy_id", "config")).strip()
    if not policy_id:
        raise ConfigurationError("policy_id must not be empty")

    selections_raw = _mapping(_require(raw, "selections", "config"), "selections")
    _reject_unknown(
        selections_raw, {"protein", "alignment_atoms", "ligand"}, "selections"
    )
    selections = Selections(
        **{
            key: str(_require(selections_raw, key, "selections")).strip()
            for key in ("protein", "alignment_atoms", "ligand")
        }
    )
    if not all(asdict(selections).values()):
        raise ConfigurationError("selection expressions must not be empty")

    reference_raw = _mapping(_require(raw, "reference", "config"), "reference")
    _reject_unknown(reference_raw, {"frame"}, "reference")
    reference_frame = int(_require(reference_raw, "frame", "reference"))
    if reference_frame != 0:
        raise ConfigurationError("version 0.1 requires reference.frame = 0")

    pbc_raw = _mapping(_require(raw, "pbc", "config"), "pbc")
    _reject_unknown(pbc_raw, {"method", "require_box_vectors"}, "pbc")
    pbc = PBC(
        method=str(_require(pbc_raw, "method", "pbc")),
        require_box_vectors=_boolean(
            _require(pbc_raw, "require_box_vectors", "pbc"),
            "pbc.require_box_vectors",
        ),
    )
    if pbc.method != SUPPORTED_PBC_METHOD:
        raise ConfigurationError(f"unsupported pbc.method: {pbc.method}")
    if not pbc.require_box_vectors:
        raise ConfigurationError("version 0.1 requires periodic box vectors")

    checkpoint_raw = _mapping(_require(raw, "checkpoint", "config"), "checkpoint")
    _reject_unknown(
        checkpoint_raw,
        {"time_ns", "window_start_ns", "window_end_ns"},
        "checkpoint",
    )
    checkpoint = Checkpoint(
        time_ns=_finite_float(
            _require(checkpoint_raw, "time_ns", "checkpoint"),
            "checkpoint.time_ns",
        ),
        window_start_ns=_finite_float(
            _require(checkpoint_raw, "window_start_ns", "checkpoint"),
            "checkpoint.window_start_ns",
        ),
        window_end_ns=_finite_float(
            _require(checkpoint_raw, "window_end_ns", "checkpoint"),
            "checkpoint.window_end_ns",
        ),
    )
    if not (
        0.0
        <= checkpoint.window_start_ns
        < checkpoint.window_end_ns
        == checkpoint.time_ns
    ):
        raise ConfigurationError(
            "checkpoint must satisfy 0 <= window_start_ns < " "window_end_ns == time_ns"
        )

    measurements_raw = _mapping(_require(raw, "measurements", "config"), "measurements")
    measurement_keys = {
        "corrected_rmsd",
        "centroid_displacement",
        "frame_relative_reorientation",
        "internal_deformation",
    }
    _reject_unknown(measurements_raw, measurement_keys, "measurements")
    measurements = Measurements(
        **{
            key: _boolean(
                _require(measurements_raw, key, "measurements"),
                f"measurements.{key}",
            )
            for key in measurement_keys
        }
    )
    if not measurements.corrected_rmsd:
        raise ConfigurationError("corrected_rmsd is required in version 0.1")

    policy_raw = _mapping(_require(raw, "policy", "config"), "policy")
    policy_keys = {
        "mode",
        "type",
        "feature_names",
        "scaler_mean",
        "scaler_scale",
        "coefficients",
        "intercept",
        "decision_threshold",
        "score_calibrated",
    }
    _reject_unknown(policy_raw, policy_keys, "policy")
    feature_value = _require(policy_raw, "feature_names", "policy")
    if not isinstance(feature_value, (list, tuple)) or not feature_value:
        raise ConfigurationError("policy.feature_names must be a non-empty list")
    feature_names = tuple(str(item) for item in feature_value)
    unsupported = sorted(set(feature_names).difference(SUPPORTED_FEATURES))
    if unsupported:
        raise ConfigurationError(f"unsupported policy features: {unsupported}")
    policy = Policy(
        mode=str(_require(policy_raw, "mode", "policy")),
        type=str(_require(policy_raw, "type", "policy")),
        feature_names=feature_names,
        scaler_mean=_float_tuple(
            _require(policy_raw, "scaler_mean", "policy"), "policy.scaler_mean"
        ),
        scaler_scale=_float_tuple(
            _require(policy_raw, "scaler_scale", "policy"), "policy.scaler_scale"
        ),
        coefficients=_float_tuple(
            _require(policy_raw, "coefficients", "policy"),
            "policy.coefficients",
        ),
        intercept=_finite_float(
            _require(policy_raw, "intercept", "policy"), "policy.intercept"
        ),
        decision_threshold=_finite_float(
            _require(policy_raw, "decision_threshold", "policy"),
            "policy.decision_threshold",
        ),
        score_calibrated=_boolean(
            _require(policy_raw, "score_calibrated", "policy"),
            "policy.score_calibrated",
        ),
    )
    if policy.mode != "shadow":
        raise ConfigurationError("version 0.1 supports only policy.mode = shadow")
    if policy.type != SUPPORTED_POLICY_TYPE:
        raise ConfigurationError(f"unsupported policy.type: {policy.type}")
    lengths = {
        len(policy.feature_names),
        len(policy.scaler_mean),
        len(policy.scaler_scale),
        len(policy.coefficients),
    }
    if len(lengths) != 1:
        raise ConfigurationError("policy vector lengths must match")
    if any(scale <= 0 for scale in policy.scaler_scale):
        raise ConfigurationError("policy.scaler_scale values must be positive")
    if not 0.0 < policy.decision_threshold < 1.0:
        raise ConfigurationError("policy.decision_threshold must be between 0 and 1")

    outcome_raw = _mapping(_require(raw, "outcome", "config"), "outcome")
    _reject_unknown(
        outcome_raw,
        {
            "window_start_ns",
            "window_end_ns",
            "retained_rmsd_threshold_angstrom",
        },
        "outcome",
    )
    outcome = Outcome(
        window_start_ns=_finite_float(
            _require(outcome_raw, "window_start_ns", "outcome"),
            "outcome.window_start_ns",
        ),
        window_end_ns=_finite_float(
            _require(outcome_raw, "window_end_ns", "outcome"),
            "outcome.window_end_ns",
        ),
        retained_rmsd_threshold_angstrom=_finite_float(
            _require(outcome_raw, "retained_rmsd_threshold_angstrom", "outcome"),
            "outcome.retained_rmsd_threshold_angstrom",
        ),
    )
    if not (checkpoint.time_ns < outcome.window_start_ns < outcome.window_end_ns):
        raise ConfigurationError(
            "outcome window must start after the checkpoint and have positive width"
        )
    if outcome.retained_rmsd_threshold_angstrom <= 0:
        raise ConfigurationError("outcome RMSD threshold must be positive")

    provenance = _mapping(_require(raw, "provenance", "config"), "provenance")
    return PoseGateConfig(
        schema_version=schema_version,
        policy_id=policy_id,
        selections=selections,
        reference_frame=reference_frame,
        pbc=pbc,
        checkpoint=checkpoint,
        measurements=measurements,
        policy=policy,
        outcome=outcome,
        provenance=dict(provenance),
        configuration_sha256=configuration_sha256,
    )


def load_config(path: str | Path) -> PoseGateConfig:
    config_path = Path(path)
    payload = config_path.read_bytes()
    try:
        raw = yaml.safe_load(payload)
    except yaml.YAMLError as exc:
        raise ConfigurationError(f"invalid YAML in {config_path}: {exc}") from exc
    mapping = _mapping(raw, "config")
    return config_from_mapping(
        mapping, configuration_sha256=hashlib.sha256(payload).hexdigest()
    )


def config_from_record(value: Mapping[str, Any]) -> PoseGateConfig:
    canonical = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return config_from_mapping(
        value, configuration_sha256=hashlib.sha256(canonical).hexdigest()
    )
