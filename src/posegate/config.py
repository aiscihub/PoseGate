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
WHOLE_LIGAND_PBC_METHOD = "whole_ligand_minimum_image_to_alignment_centroid_v2"
SUPPORTED_PBC_METHODS = {SUPPORTED_PBC_METHOD, WHOLE_LIGAND_PBC_METHOD}
SUPPORTED_POLICY_TYPES = ("standardized_logistic", "threshold_rule")
SUPPORTED_STOP_DIRECTION = "greater_than"
SUPPORTED_TIE_RULE = "no_stop"
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
class Timing:
    """Nominal saving interval; floating DCD timestamps do not select frames."""

    frame_interval_ns: float


@dataclass(frozen=True)
class Measurements:
    corrected_rmsd: bool
    centroid_displacement: bool
    frame_relative_reorientation: bool
    internal_deformation: bool


@dataclass(frozen=True)
class Policy:
    """A frozen forecasting rule.

    ``standardized_logistic`` uses the scaler/coefficient fields and emits an
    uncalibrated retention score. ``threshold_rule`` uses the one-sided
    ``stop_threshold_angstrom`` fields and emits no score at all. Fields that do
    not apply to the declared type stay ``None`` and are pruned from
    :meth:`PoseGateConfig.scientific_dict`, so a policy of one type hashes
    exactly as it did before the other type existed.
    """

    mode: str
    type: str
    feature_names: tuple[str, ...]
    score_calibrated: bool
    scaler_mean: tuple[float, ...] | None = None
    scaler_scale: tuple[float, ...] | None = None
    coefficients: tuple[float, ...] | None = None
    intercept: float | None = None
    decision_threshold: float | None = None
    stop_threshold_angstrom: float | None = None
    stop_direction: str | None = None
    tie_rule: str | None = None


@dataclass(frozen=True)
class Applicability:
    """Optional domain-profile limits checked before a policy is applied."""

    max_ligand_diameter_box_fraction: float


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
    applicability: Applicability | None = None
    timing: Timing | None = None

    def scientific_dict(self) -> dict[str, Any]:
        value = asdict(self)
        value.pop("configuration_sha256")
        value["reference"] = {"frame": value.pop("reference_frame")}
        value["policy"] = {
            key: item for key, item in value["policy"].items() if item is not None
        }
        if value.get("applicability") is None:
            value.pop("applicability", None)
        if value.get("timing") is None:
            value.pop("timing", None)
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
            "applicability",
            "timing",
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
        raise ConfigurationError("reference.frame must be 0")

    pbc_raw = _mapping(_require(raw, "pbc", "config"), "pbc")
    _reject_unknown(pbc_raw, {"method", "require_box_vectors"}, "pbc")
    pbc = PBC(
        method=str(_require(pbc_raw, "method", "pbc")),
        require_box_vectors=_boolean(
            _require(pbc_raw, "require_box_vectors", "pbc"),
            "pbc.require_box_vectors",
        ),
    )
    if pbc.method not in SUPPORTED_PBC_METHODS:
        raise ConfigurationError(f"unsupported pbc.method: {pbc.method}")
    if not pbc.require_box_vectors:
        raise ConfigurationError("periodic box vectors are required")

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
        raise ConfigurationError("corrected_rmsd is required")

    policy_raw = _mapping(_require(raw, "policy", "config"), "policy")
    shared_keys = {"mode", "type", "feature_names", "score_calibrated"}
    logistic_keys = {
        "scaler_mean",
        "scaler_scale",
        "coefficients",
        "intercept",
        "decision_threshold",
    }
    threshold_keys = {"stop_threshold_angstrom", "stop_direction", "tie_rule"}
    policy_type = str(_require(policy_raw, "type", "policy"))
    if policy_type not in SUPPORTED_POLICY_TYPES:
        raise ConfigurationError(f"unsupported policy.type: {policy_type}")
    type_keys = (
        logistic_keys
        if policy_type == "standardized_logistic"
        else threshold_keys
    )
    _reject_unknown(policy_raw, shared_keys | type_keys, f"policy ({policy_type})")

    feature_value = _require(policy_raw, "feature_names", "policy")
    if not isinstance(feature_value, (list, tuple)) or not feature_value:
        raise ConfigurationError("policy.feature_names must be a non-empty list")
    feature_names = tuple(str(item) for item in feature_value)
    unsupported = sorted(set(feature_names).difference(SUPPORTED_FEATURES))
    if unsupported:
        raise ConfigurationError(f"unsupported policy features: {unsupported}")

    common = {
        "mode": str(_require(policy_raw, "mode", "policy")),
        "type": policy_type,
        "feature_names": feature_names,
        "score_calibrated": _boolean(
            _require(policy_raw, "score_calibrated", "policy"),
            "policy.score_calibrated",
        ),
    }
    if policy_type == "standardized_logistic":
        policy = Policy(
            **common,
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
        )
    else:
        policy = Policy(
            **common,
            stop_threshold_angstrom=_finite_float(
                _require(policy_raw, "stop_threshold_angstrom", "policy"),
                "policy.stop_threshold_angstrom",
            ),
            stop_direction=str(_require(policy_raw, "stop_direction", "policy")),
            tie_rule=str(_require(policy_raw, "tie_rule", "policy")),
        )

    if policy.mode != "shadow":
        raise ConfigurationError("only policy.mode = shadow is supported")

    if policy.type == "standardized_logistic":
        lengths = {
            len(policy.feature_names),
            len(policy.scaler_mean or ()),
            len(policy.scaler_scale or ()),
            len(policy.coefficients or ()),
        }
        if len(lengths) != 1:
            raise ConfigurationError("policy vector lengths must match")
        if any(scale <= 0 for scale in policy.scaler_scale or ()):
            raise ConfigurationError("policy.scaler_scale values must be positive")
        if not 0.0 < (policy.decision_threshold or 0.0) < 1.0:
            raise ConfigurationError(
                "policy.decision_threshold must be between 0 and 1"
            )
    else:
        if len(policy.feature_names) != 1:
            raise ConfigurationError(
                "a threshold_rule policy must declare exactly one feature"
            )
        if (policy.stop_threshold_angstrom or 0.0) <= 0:
            raise ConfigurationError("policy.stop_threshold_angstrom must be positive")
        if policy.stop_direction != SUPPORTED_STOP_DIRECTION:
            raise ConfigurationError(
                f"policy.stop_direction must be {SUPPORTED_STOP_DIRECTION}"
            )
        if policy.tie_rule != SUPPORTED_TIE_RULE:
            raise ConfigurationError(f"policy.tie_rule must be {SUPPORTED_TIE_RULE}")
        if policy.score_calibrated:
            raise ConfigurationError(
                "a threshold_rule policy emits no score and cannot be calibrated"
            )

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

    applicability: Applicability | None = None
    if raw.get("applicability") is not None:
        applicability_raw = _mapping(raw["applicability"], "applicability")
        _reject_unknown(
            applicability_raw, {"max_ligand_diameter_box_fraction"}, "applicability"
        )
        fraction = _finite_float(
            _require(
                applicability_raw,
                "max_ligand_diameter_box_fraction",
                "applicability",
            ),
            "applicability.max_ligand_diameter_box_fraction",
        )
        if not 0.0 < fraction < 1.0:
            raise ConfigurationError(
                "applicability.max_ligand_diameter_box_fraction must be between 0 and 1"
            )
        applicability = Applicability(max_ligand_diameter_box_fraction=fraction)

    timing = None
    if raw.get("timing") is not None:
        timing_raw = _mapping(raw["timing"], "timing")
        _reject_unknown(timing_raw, {"frame_interval_ns"}, "timing")
        interval = _finite_float(_require(timing_raw, "frame_interval_ns", "timing"),
                                 "timing.frame_interval_ns")
        if interval <= 0:
            raise ConfigurationError("timing.frame_interval_ns must be positive")
        timing = Timing(frame_interval_ns=interval)
        for boundary in (checkpoint.window_start_ns, checkpoint.window_end_ns,
                         outcome.window_start_ns, outcome.window_end_ns):
            if abs(round(boundary / interval) * interval - boundary) > 1e-8:
                raise ConfigurationError("window endpoints must lie on the nominal frame grid")
    if pbc.method == WHOLE_LIGAND_PBC_METHOD:
        if timing is None:
            raise ConfigurationError("whole-ligand v2 requires explicit nominal timing")
        if applicability is None or applicability.max_ligand_diameter_box_fraction >= 0.5:
            raise ConfigurationError("whole-ligand v2 requires a ligand diameter limit below 0.5")
    elif timing is not None:
        raise ConfigurationError("legacy v1 timing must remain unchanged")

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
        applicability=applicability,
        timing=timing,
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
