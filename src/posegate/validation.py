"""Close the loop between an immutable forecast and a completed trajectory."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .config import config_from_record
from .exceptions import ImmutableRecordError
from .policy import evaluate_policy
from .records import sha256_file, sha256_json
from .trajectory import measure_late_outcome, measure_prefix


#: Forecast x outcome, in the vocabulary of the frozen policy audit.
AUDIT_CATEGORIES = {
    ("POSE_NONRETAINED", False): "successful_stop",
    ("POSE_NONRETAINED", True): "false_stop",
    ("POSE_RETAINED", True): "safe_continue",
    ("POSE_RETAINED", False): "false_continue",
    ("DEFER", True): "no_action_late_retained",
    ("DEFER", False): "no_action_late_lost",
}


def audit_category(forecast: str, pose_retained: bool) -> str:
    try:
        return AUDIT_CATEGORIES[(forecast, bool(pose_retained))]
    except KeyError as exc:
        raise ImmutableRecordError(f"unsupported forecast value: {forecast}") from exc


def _recompute_forecast(
    record: Mapping[str, Any], config: Any, stored_measurement: Mapping[str, Any]
) -> dict[str, Any]:
    """Re-apply the sealed policy to the sealed measurement.

    This is the tamper check behind the report's "original forecast modified"
    row: it catches a record whose stored forecast, score, or margin no longer
    follows from its own stored measurement under its own embedded policy.
    """
    features = {
        name: float(stored_measurement[name])
        for name in config.policy.feature_names
        if name in stored_measurement
    }
    missing = sorted(set(config.policy.feature_names).difference(features))
    if missing:
        raise ImmutableRecordError(
            f"record measurement is missing policy features: {missing}"
        )
    recomputed = evaluate_policy(config.policy, features)
    fields = {
        "forecast": recomputed.forecast,
        "recommendation": recomputed.recommendation,
        "uncalibrated_retention_score": recomputed.uncalibrated_retention_score,
        "signed_stop_margin_angstrom": recomputed.signed_stop_margin_angstrom,
    }
    mismatches = sorted(
        name
        for name, value in fields.items()
        if name in record and record.get(name) != value
    )
    return {
        "forecast_recomputed_from_sealed_measurement": recomputed.forecast,
        "sealed_forecast_still_follows_from_sealed_measurement": not mismatches,
        "recomputation_mismatched_fields": mismatches,
    }


def validate_shadow_record(
    record: Mapping[str, Any],
    *,
    topology: str | Path,
    trajectory: str | Path,
    record_path: str | Path | None = None,
) -> dict[str, Any]:
    if record.get("late_outcome_accessed") is not False:
        raise ImmutableRecordError(
            "record does not attest that the late outcome was inaccessible"
        )
    if record.get("scope") != "FORECASTS_THIS_TRAJECTORY_ONLY":
        raise ImmutableRecordError("record has an unsupported forecast scope")
    raw_config = record.get("configuration")
    if not isinstance(raw_config, Mapping):
        raise ImmutableRecordError("record does not embed its scientific configuration")
    observed_config_hash = sha256_json(dict(raw_config))
    if observed_config_hash != record.get("configuration_scientific_sha256"):
        raise ImmutableRecordError(
            "embedded scientific configuration does not match its sealed hash"
        )
    config = config_from_record(raw_config)

    topology_path = Path(topology).resolve()
    trajectory_path = Path(trajectory).resolve()
    topology_hash = sha256_file(topology_path)
    topology_matches = topology_hash == record.get("topology_sha256")
    if not topology_matches:
        raise ImmutableRecordError(
            "validation topology does not match the sealed forecast"
        )

    prefix = measure_prefix(
        topology_path, trajectory_path, config, require_pre_outcome=False
    )
    stored_measurement = record.get("measurement")
    if not isinstance(stored_measurement, Mapping):
        raise ImmutableRecordError("record measurement is missing")
    stored_prefix_hash = record.get("trajectory_prefix_coordinates_sha256")
    prefix_matches = prefix.prefix_coordinates_sha256 == stored_prefix_hash
    if not prefix_matches:
        raise ImmutableRecordError(
            "validation trajectory prefix does not match the sealed forecast"
        )
    early_rmsd_delta = abs(
        prefix.corrected_pose_rmsd_mean_angstrom
        - float(stored_measurement["corrected_pose_rmsd_mean_angstrom"])
    )
    early_centroid_delta = abs(
        prefix.corrected_centroid_displacement_mean_angstrom
        - float(stored_measurement["corrected_centroid_displacement_mean_angstrom"])
    )
    # Legacy records use DCD-derived times; current records use nominal indices.
    # Preserve the sealed contract and expose the window used on revalidation.
    window_cross_check = {
        "frames_in_record": stored_measurement.get("window_frame_count"),
        "frames_on_revalidation": prefix.window_frame_count,
        "frame_interval_ns_in_record": stored_measurement.get("frame_interval_ns"),
        "frame_interval_ns_on_revalidation": prefix.frame_interval_ns,
    }
    window_cross_check["matches_record"] = (
        window_cross_check["frames_in_record"]
        == window_cross_check["frames_on_revalidation"]
    )

    integrity = _recompute_forecast(record, config, stored_measurement)
    if record_path is not None:
        sealed = Path(record_path)
        integrity["record_file_sha256"] = sha256_file(sealed)
        integrity["record_file_read_only"] = not bool(sealed.stat().st_mode & 0o222)

    outcome = measure_late_outcome(topology_path, trajectory_path, config)
    forecast = str(record["forecast"])
    category = audit_category(forecast, outcome.pose_retained)
    # A deferral makes no retention claim, so it is neither right nor wrong.
    decision_correct: bool | None = None
    if forecast != "DEFER":
        decision_correct = (forecast == "POSE_RETAINED") == outcome.pose_retained
    return {
        "schema_version": "1.1",
        "record_id": record["record_id"],
        "policy_id": record["policy_id"],
        "policy_type": record.get("policy_type", "standardized_logistic"),
        "run_identity": record.get("run_identity"),
        "topology_sha256": topology_hash,
        "topology_matches_record": topology_matches,
        "prefix_coordinates_match_record": prefix_matches,
        "early_rmsd_cross_check_delta_angstrom": early_rmsd_delta,
        "early_centroid_cross_check_delta_angstrom": early_centroid_delta,
        "feature_window_cross_check": window_cross_check,
        "forecast": forecast,
        "recommendation": record.get("recommendation"),
        "uncalibrated_retention_score": record["uncalibrated_retention_score"],
        "signed_stop_margin_angstrom": record.get("signed_stop_margin_angstrom"),
        "record_integrity": integrity,
        "outcome": outcome.as_dict(),
        "audit_category": category,
        "decision_was_correct": decision_correct,
        "scope": "FORECASTS_THIS_TRAJECTORY_ONLY",
    }
