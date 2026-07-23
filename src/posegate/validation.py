"""Close the loop between an immutable forecast and a completed trajectory."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Mapping

from .config import config_from_record
from .exceptions import ImmutableRecordError
from .records import sha256_file, sha256_json
from .trajectory import measure_late_outcome, measure_prefix


def validate_shadow_record(
    record: Mapping[str, Any],
    *,
    topology: str | Path,
    trajectory: str | Path,
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

    outcome = measure_late_outcome(topology_path, trajectory_path, config)
    predicted_retained = record.get("forecast") == "POSE_RETAINED"
    decision_correct = predicted_retained == outcome.pose_retained
    return {
        "schema_version": "1.0",
        "record_id": record["record_id"],
        "policy_id": record["policy_id"],
        "topology_sha256": topology_hash,
        "topology_matches_record": topology_matches,
        "prefix_coordinates_match_record": prefix_matches,
        "early_rmsd_cross_check_delta_angstrom": early_rmsd_delta,
        "early_centroid_cross_check_delta_angstrom": early_centroid_delta,
        "forecast": record["forecast"],
        "uncalibrated_retention_score": record["uncalibrated_retention_score"],
        "outcome": outcome.as_dict(),
        "decision_was_correct": decision_correct,
        "scope": "FORECASTS_THIS_TRAJECTORY_ONLY",
    }
