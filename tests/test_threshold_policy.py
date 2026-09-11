"""One-sided threshold_rule policies."""

from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from posegate.config import config_from_mapping, config_from_record
from posegate.exceptions import ConfigurationError
from posegate.policy import evaluate_policy

ROOT = Path(__file__).resolve().parents[1]
THRESHOLD_A = 2.7634173197224974


def threshold_mapping() -> dict:
    raw = yaml.safe_load(
        (ROOT / "tests/fixtures/legacy/posegate_5ns_v1.yaml").read_text()
    )
    raw = deepcopy(raw)
    raw["policy_id"] = "example_5ns_rmsd_threshold"
    raw["policy"] = {
        "mode": "shadow",
        "type": "threshold_rule",
        "feature_names": ["corrected_pose_rmsd_mean_angstrom"],
        "stop_threshold_angstrom": THRESHOLD_A,
        "stop_direction": "greater_than",
        "tie_rule": "no_stop",
        "score_calibrated": False,
    }
    raw["applicability"] = {"max_ligand_diameter_box_fraction": 0.45}
    return raw


def test_threshold_rule_stops_above_and_defers_below() -> None:
    config = config_from_mapping(threshold_mapping(), configuration_sha256="test")

    above = evaluate_policy(
        config.policy, {"corrected_pose_rmsd_mean_angstrom": 4.168600883739357}
    )
    assert above.forecast == "POSE_NONRETAINED"
    assert above.recommendation == "CANDIDATE_FOR_PAUSE"
    assert above.signed_stop_margin_angstrom == pytest.approx(
        4.168600883739357 - THRESHOLD_A, abs=1e-15
    )
    assert above.uncalibrated_retention_score is None

    below = evaluate_policy(
        config.policy, {"corrected_pose_rmsd_mean_angstrom": 1.5692652937833238}
    )
    # Not meeting a one-sided stop condition is not a retention prediction.
    assert below.forecast == "DEFER"
    assert below.recommendation == "CONTINUE"
    assert below.signed_stop_margin_angstrom < 0


def test_threshold_exactly_equal_does_not_stop() -> None:
    config = config_from_mapping(threshold_mapping(), configuration_sha256="test")
    decision = evaluate_policy(
        config.policy, {"corrected_pose_rmsd_mean_angstrom": THRESHOLD_A}
    )
    assert decision.forecast == "DEFER"
    assert decision.signed_stop_margin_angstrom == 0.0


def test_threshold_policy_round_trips_through_a_record() -> None:
    config = config_from_mapping(threshold_mapping(), configuration_sha256="test")
    embedded = config.scientific_dict()
    assert "scaler_mean" not in embedded["policy"]
    assert embedded["applicability"] == {"max_ligand_diameter_box_fraction": 0.45}
    assert config_from_record(embedded).policy == config.policy


def test_logistic_fields_are_rejected_in_a_threshold_policy() -> None:
    raw = threshold_mapping()
    raw["policy"]["scaler_mean"] = [1.0]
    with pytest.raises(ConfigurationError, match="unknown fields"):
        config_from_mapping(raw, configuration_sha256="test")


def test_threshold_policy_cannot_claim_calibration() -> None:
    raw = threshold_mapping()
    raw["policy"]["score_calibrated"] = True
    with pytest.raises(ConfigurationError, match="cannot be calibrated"):
        config_from_mapping(raw, configuration_sha256="test")


def test_threshold_policy_requires_exactly_one_feature() -> None:
    raw = threshold_mapping()
    raw["policy"]["feature_names"] = [
        "corrected_pose_rmsd_mean_angstrom",
        "corrected_centroid_displacement_mean_angstrom",
    ]
    with pytest.raises(ConfigurationError, match="exactly one feature"):
        config_from_mapping(raw, configuration_sha256="test")


def test_tie_rule_and_direction_are_not_free_text() -> None:
    raw = threshold_mapping()
    raw["policy"]["stop_direction"] = "less_than"
    with pytest.raises(ConfigurationError, match="stop_direction"):
        config_from_mapping(raw, configuration_sha256="test")

    raw = threshold_mapping()
    raw["policy"]["tie_rule"] = "stop"
    with pytest.raises(ConfigurationError, match="tie_rule"):
        config_from_mapping(raw, configuration_sha256="test")
