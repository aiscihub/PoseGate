from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from posegate.config import config_from_mapping, config_from_record, load_config
from posegate.exceptions import ConfigurationError


ROOT = Path(__file__).resolve().parents[1]


def test_bundled_configs_are_strict_and_round_trip() -> None:
    five = load_config(ROOT / "configs/posegate_5ns_v1.yaml")
    twenty = load_config(ROOT / "configs/posegate_20ns_v1.yaml")

    assert five.checkpoint.time_ns == 5.0
    assert five.policy.feature_names == ("corrected_pose_rmsd_mean_angstrom",)
    assert twenty.checkpoint.time_ns == 20.0
    assert twenty.policy.feature_names == (
        "corrected_pose_rmsd_mean_angstrom",
        "corrected_centroid_displacement_mean_angstrom",
    )
    assert config_from_record(twenty.scientific_dict()).policy == twenty.policy


def test_missing_scientific_field_fails_loudly() -> None:
    raw = yaml.safe_load((ROOT / "configs/posegate_5ns_v1.yaml").read_text())
    broken = deepcopy(raw)
    del broken["selections"]["ligand"]
    with pytest.raises(ConfigurationError, match="selections.ligand"):
        config_from_mapping(broken, configuration_sha256="test")


def test_unknown_scientific_field_fails_loudly() -> None:
    raw = yaml.safe_load((ROOT / "configs/posegate_5ns_v1.yaml").read_text())
    broken = deepcopy(raw)
    broken["checkpoint"]["silent_fallback"] = True
    with pytest.raises(ConfigurationError, match="unknown fields"):
        config_from_mapping(broken, configuration_sha256="test")


def test_policy_vector_lengths_must_match() -> None:
    raw = yaml.safe_load((ROOT / "configs/posegate_20ns_v1.yaml").read_text())
    broken = deepcopy(raw)
    broken["policy"]["coefficients"] = [-1.0]
    with pytest.raises(ConfigurationError, match="vector lengths"):
        config_from_mapping(broken, configuration_sha256="test")
