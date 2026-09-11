from pathlib import Path

import pytest

from posegate.config import load_config
from posegate.policy import evaluate_policy

ROOT = Path(__file__).resolve().parents[1]


def test_5ns_portable_policy_reproduces_frozen_shadow_score() -> None:
    config = load_config(ROOT / "tests/fixtures/legacy/posegate_5ns_v1.yaml")
    decision = evaluate_policy(
        config.policy,
        {
            "corrected_pose_rmsd_mean_angstrom": 1.5692652937833238,
            "corrected_centroid_displacement_mean_angstrom": 1.3013851900927742,
        },
    )
    assert decision.uncalibrated_retention_score == pytest.approx(
        0.7172454274964617, abs=1e-14
    )
    assert decision.forecast == "POSE_RETAINED"
    assert decision.recommendation == "CONTINUE"


def test_20ns_portable_policy_reproduces_frozen_shadow_score() -> None:
    config = load_config(ROOT / "tests/fixtures/legacy/posegate_20ns_v1.yaml")
    decision = evaluate_policy(
        config.policy,
        {
            "corrected_pose_rmsd_mean_angstrom": 2.051554013411634,
            "corrected_centroid_displacement_mean_angstrom": 1.8469326657070295,
        },
    )
    assert decision.uncalibrated_retention_score == pytest.approx(
        0.7199208871876361, abs=1e-14
    )
    assert decision.forecast == "POSE_RETAINED"
