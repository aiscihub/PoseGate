import json
import os
from pathlib import Path

import pytest

from posegate.config import load_config
from posegate.policy import evaluate_policy
from posegate.trajectory import measure_prefix


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.integration
def test_real_pocket12_20ns_matches_sealed_legacy_capture() -> None:
    run_dir_value = os.environ.get("POSEGATE_POCKET12_RUN_DIR")
    capture_value = os.environ.get("POSEGATE_POCKET12_CAPTURE")
    if not run_dir_value or not capture_value:
        pytest.skip("set POSEGATE_POCKET12_RUN_DIR and POSEGATE_POCKET12_CAPTURE")

    run_dir = Path(run_dir_value)
    topology = [
        path
        for path in run_dir.glob("*_explicit_initial_frame.pdb")
        if not path.name.endswith("_stripped_initial_frame.pdb")
    ]
    trajectories = list(run_dir.glob("*_explicit_trajectory.dcd"))
    assert len(topology) == 1
    assert len(trajectories) == 1

    legacy = json.loads(Path(capture_value).read_text())
    config = load_config(ROOT / "configs/posegate_20ns_v1.yaml")
    measured = measure_prefix(topology[0], trajectories[0], config)
    decision = evaluate_policy(config.policy, measured.policy_measurements())

    assert measured.corrected_pose_rmsd_mean_angstrom == pytest.approx(
        legacy["corrected_pose_rmsd_mean_18_20_A"], abs=1e-12
    )
    assert measured.corrected_centroid_displacement_mean_angstrom == pytest.approx(
        legacy["corrected_centroid_displacement_mean_18_20_A"], abs=1e-12
    )
    assert measured.prefix_coordinates_sha256 == legacy["prefix_coordinates_sha256"]
    assert decision.uncalibrated_retention_score == pytest.approx(
        legacy["uncalibrated_retention_score"], abs=1e-14
    )
    assert decision.recommendation.lower() == legacy["stop_continue_decision"]
