"""HTML rendering of sealed records."""

from pathlib import Path

import pytest

from posegate.config import config_from_mapping, load_config
from posegate.policy import evaluate_policy
from posegate.records import build_run_identity, build_shadow_record
from posegate.report import render_record_html, write_record_report

from test_record_schema import synthetic_measurement
from test_threshold_policy import threshold_mapping


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture()
def logistic_record(tmp_path: Path) -> dict:
    topology = tmp_path / "system.pdb"
    trajectory = tmp_path / "trajectory.dcd"
    topology.write_bytes(b"REMARK synthetic\n")
    trajectory.write_bytes(b"synthetic")
    config = load_config(ROOT / "configs/posegate_5ns_v1.yaml")
    measurement = synthetic_measurement(topology, trajectory)
    decision = evaluate_policy(config.policy, measurement.policy_measurements())
    return build_shadow_record(
        run_id="CIMG_06197|Milbemycin|pocket15|replica_1",
        config=config,
        measurement=measurement,
        decision=decision,
        run_identity=build_run_identity(
            protein="CIMG_06197",
            ligand="Milbemycin",
            pocket="pocket15",
            replica="replica_1",
        ),
    )


def test_report_is_self_contained(logistic_record: dict) -> None:
    html = render_record_html(logistic_record)
    # A sealed record must stay readable from an archive with no network.
    for forbidden in ("http://", "https://", "<script", "src="):
        assert forbidden not in html
    assert "<title>" in html and "<style>" in html


def test_report_shows_identity_window_and_hashes(logistic_record: dict) -> None:
    html = render_record_html(logistic_record)
    assert "CIMG_06197 : pocket15" in html
    assert "Milbemycin" in html
    assert "Outcome-blinded shadow" in html
    assert logistic_record["trajectory_prefix_coordinates_sha256"][:16] in html
    assert "Uncalibrated score" in html
    assert "ranking score, not a probability" in html


def test_report_plots_every_prefix_frame(logistic_record: dict) -> None:
    html = render_record_html(logistic_record)
    points = html.split('points="')[1].split('"')[0].split()
    assert len(points) == len(logistic_record["measurement"]["series"]["time_ns"])


def test_threshold_record_reports_margin_not_score(tmp_path: Path) -> None:
    topology = tmp_path / "system.pdb"
    trajectory = tmp_path / "trajectory.dcd"
    topology.write_bytes(b"REMARK synthetic\n")
    trajectory.write_bytes(b"synthetic")
    config = config_from_mapping(threshold_mapping(), configuration_sha256="test")
    measurement = synthetic_measurement(topology, trajectory)
    decision = evaluate_policy(config.policy, measurement.policy_measurements())
    record = build_shadow_record(
        run_id="example", config=config, measurement=measurement, decision=decision
    )
    html = render_record_html(record)
    assert "Stop threshold" in html
    assert "Signed margin" in html
    assert "Uncalibrated score" not in html
    assert "not a prediction that the pose is retained" in html
    assert "2.7634" in html


def test_report_renders_outcome_and_tamper_row(
    logistic_record: dict, tmp_path: Path
) -> None:
    validation = {
        "audit_category": "false_continue",
        "decision_was_correct": False,
        "outcome": {
            "late_pose_rmsd_median_angstrom": 4.892,
            "retained_rmsd_threshold_angstrom": 3.0,
            "late_window_start_ns": 70.0,
            "late_window_end_ns": 100.0,
            "pose_retained": False,
        },
        "record_integrity": {
            "sealed_forecast_still_follows_from_sealed_measurement": True,
            "record_file_read_only": True,
        },
    }
    path = write_record_report(tmp_path / "run.html", logistic_record, validation)
    html = path.read_text()
    assert "false_continue" in html
    assert "Non-retained" in html
    assert "Original forecast modified" in html
    assert "read-only" in html


def test_schema_1_0_record_without_series_still_renders(logistic_record: dict) -> None:
    legacy = dict(logistic_record)
    legacy["schema_version"] = "1.0"
    legacy["quality_control"] = "PASS"
    measurement = dict(legacy["measurement"])
    measurement.pop("series")
    legacy["measurement"] = measurement
    html = render_record_html(legacy)
    assert "predates per-frame series capture" in html
    assert "Schema 1.0 record" in html


def test_report_flags_a_feature_window_that_moved(logistic_record: dict) -> None:
    validation = {
        "audit_category": "successful_stop",
        "decision_was_correct": True,
        "outcome": {
            "late_pose_rmsd_median_angstrom": 8.427,
            "retained_rmsd_threshold_angstrom": 3.0,
            "late_window_start_ns": 70.0,
            "late_window_end_ns": 100.0,
            "pose_retained": False,
        },
        "feature_window_cross_check": {
            "frames_in_record": 10,
            "frames_on_revalidation": 11,
            "frame_interval_ns_in_record": 0.19999999612625977,
            "frame_interval_ns_on_revalidation": 0.20000000059628115,
            "matches_record": False,
        },
        "record_integrity": {
            "sealed_forecast_still_follows_from_sealed_measurement": True
        },
    }
    html = render_record_html(logistic_record, validation)
    assert "Feature window on revalidation" in html
    assert "10 &rarr; 11 frames" in html
