"""Schema 1.1 additions: identity, named QC checks, and hash stability."""

from pathlib import Path

import pytest

from posegate.config import load_config
from posegate.exceptions import ImmutableRecordError, InputValidationError
from posegate.policy import evaluate_policy
from posegate.records import build_run_identity, build_shadow_record, sha256_json
from posegate.trajectory import DomainProfile, PrefixMeasurement, PrefixSeries

from test_threshold_policy import threshold_mapping  # noqa: F401


ROOT = Path(__file__).resolve().parents[1]

# Frozen so a later schema change cannot silently alter how a bundled v1 policy
# hashes. These are the values the v1.0 code produced.
V1_SCIENTIFIC_HASHES = {
    "posegate_5ns_v1": (
        "33ff80cdc85fc9ad7d336376191c8f44588f3974dbe55a21fa1a5406c1640a5e"
    ),
    "posegate_20ns_v1": (
        "61bce3150a5b15fb3274f50f35dc1d013e353515f1ef23e6500d1dc82ad12386"
    ),
}


def synthetic_measurement(
    topology: Path, trajectory: Path, *, diameter_fraction: float = 0.2
) -> PrefixMeasurement:
    times = tuple(round(0.2 * (index + 1), 4) for index in range(25))
    rmsd = tuple(0.5 + 0.14 * index for index in range(25))
    window = tuple(3.0 < time <= 5.0001 for time in times)
    return PrefixMeasurement(
        topology_path=str(topology),
        trajectory_path=str(trajectory),
        alignment_selection="protein and name CA",
        ligand_selection="resname UNK and not name H*",
        alignment_atom_count=1343,
        ligand_atom_count=38,
        trajectory_frames_available=25,
        prefix_frames_used=25,
        frame_interval_ns=0.2,
        max_available_time_ns=5.0,
        checkpoint_ns=5.0,
        window_start_ns=3.0,
        window_end_ns=5.0,
        window_frame_count=sum(window),
        prefix_coordinates_sha256="0" * 64,
        corrected_pose_rmsd_mean_angstrom=4.168600883739357,
        corrected_centroid_displacement_mean_angstrom=3.6202337229997967,
        frame_relative_reorientation_mean_degrees=12.5,
        internal_deformation_mean_angstrom=0.4,
        domain_profile=DomainProfile(
            max_ligand_diameter_angstrom=13.199920499173938,
            min_box_length_angstrom=13.199920499173938 / diameter_fraction,
            ligand_diameter_box_fraction=diameter_fraction,
        ),
        series=PrefixSeries(
            time_ns=times,
            corrected_pose_rmsd_angstrom=rmsd,
            centroid_displacement_angstrom=tuple(value * 0.8 for value in rmsd),
            frame_relative_reorientation_degrees=tuple(10.0 for _ in rmsd),
            internal_deformation_angstrom=tuple(0.4 for _ in rmsd),
            in_feature_window=window,
        ),
    )


@pytest.fixture()
def sealed_inputs(tmp_path: Path) -> tuple[Path, Path]:
    topology = tmp_path / "system.pdb"
    trajectory = tmp_path / "trajectory.dcd"
    topology.write_bytes(b"REMARK synthetic\n")
    trajectory.write_bytes(b"synthetic")
    return topology, trajectory


def test_bundled_v1_policies_hash_exactly_as_before() -> None:
    for name, expected in V1_SCIENTIFIC_HASHES.items():
        config = load_config(ROOT / f"configs/{name}.yaml")
        assert sha256_json(config.scientific_dict()) == expected


def test_record_carries_identity_and_named_quality_checks(
    sealed_inputs: tuple[Path, Path],
) -> None:
    topology, trajectory = sealed_inputs
    config = load_config(ROOT / "configs/posegate_5ns_v1.yaml")
    measurement = synthetic_measurement(topology, trajectory)
    decision = evaluate_policy(config.policy, measurement.policy_measurements())
    record = build_shadow_record(
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

    assert record["schema_version"] == "1.1"
    assert record["run_identity"]["protein"] == "CIMG_06197"
    # Every identity key is present even when it was not supplied.
    assert record["run_identity"]["launched_utc"] is None
    assert record["policy_type"] == "standardized_logistic"

    quality = record["quality_control"]
    assert quality["overall"] == "PASS"
    assert quality["failed_checks"] == []
    assert quality["checks"]["future_frames_accessed"] is False
    # v1 declares no applicability block, so the domain check is not asserted.
    assert quality["checks"]["ligand_within_declared_domain"] is None
    assert quality["domain_profile"]["max_ligand_diameter_box_fraction"] is None

    assert len(record["measurement"]["series"]["time_ns"]) == 25


def test_identity_rejects_unknown_fields() -> None:
    with pytest.raises(ImmutableRecordError, match="unknown run identity"):
        build_run_identity(protein="A", organism="B")


def test_blank_identity_values_become_null() -> None:
    identity = build_run_identity(protein="  ", ligand="Milbemycin")
    assert identity["protein"] is None
    assert identity["ligand"] == "Milbemycin"


def test_domain_limit_is_recorded_when_configured(
    sealed_inputs: tuple[Path, Path],
) -> None:
    from posegate.config import config_from_mapping

    topology, trajectory = sealed_inputs
    config = config_from_mapping(threshold_mapping(), configuration_sha256="test")
    measurement = synthetic_measurement(topology, trajectory, diameter_fraction=0.2)
    decision = evaluate_policy(config.policy, measurement.policy_measurements())
    record = build_shadow_record(
        run_id="example", config=config, measurement=measurement, decision=decision
    )
    quality = record["quality_control"]
    assert quality["checks"]["ligand_within_declared_domain"] is True
    assert quality["domain_profile"]["max_ligand_diameter_box_fraction"] == 0.45
    assert record["policy_type"] == "threshold_rule"
    assert record["stop_threshold_angstrom"] == pytest.approx(2.7634173197224974)
    assert record["signed_stop_margin_angstrom"] > 0
    assert record["uncalibrated_retention_score"] is None
