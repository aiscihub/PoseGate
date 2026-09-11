"""Immutable, self-contained shadow forecast records."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
from typing import Any, Mapping

from . import __version__
from .config import PoseGateConfig
from .exceptions import ImmutableRecordError
from .policy import PolicyDecision
from .trajectory import PrefixMeasurement


SCOPE_WARNING = (
    "This shadow assessment applies only to the trajectory already observed. "
    "Advance warning before observable deviation and transfer to an independently "
    "launched simulation have not been established."
)


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_json(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def _slug(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "-", value).strip("-")
    if not slug:
        raise ImmutableRecordError("run_id must contain a filename-safe character")
    return slug


def capture_path(registry: str | Path, run_id: str, policy_id: str) -> Path:
    return Path(registry) / f"{_slug(run_id)}__{_slug(policy_id)}.json"


def write_immutable_json(path: str | Path, value: dict[str, Any]) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL
    try:
        fd = os.open(output, flags, 0o644)
    except FileExistsError as exc:
        raise ImmutableRecordError(
            f"immutable record already exists: {output}"
        ) from exc
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
    except Exception:
        output.unlink(missing_ok=True)
        raise
    output.chmod(0o444)


def _git_snapshot(start: Path) -> dict[str, Any]:
    current = start.resolve()
    repository: Path | None = None
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists():
            repository = candidate
            break
    if repository is None:
        return {"git_commit": None, "repository_dirty": None}
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repository,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout.strip()
        status = subprocess.run(
            ["git", "status", "--porcelain=v1", "--untracked-files=all"],
            cwd=repository,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return {"git_commit": None, "repository_dirty": None}
    return {
        "git_commit": commit,
        "repository_dirty": bool(status),
        "git_status_sha256": hashlib.sha256(status.encode()).hexdigest(),
    }


IDENTITY_FIELDS = ("protein", "ligand", "pocket", "replica", "launched_utc")


def build_run_identity(**values: str | None) -> dict[str, str | None]:
    """Return the structured launch identity, with every key always present.

    ``run_id`` stays the addressable primary key. This block exists so a
    registry can group and label records without parsing that free-text string.
    """
    unknown = sorted(set(values).difference(IDENTITY_FIELDS))
    if unknown:
        raise ImmutableRecordError(f"unknown run identity fields: {unknown}")
    identity: dict[str, str | None] = {}
    for field in IDENTITY_FIELDS:
        value = values.get(field)
        if value is None:
            identity[field] = None
            continue
        text = str(value).strip()
        identity[field] = text or None
    return identity


def _quality_control(
    config: PoseGateConfig, measurement: PrefixMeasurement
) -> dict[str, Any]:
    """Name the checks the measurement path already enforced.

    Every entry is a precondition of reaching this point: the record fails to
    exist at all if one of them is violated. Naming them individually replaces
    an opaque ``"PASS"`` string and lets a reader see which checks were in
    force under this configuration.
    """
    domain = measurement.domain_profile
    limit = config.applicability
    checks: dict[str, Any] = {
        "periodic_box_vectors_valid": True,
        "protein_ligand_selections_disjoint": True,
        "checkpoint_coverage_complete": (
            measurement.prefix_frames_used == len(measurement.series.time_ns)
        ),
        "feature_window_populated": measurement.window_frame_count > 0,
        "future_frames_accessed": False,
    }
    if limit is None:
        checks["ligand_within_declared_domain"] = None
    else:
        checks["ligand_within_declared_domain"] = (
            domain.ligand_diameter_box_fraction
            <= limit.max_ligand_diameter_box_fraction
        )
    failed = sorted(
        name
        for name, value in checks.items()
        if value is False and name != "future_frames_accessed"
    )
    return {
        "overall": "FAIL" if failed else "PASS",
        "failed_checks": failed,
        "checks": checks,
        "domain_profile": {
            "max_ligand_diameter_angstrom": domain.max_ligand_diameter_angstrom,
            "min_box_length_angstrom": domain.min_box_length_angstrom,
            "ligand_diameter_box_fraction": domain.ligand_diameter_box_fraction,
            "max_ligand_diameter_box_fraction": (
                None if limit is None else limit.max_ligand_diameter_box_fraction
            ),
        },
    }


def build_shadow_record(
    *,
    run_id: str,
    config: PoseGateConfig,
    measurement: PrefixMeasurement,
    decision: PolicyDecision,
    checkpoint_file: str | Path | None = None,
    run_identity: Mapping[str, str | None] | None = None,
) -> dict[str, Any]:
    topology = Path(measurement.topology_path)
    trajectory = Path(measurement.trajectory_path)
    embedded_config = config.scientific_dict()
    checkpoint_path = Path(checkpoint_file).resolve() if checkpoint_file else None
    if checkpoint_path is not None and not checkpoint_path.is_file():
        raise ImmutableRecordError(f"checkpoint file does not exist: {checkpoint_path}")
    record = {
        "schema_version": "1.1",
        "record_id": f"posegate-{_slug(run_id)}-{_slug(config.policy_id)}",
        "created_utc": utc_now(),
        "run_id": run_id,
        "run_identity": (
            build_run_identity() if run_identity is None else dict(run_identity)
        ),
        "software": {"name": "posegate-md", "version": __version__},
        "configuration": embedded_config,
        "configuration_file_sha256": config.configuration_sha256,
        "configuration_scientific_sha256": sha256_json(embedded_config),
        "policy_id": config.policy_id,
        "policy_mode": config.policy.mode,
        "topology_path": str(topology),
        "topology_sha256": sha256_file(topology),
        "trajectory_path": str(trajectory),
        "trajectory_size_bytes_at_capture": trajectory.stat().st_size,
        "trajectory_prefix_coordinates_sha256": (measurement.prefix_coordinates_sha256),
        "checkpoint_file_path": str(checkpoint_path) if checkpoint_path else None,
        "checkpoint_file_sha256": (
            sha256_file(checkpoint_path) if checkpoint_path else None
        ),
        "measurement": measurement.as_dict(),
        "quality_control": _quality_control(config, measurement),
        "forecast": decision.forecast,
        "recommendation": decision.recommendation,
        "action_applied": False,
        "policy_type": decision.policy_type,
        "feature_values": dict(decision.feature_values),
        "uncalibrated_retention_score": (decision.uncalibrated_retention_score),
        "score_calibrated": decision.score_calibrated,
        "decision_threshold": decision.decision_threshold,
        "linear_predictor": decision.linear_predictor,
        "standardized_features": (
            None
            if decision.standardized_features is None
            else dict(decision.standardized_features)
        ),
        "stop_threshold_angstrom": decision.stop_threshold_angstrom,
        "signed_stop_margin_angstrom": decision.signed_stop_margin_angstrom,
        "late_outcome_accessed": False,
        "outcome_embargo_until_ns": config.outcome.window_start_ns,
        "scope": "FORECASTS_THIS_TRAJECTORY_ONLY",
        "scope_warning": SCOPE_WARNING,
        "provenance": dict(config.provenance),
        **_git_snapshot(Path.cwd()),
    }
    return record


def create_shadow_record(
    *,
    registry: str | Path,
    run_id: str,
    config: PoseGateConfig,
    measurement: PrefixMeasurement,
    decision: PolicyDecision,
    checkpoint_file: str | Path | None = None,
    run_identity: Mapping[str, str | None] | None = None,
) -> tuple[Path, dict[str, Any]]:
    path = capture_path(registry, run_id, config.policy_id)
    record = build_shadow_record(
        run_id=run_id,
        config=config,
        measurement=measurement,
        decision=decision,
        checkpoint_file=checkpoint_file,
        run_identity=run_identity,
    )
    write_immutable_json(path, record)
    return path, record


def load_record(path: str | Path) -> dict[str, Any]:
    record_path = Path(path)
    try:
        value = json.loads(record_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ImmutableRecordError(
            f"could not read shadow record {record_path}: {exc}"
        ) from exc
    if not isinstance(value, dict) or "record_id" not in value:
        raise ImmutableRecordError(f"invalid shadow record: {record_path}")
    return value
