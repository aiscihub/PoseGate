"""Immutable, self-contained shadow forecast records."""

from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time
from typing import Any

from . import __version__
from .config import PoseGateConfig
from .exceptions import ImmutableRecordError
from .policy import PolicyDecision
from .trajectory import PrefixMeasurement


SCOPE_WARNING = (
    "This forecast applies only to continuation of the trajectory already "
    "observed. It has not been validated for predicting a separately launched "
    "simulation."
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


def build_shadow_record(
    *,
    run_id: str,
    config: PoseGateConfig,
    measurement: PrefixMeasurement,
    decision: PolicyDecision,
    checkpoint_file: str | Path | None = None,
) -> dict[str, Any]:
    topology = Path(measurement.topology_path)
    trajectory = Path(measurement.trajectory_path)
    embedded_config = config.scientific_dict()
    checkpoint_path = Path(checkpoint_file).resolve() if checkpoint_file else None
    if checkpoint_path is not None and not checkpoint_path.is_file():
        raise ImmutableRecordError(f"checkpoint file does not exist: {checkpoint_path}")
    record = {
        "schema_version": "1.0",
        "record_id": f"posegate-{_slug(run_id)}-{_slug(config.policy_id)}",
        "created_utc": utc_now(),
        "run_id": run_id,
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
        "quality_control": "PASS",
        "forecast": decision.forecast,
        "recommendation": decision.recommendation,
        "action_applied": False,
        "uncalibrated_retention_score": (decision.uncalibrated_retention_score),
        "score_calibrated": decision.score_calibrated,
        "decision_threshold": decision.decision_threshold,
        "linear_predictor": decision.linear_predictor,
        "standardized_features": dict(decision.standardized_features),
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
) -> tuple[Path, dict[str, Any]]:
    path = capture_path(registry, run_id, config.policy_id)
    record = build_shadow_record(
        run_id=run_id,
        config=config,
        measurement=measurement,
        decision=decision,
        checkpoint_file=checkpoint_file,
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
