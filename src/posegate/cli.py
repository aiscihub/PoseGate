"""PoseGate-MD command-line interface."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys
import time
from typing import Any

from . import __version__
from .config import load_config
from .exceptions import (
    CheckpointNotReady,
    ImmutableRecordError,
    InputValidationError,
    PoseGateError,
)
from .policy import evaluate_policy
from .records import (
    capture_path,
    create_shadow_record,
    load_record,
    sha256_file,
)
from .trajectory import inspect_inputs, measure_prefix
from .validation import validate_shadow_record


def _print_json(value: Any) -> None:
    print(json.dumps(value, indent=2, sort_keys=True))


def _write_json(path: str | Path, value: Any) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(output)


def _add_inputs(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--topology", type=Path, required=True)
    parser.add_argument("--trajectory", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="posegate",
        description=(
            "Coordinate-corrected, same-trajectory geometric forecasting for MD"
        ),
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    inspect_parser = subparsers.add_parser(
        "inspect", help="verify inputs, selections, PBC, and checkpoint coverage"
    )
    _add_inputs(inspect_parser)

    measure_parser = subparsers.add_parser(
        "measure", help="measure the configured coordinate prefix"
    )
    _add_inputs(measure_parser)
    measure_parser.add_argument("--output", type=Path)

    shadow_parser = subparsers.add_parser(
        "shadow", help="create one immutable forecast before late frames exist"
    )
    _add_inputs(shadow_parser)
    shadow_parser.add_argument("--run-id", required=True)
    shadow_parser.add_argument("--registry", type=Path, required=True)
    shadow_parser.add_argument("--checkpoint-file", type=Path)

    watch_parser = subparsers.add_parser(
        "watch", help="wait for a growing trajectory and capture once"
    )
    _add_inputs(watch_parser)
    watch_parser.add_argument("--run-id", required=True)
    watch_parser.add_argument("--registry", type=Path, required=True)
    watch_parser.add_argument("--checkpoint-file", type=Path)
    watch_parser.add_argument("--poll-seconds", type=float, default=30.0)

    validate_parser = subparsers.add_parser(
        "validate", help="compare a sealed forecast with a completed trajectory"
    )
    validate_parser.add_argument("--record", type=Path, required=True)
    validate_parser.add_argument("--topology", type=Path, required=True)
    validate_parser.add_argument("--trajectory", type=Path, required=True)
    validate_parser.add_argument("--output", type=Path)
    return parser


def _capture(args: argparse.Namespace) -> tuple[Path, dict[str, Any]]:
    config = load_config(args.config)
    existing = capture_path(args.registry, args.run_id, config.policy_id)
    if existing.exists():
        record = load_record(existing)
        expected = {
            "run_id": args.run_id,
            "policy_id": config.policy_id,
            "configuration_file_sha256": config.configuration_sha256,
            "topology_sha256": sha256_file(args.topology),
        }
        mismatches = {
            key: {"expected": value, "observed": record.get(key)}
            for key, value in expected.items()
            if record.get(key) != value
        }
        if mismatches:
            raise ImmutableRecordError(
                f"existing capture does not match this request: {mismatches}"
            )
        return existing, record
    measurement = measure_prefix(
        args.topology,
        args.trajectory,
        config,
        require_pre_outcome=True,
    )
    decision = evaluate_policy(config.policy, measurement.policy_measurements())
    return create_shadow_record(
        registry=args.registry,
        run_id=args.run_id,
        config=config,
        measurement=measurement,
        decision=decision,
        checkpoint_file=args.checkpoint_file,
    )


def run(args: argparse.Namespace) -> int:
    if args.command == "inspect":
        config = load_config(args.config)
        inspection = inspect_inputs(args.topology, args.trajectory, config)
        _print_json(asdict(inspection))
        return 0

    if args.command == "measure":
        config = load_config(args.config)
        measurement = measure_prefix(args.topology, args.trajectory, config)
        result = measurement.as_dict()
        if args.output:
            _write_json(args.output, result)
        _print_json(result)
        return 0

    if args.command == "shadow":
        path, record = _capture(args)
        _print_json(
            {
                "capture_path": str(path),
                "capture_sha256": sha256_file(path),
                "record": record,
            }
        )
        return 0

    if args.command == "watch":
        if args.poll_seconds <= 0:
            raise InputValidationError("--poll-seconds must be positive")
        while True:
            try:
                path, record = _capture(args)
                _print_json(
                    {
                        "capture_path": str(path),
                        "capture_sha256": sha256_file(path),
                        "record": record,
                    }
                )
                return 0
            except CheckpointNotReady as exc:
                print(f"waiting: {exc}", file=sys.stderr, flush=True)
                time.sleep(args.poll_seconds)
            except InputValidationError as exc:
                message = str(exc)
                if "could not read topology/trajectory" not in message:
                    raise
                print(f"waiting: {message}", file=sys.stderr, flush=True)
                time.sleep(args.poll_seconds)

    if args.command == "validate":
        record = load_record(args.record)
        result = validate_shadow_record(
            record, topology=args.topology, trajectory=args.trajectory
        )
        if args.output:
            _write_json(args.output, result)
        _print_json(result)
        return 0

    raise AssertionError(f"unhandled command: {args.command}")


def main() -> None:
    try:
        raise SystemExit(run(build_parser().parse_args()))
    except ImmutableRecordError as exc:
        print(f"immutable record error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    except PoseGateError as exc:
        print(f"posegate error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc


if __name__ == "__main__":
    main()
