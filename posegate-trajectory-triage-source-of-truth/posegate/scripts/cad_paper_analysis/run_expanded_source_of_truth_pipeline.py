#!/usr/bin/env python3
"""Single entrypoint for the expanded trajectory source-of-truth workflow."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
DEFAULT_RAW_ROOT = Path("/media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md")
DEFAULT_RESULT_DIR = (
    REPO_ROOT
    / "pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth"
)


def run(command: list[str], env: dict[str, str] | None = None) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, env=env)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--reuse-measurements", action="store_true")
    parser.add_argument("--skip-ai", action="store_true")
    args = parser.parse_args()

    measurements = args.result_dir / "full_100ns_corrected_measurements.csv"
    master_dir = args.result_dir / "master"
    master_table = master_dir / "expanded_trajectory_source_of_truth.csv"
    primary_dir = args.result_dir / "primary_analysis"
    ai_dir = args.result_dir / "physics_guided_ai_audit"
    generator = HERE / "28_generate_full_100ns_measurements.py"

    args.result_dir.mkdir(parents=True, exist_ok=True)
    if not args.reuse_measurements:
        run(
            [
                sys.executable, str(generator), "--root", str(args.raw_root),
                "--output", str(measurements), "--workers", str(args.workers),
            ]
        )
    elif not measurements.exists():
        raise FileNotFoundError(f"--reuse-measurements requested but missing: {measurements}")

    run(
        [
            sys.executable, str(HERE / "31_build_expanded_source_of_truth.py"),
            "--measurements", str(measurements), "--output-dir", str(master_dir),
            "--generator", str(generator),
        ]
    )
    run(
        [
            sys.executable, str(HERE / "29_analyze_expanded_100ns_cohort.py"),
            "--input", str(master_table), "--output-dir", str(primary_dir),
        ],
        env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
    )
    if not args.skip_ai:
        run(
            [
                sys.executable, str(HERE / "30_physics_guided_ai_audit.py"),
                "--input", str(master_table), "--output-dir", str(ai_dir),
            ],
            env={**os.environ, "OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1"},
        )

    print(f"Expanded source-of-truth workflow complete: {args.result_dir}")


if __name__ == "__main__":
    main()
