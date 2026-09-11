#!/usr/bin/env python3
"""Export the current manuscript and its numerical inputs without changing archives."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def export(source: Path, destination: Path) -> None:
    manifest = {}

    def copy(relative: str, target: str) -> None:
        origin, output = source / relative, destination / target
        payload = origin.read_bytes()
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_bytes(payload)
        manifest[target] = {
            "source": relative,
            "sha256": hashlib.sha256(payload).hexdigest(),
        }

    def generated(target, value, provenance):
        output = destination / target
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(value, indent=2) + "\n")
        manifest[target] = {
            "source": provenance,
            "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
        }

    doc = "doc/YAU_competition_2026"
    tables = {}
    figures = {}
    source_documents = {}
    for name in ("yau_aj", "yau_supplementary_materials"):
        for suffix in ("tex", "pdf"):
            relative = f"{doc}/{name}.{suffix}"
            source_documents[relative] = {
                "sha256": hashlib.sha256((source / relative).read_bytes()).hexdigest(),
                "distributed": False,
            }
        text = (source / doc / f"{name}.tex").read_text()
        text = re.sub(r"(?<!\\)%[^\n]*", "", text)
        for asset in re.findall(
            r"\\includegraphics(?:\[[^\]]*\])?\s*\{([^}]+)\}", text
        ):
            copy(f"{doc}/{asset}", asset)
        for match in re.finditer(r"\\begin\{(figure\*?)\}.*?\\end\{\1\}", text, re.S):
            labels = re.findall(r"\\label\{([^}]+)\}", match.group())
            assets = re.findall(
                r"\\includegraphics(?:\[[^\]]*\])?\s*\{([^}]+)\}", match.group()
            )
            if labels:
                figures[labels[0]] = {"source": f"{name}.tex", "assets": assets}
        for match in re.finditer(
            r"\\begin\{(table\*?|longtable)\}.*?\\end\{\1\}", text, re.S
        ):
            labels = re.findall(r"\\label\{([^}]+)\}", match.group())
            if labels:
                tables[labels[0]] = {
                    "source": f"{name}.tex",
                    "aliases": labels[1:],
                    "latex": match.group(),
                }
    if len(tables) != 13:
        raise ValueError(
            f"Expected 6 main and 7 supplementary tables, found {len(tables)}"
        )

    generated(
        "figure_index.json",
        figures,
        "Figure labels and assets extracted from the source manuscript",
    )
    directories = {
        "current_statistics_v2": [
            "analysis_rows.csv",
            "statistics.json",
            "crossligand_rows.csv",
            "crossligand_statistics.json",
            "local_pocket_rows.csv",
            "local_pocket_statistics.json",
            "crossing_times.csv",
            "crossing_summary.json",
        ],
        "final_protein_frame_control": [
            "local_checkpoint_measurements.csv",
            "primary_temporal_statistics.json",
            "protein_side_stability_statistics.json",
            "protein_side_stability_summary.csv",
        ],
        "current_separate_statistics_v2": [
            "current_separate_paired_rows.csv",
            "separate_launch_statistics.json",
        ],
        "current_boundary_sensitivity": [
            "boundary_sensitivity.csv",
            "boundary_sensitivity.json",
        ],
        "current_primary_v2": ["primary_v2_statistics.json"],
    }
    for directory, files in directories.items():
        for filename in files:
            copy(
                f"yau_competition/reanalysis/{directory}/{filename}",
                f"data/{directory}/{filename}",
            )
    for cohort, directory, prefix in (
        (
            "prospective_validation",
            "fresh100ns_preliminary_analysis_n24",
            "confirmation_v1_fresh100ns_preliminary_n24",
        ),
        ("prospective_pilot", "pilot20ns_analysis", "confirmation_v1_pilot20ns"),
    ):
        for kind in (
            "measurements.csv",
            "metrics.csv",
            "correlations.csv",
            "summary.json",
        ):
            copy(
                f"yau_competition/prospective_shadow/confirmation_v1/{directory}/{prefix}_{kind}",
                f"data/{cohort}/{kind}",
            )
    archived = "pipeline/training/outputs/yau_recapture_screen/cad_states/expanded_source_of_truth_frozen_n81"
    copy(
        f"{archived}/contact_retention_outcome_n81.csv",
        "data/secondary/contact_retention.csv",
    )
    for name in (
        "consolidated_100ns_extension_n81.csv",
        "local_pocket_control_measurements_n81_fixed.csv",
        "local_pocket_control_results_n81_fixed.csv",
        "local_pocket_control_theta_agreement.csv",
        "two_fate_example_prpd_windows.csv",
        "master/expanded_trajectory_source_of_truth.csv",
        "protein_side_stability/protein_side_stability_frames.csv",
    ):
        copy(f"{archived}/{name}", f"data/secondary/{name}")
    for name in (
        "current_100ns_manifest.json",
        "current_crossligand_manifest.json",
        "current_separate_manifest.json",
    ):
        copy(f"yau_competition/reanalysis/{name}", f"cohorts/{name}")
    cases = json.loads(
        (source / "yau_competition/reanalysis/current_100ns_manifest.json").read_text()
    )["cases"]
    if len(cases) != 81 or sum(c["primary"] for c in cases) != 52:
        raise ValueError("Frozen cohort membership changed")
    final_cases = []
    for case in cases:
        key = case["cohort"] + "__" + case["complex_id"].replace("|", "__")
        prefix = f"yau_competition/reanalysis/current_frames_v2/{key}"
        for name in ("frames.csv", "measurement.json"):
            copy(f"{prefix}/{name}", f"data/current_frames_v2/{key}/{name}")
        measurement = json.loads((source / prefix / "measurement.json").read_text())
        if measurement["status"] != "ok":
            raise ValueError(f"Incomplete measurement: {key}")
        final_cases.append(
            {
                **{
                    k: case[k]
                    for k in (
                        "complex_id",
                        "protein",
                        "ligand",
                        "pocket",
                        "replica",
                        "primary",
                        "topology_path",
                        "trajectory_path",
                    )
                },
                "measurements": f"data/current_frames_v2/{key}/measurement.json",
                "frames": f"data/current_frames_v2/{key}/frames.csv",
                "topology_sha256": measurement["topology_sha256"],
                "selected_coordinates_sha256": measurement[
                    "selected_coordinates_sha256"
                ],
                "all_measured_selected_atoms_sha256": measurement[
                    "all_measured_selected_atoms_sha256"
                ],
                "frames_csv_sha256": measurement["frames_csv_sha256"],
                "X5_A": measurement["early"]["corrected_A"],
                "L100_A": measurement["late"]["corrected_A"],
                "retained": measurement["retained"],
            }
        )
    generated(
        "cohorts/final_expanded_81.json",
        {
            "scope": "Final whole-ligand, nominal-window measurements; 81 trajectories, 52 primary, 16 proteins",
            "raw_availability": "Complete trajectories are not distributed because of file size.",
            "hash_scope": "Coordinate digests cover selected atoms, boxes and indices, not full raw DCD bytes.",
            "cases": final_cases,
        },
        "current_100ns_manifest.json joined to current_frames_v2 measurement.json by complex_id",
    )
    scripts = (
        "summarize_current_frames.py",
        "reanalyze_current_primary.py",
        "shadow_v2_core.py",
        "remeasure_current_frames.py",
        "remeasure_current_crossligand_frames.py",
        "summarize_current_crossligand.py",
        "remeasure_current_separate_launch.py",
        "summarize_current_separate_launch.py",
        "remeasure_current_local_pocket.py",
        "finalize_protein_frame_control.py",
        "finalize_supplement_temporal.py",
        "75_protein_side_stability_audit.py",
        "current_boundary_sensitivity.py",
        "77_analyze_confirmation_v1_fresh100ns_preliminary.py",
        "make_current_version_figures.py",
        "make_figure3_early_to_late_and_checkpoints.py",
        "make_figure3_physical_attribution_combined.py",
        "make_figure4_early_detection_redesign.py",
        "common.py",
        "35_local_pocket_alignment_control.py",
        "36_local_pocket_control_evaluation.py",
        "48b_priority2_horizon_and_contact_analysis_n81.py",
        "48_priority2_horizon_and_contact_analysis.py",
        "43_priority1_reanalysis.py",
        "26_plot_locked_claim_results.py",
    )
    for name in (
        "76_analyze_confirmation_v1_pilot20ns.py",
        "55_select_shadow_threshold_v2.py",
    ):
        copy(
            f"yau_competition/scripts/cad_paper_analysis/{name}",
            f"source_analysis/{name}",
        )
    for name in (
        "frozen_policy/same_trajectory_5ns_rmsd_threshold_v2.json",
        "audits/threshold_selection_v2.json",
        "audits/strict_remeasurement_v2/strict_development_rows_v2.csv",
        "audits/strict_remeasurement_v2/strict_remeasurement_summary_v2.json",
    ):
        copy(
            f"yau_competition/prospective_shadow/v2/{name}",
            "policy/" + name.split("/", 1)[1],
        )
    copy(
        "yau_competition/prospective_shadow/confirmation_v1/PROSPECTIVE_CONFIRMATION_V1_MANIFEST.json",
        "cohorts/PROSPECTIVE_CONFIRMATION_V1_MANIFEST.json",
    )
    for name in scripts:
        copy(
            f"yau_competition/scripts/cad_paper_analysis/{name}",
            f"source_analysis/{name}",
        )
    for name in (
        "checkpoint_plot_data.csv",
        "descriptor_plot_data.csv",
        "descriptor_paired_plot_data.csv",
    ):
        copy(f"{doc}/figures/current_v2/{name}", f"figures/current_v2/{name}")
    for name in (
        "select_and_extract_frames.py",
        "contact_frequency.py",
        "METHODS.md",
        "v3/compose_two_fate_figure_v3.py",
        "v3/render_two_fate_pocket_poses_v3.py",
        "v3/renders/render_meta_retained.json",
        "v3/renders/render_meta_nonretained.json",
    ):
        copy(
            f"{doc}/figures/two_fate_pocket_poses/{name}",
            f"source_visualization/{name}",
        )
    for tag in ("retained", "nonretained"):
        for name in (
            "selection_meta.json",
            "contact_table.json",
            "frame_init.pdb",
            "frame_5ns.pdb",
            "frame_late.pdb",
        ):
            copy(
                f"{doc}/figures/data/{tag}/{name}",
                f"source_visualization/data/{tag}/{name}",
            )
    (destination / "table_index.json").write_text(json.dumps(tables, indent=2) + "\n")
    manifest["table_index.json"] = {
        "source": "Extracted from the two included TeX sources; no numerical rewriting",
        "sha256": hashlib.sha256(
            (destination / "table_index.json").read_bytes()
        ).hexdigest(),
    }
    (destination / "manifest.json").write_text(
        json.dumps(
            {
                "scope": "Current manuscript table/figure extracts and source results; not new prospective observations",
                "files": manifest,
                "source_documents": source_documents,
                "notes": [
                    "Full manuscript PDF and TeX files are not distributed; source hashes bind the extracted tables and figure index.",
                    "Historical filenames and source metadata are preserved, not renamed scientific cohorts.",
                    "The current validation table counts nominal remaining time from 5 ns (855 ns), not the archived 6 ns accounting (846 ns).",
                    "The contact-neighborhood check remains an archived secondary analysis.",
                    "Source-analysis scripts preserve original implementations; the portable verifier calls their statistical helpers only.",
                ],
            },
            indent=2,
        )
        + "\n"
    )
    print(f"Exported {len(manifest)} files and {len(tables)} tables to {destination}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    args = parser.parse_args()
    export(args.source_root.resolve(), ROOT / "paper/current")
