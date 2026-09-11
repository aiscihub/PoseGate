#!/usr/bin/env python3
"""Export numerical table CSVs and manuscript mapping from frozen source results."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper/current"
DATA = PAPER / "data"
OUT = PAPER / "tables"


def read(relative):
    return json.loads((DATA / relative).read_text())


def write_csv(path, rows):
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def export_tables():
    OUT.mkdir(parents=True, exist_ok=True)
    tables = json.loads((PAPER / "table_index.json").read_text())
    main_n = supp_n = 0
    for entry in tables.values():
        if entry["source"] == "yau_aj.tex":
            main_n += 1
            name = f"main_table_{main_n}"
        else:
            supp_n += 1
            name = f"supp_table_S{supp_n}"
        (OUT / (name + ".tex")).write_text(entry["latex"] + "\n")
    stats = read("current_statistics_v2/statistics.json")
    rows = pd.read_csv(DATA / "current_statistics_v2/analysis_rows.csv")
    cross = pd.read_csv(DATA / "current_statistics_v2/crossligand_rows.csv")
    pilot = pd.read_csv(DATA / "prospective_pilot/measurements.csv")
    validation = pd.read_csv(DATA / "prospective_validation/measurements.csv")
    paired = pd.read_csv(
        DATA / "current_separate_statistics_v2/current_separate_paired_rows.csv"
    )
    cohort_rows = []
    for name, data, length, late, unit in (
        (
            "Primary milbemycin",
            rows.loc[rows.primary],
            "100",
            "(70,100]",
            "trajectories",
        ),
        ("Expanded milbemycin", rows, "100", "(70,100]", "trajectories"),
        ("Prospective validation", validation, "100", "(70,100]", "trajectories"),
        ("Prospective pilot", pilot, "20", "(14,20]", "trajectories"),
        (
            "Separate-launch control",
            paired,
            "20 source; 100 target",
            "(70,100]",
            "matched pairs",
        ),
        (
            "Beauvericin",
            cross.loc[cross.cohort == "beauvericin"],
            "20",
            "(14,20]",
            "trajectories",
        ),
        (
            "Verapamil",
            cross.loc[cross.cohort == "verapamil"],
            "20",
            "(14,20]",
            "trajectories",
        ),
    ):
        cohort_rows.append(
            dict(
                cohort=name,
                n=len(data),
                proteins=data.protein.nunique(),
                unit=unit,
                trajectory_length_ns=length,
                early_window_ns="(3,5]",
                late_window_ns=late,
            )
        )
    write_csv(OUT / "main_table_1.csv", cohort_rows)
    result = []
    for score, name in (
        ("corrected_A", "Corrected protein-relative RMSD"),
        ("pbc_naive_A", "PBC-naive RMSD"),
        ("ligand_self_A", "Ligand-self-aligned RMSD"),
    ):
        late = rows.loc[rows.primary, "late_" + score]
        row = dict(
            representation=name,
            late_min_A=late.min(),
            late_max_A=late.max(),
            pooled_n=52,
            within_n=50,
            within_proteins=13,
        )
        for kind in ("pooled", "within"):
            row.update(
                {
                    kind + "_" + k: v
                    for k, v in stats["correlations"][score][kind].items()
                }
            )
        result.append(row)
    write_csv(OUT / "main_table_4.csv", result)
    vm = pd.read_csv(DATA / "prospective_validation/metrics.csv").set_index("metric_id")
    vc = pd.read_csv(DATA / "prospective_validation/correlations.csv")
    assert len(vc) == 1
    association = vc.iloc[0]
    correlation_ci = json.loads(association.ci)
    result = [
        dict(
            panel="A",
            analysis="Full-cohort Spearman X5 vs L100",
            estimate=association.rho,
            ci_low=correlation_ci[0],
            ci_high=correlation_ci[1],
            n=24,
            proteins=16,
            retained=9,
            nonretained=15,
        )
    ]
    for key, label in (
        ("corrected_X5", "Full-cohort AUROC"),
        ("corrected_X5_strict_rolling_prefix_clean", "Rolling-prefix-clean AUROC"),
    ):
        metric = vm.loc[key]
        result.append(
            dict(
                panel="A",
                analysis=label,
                estimate=metric.auroc,
                ci_low=metric.ci_lo,
                ci_high=metric.ci_hi,
                n=metric.n,
                proteins=metric.n_proteins,
                retained=metric.n_retained,
                nonretained=metric.n_nonretained,
            )
        )
    flagged = validation.corrected_rmsd_X5_mean_3_5_A > 2.7634173197224974
    y = validation.late_nonretained_70_100
    for flag, label in ((flagged, "Flagged"), (~flagged, "Unflagged")):
        result.append(
            dict(
                panel="B",
                analysis=label,
                retained=int((flag & ~y).sum()),
                nonretained=int((flag & y).sum()),
                n=int(flag.sum()),
            )
        )
    result.append(
        dict(
            panel="B",
            analysis="Total",
            retained=int((~y).sum()),
            nonretained=int(y.sum()),
            n=len(y),
        )
    )
    for label, value, unit in (
        ("Potential remaining time", flagged.sum() * 95, "ns"),
        (
            "Potential remaining time from successful stops",
            (flagged & y).sum() * 95,
            "ns",
        ),
        ("Potential remaining time from false stops", (flagged & ~y).sum() * 95, "ns"),
        (
            "Potential fraction of total trajectory time",
            flagged.sum() * 95 / 2400 * 100,
            "percent",
        ),
        ("Actual saved time", 0, "ns"),
    ):
        result.append(dict(panel="B", analysis=label, estimate=float(value), unit=unit))
    write_csv(OUT / "main_table_5.csv", result)
    result = []
    for name, data, late in (
        ("Primary milbemycin", stats["primary"]["metrics"], "(70,100]"),
        ("Prospective pilot", None, "(14,20]"),
        ("Prospective validation", None, "(70,100]"),
        ("Beauvericin", None, "(14,20]"),
        ("Verapamil", None, "(14,20]"),
    ):
        if name.startswith("Prospective"):
            folder = (
                "prospective_pilot"
                if name.endswith("pilot")
                else "prospective_validation"
            )
            metrics = pd.read_csv(DATA / folder / "metrics.csv").set_index("metric_id")
            values = [
                metrics.loc[k, "auroc"]
                for k in ("corrected_X5", "pbc_naive_X5", "ligand_self_X5")
            ]
            n = 24
        elif data is not None:
            values = [
                data[k]["estimate"]
                for k in ("corrected_A", "pbc_naive_A", "ligand_self_A")
            ]
            n = 52
        else:
            from verify_paper_results import auc

            subset = cross.loc[cross.cohort == name.lower()]
            values = [
                auc(subset.nonretained, subset[k])
                for k in ("corrected_A", "pbc_naive_A", "ligand_self_A")
            ]
            n = len(subset)
        result.append(
            dict(
                cohort=name,
                n=n,
                late_window_ns=late,
                corrected_auroc=values[0],
                pbc_naive_auroc=values[1],
                ligand_self_auroc=values[2],
            )
        )
    write_csv(OUT / "main_table_6.csv", result)
    temporal = read("final_protein_frame_control/primary_temporal_statistics.json")
    result = []
    for name, data, subsets in (
        ("Expanded milbemycin", rows, stats["subsets"]),
        ("Primary milbemycin", rows.loc[rows.primary], temporal["subsets"]),
    ):
        result.append(
            dict(
                cohort=name,
                risk_set="Full cohort",
                n=len(data),
                proteins=data.protein.nunique(),
                retained=int((~data.nonretained.astype(bool)).sum()),
                nonretained=int(data.nonretained.sum()),
            )
        )
        for risk, statistics in subsets.items():
            subset = data.loc[data[risk]]
            row = dict(
                cohort=name,
                risk_set=risk,
                n=len(subset),
                proteins=subset.protein.nunique(),
                retained=int((~subset.nonretained.astype(bool)).sum()),
                nonretained=int(subset.nonretained.sum()),
            )
            for metric in ("corrected_A", "reorientation_deg"):
                row.update(
                    {
                        metric + "_" + k: v
                        for k, v in statistics["metrics"][metric].items()
                        if k in ("estimate", "ci_low", "ci_high")
                    }
                )
            result.append(row)
    write_csv(OUT / "supp_table_S2.csv", result)
    write_csv(
        OUT / "supp_table_S3.csv",
        [dict(metric=k, value=v) for k, v in temporal["primary_timing"].items()],
    )
    result = []
    for checkpoint, value in read("current_statistics_v2/local_pocket_statistics.json")[
        "checkpoints"
    ].items():
        row = dict(checkpoint_ns=checkpoint)
        for frame in ("whole", "local"):
            row.update({frame + "_" + k: v for k, v in value[frame].items()})
        row.update(
            {k: value[k] for k in ("delta_local_minus_whole", "pearson", "spearman")}
        )
        result.append(row)
    write_csv(OUT / "supp_table_S4.csv", result)
    result = []
    for key, value in read(
        "final_protein_frame_control/protein_side_stability_statistics.json"
    )["metrics"].items():
        if key == "protein_ca_rmsd_mean_8_10_A":
            continue
        row = dict(quantity=key)
        for k, v in value.items():
            if isinstance(v, dict):
                row.update({k + "_" + field: n for field, n in v.items()})
            else:
                row[k] = v
        result.append(row)
    write_csv(OUT / "supp_table_S5.csv", result)
    for source, target in (
        ("secondary/contact_retention.csv", "supp_table_S6.csv"),
        ("current_boundary_sensitivity/boundary_sensitivity.csv", "supp_table_S7.csv"),
    ):
        (OUT / target).write_bytes((DATA / source).read_bytes())


def export_mapping():
    p = "paper/current/"
    s = p + "source_analysis/"
    d = p + "data/"
    t = p + "tables/"
    entries = []

    def add(item, label, scripts, inputs, outputs, scope):
        entries.append(
            dict(
                item=item,
                label=label,
                scripts=scripts,
                inputs=inputs,
                outputs=outputs,
                audit_scope=scope,
            )
        )

    add(
        "Table 1",
        "tab:early-late-link-summary",
        "tools/export_manuscript_tables.py",
        p
        + "cohorts/final_expanded_81.json; "
        + d
        + "current_statistics_v2/analysis_rows.csv; "
        + d
        + "prospective_validation/measurements.csv; "
        + d
        + "prospective_pilot/measurements.csv; "
        + d
        + "current_separate_statistics_v2/current_separate_paired_rows.csv; "
        + d
        + "current_statistics_v2/crossligand_rows.csv",
        t + "main_table_1.csv; " + t + "main_table_1.tex",
        "Composition is counted from cohort rows; evidential roles and screened-candidate note are authored in the TeX.",
    )
    add(
        "Table 2",
        "tab:early-late-measurements",
        "src/posegate/trajectory.py",
        p
        + "tables/main_table_2.tex; configs/posegate_5ns_v2.yaml; configs/posegate_pilot20ns_v2.yaml",
        t + "main_table_2.tex",
        "Definition table, not an estimated result; no numerical CSV.",
    )
    add(
        "Table 3",
        "tab:early-clean-definitions",
        s + "remeasure_current_frames.py; " + s + "summarize_current_frames.py",
        p + "tables/main_table_3.tex; " + d + "current_frames_v2/",
        t + "main_table_3.tex",
        "Authored eligibility definitions; derived flags in analysis_rows.csv and estimates in S2.",
    )
    add(
        "Table 4",
        "tab:representation-matched",
        s + "summarize_current_frames.py",
        d + "current_statistics_v2/analysis_rows.csv",
        d
        + "current_statistics_v2/statistics.json; "
        + t
        + "main_table_4.csv; "
        + t
        + "main_table_4.tex",
        "Recompute pooled and within-protein correlations from early/late rows; retain clustered CI provenance.",
    )
    add(
        "Table 5",
        "tab:prospective-validation",
        s
        + "77_analyze_confirmation_v1_fresh100ns_preliminary.py; tools/export_manuscript_tables.py; src/posegate/policy.py",
        d + "prospective_validation/measurements.csv; configs/posegate_5ns_v2.yaml",
        d
        + "prospective_validation/metrics.csv; "
        + d
        + "prospective_validation/correlations.csv; "
        + t
        + "main_table_5.csv; "
        + t
        + "main_table_5.tex",
        "Current table uses 5 ns resource accounting (855 ns); preserved source summary uses 6 ns (846 ns). Actual savings zero; unflagged is unresolved.",
    )
    add(
        "Table 6",
        "tab:cross_ligand_results",
        s
        + "summarize_current_frames.py; "
        + s
        + "summarize_current_crossligand.py; "
        + s
        + "77_analyze_confirmation_v1_fresh100ns_preliminary.py",
        d
        + "current_statistics_v2/analysis_rows.csv; "
        + d
        + "current_statistics_v2/crossligand_rows.csv; "
        + d
        + "prospective_pilot/measurements.csv; "
        + d
        + "prospective_validation/measurements.csv",
        t + "main_table_6.csv; " + t + "main_table_6.tex",
        "AUROC uses the same corrected late class while the early coordinate representation changes.",
    )
    add(
        "Table S1",
        "tab:simulation-protocol",
        "Authored preparation/protocol record; no statistical generator",
        p + "tables/supp_table_S1.tex",
        t + "supp_table_S1.tex",
        "Protocol settings and references are documented, not inferred from result tables; no numerical CSV.",
    )
    add(
        "Table S2",
        "tab:prefix-clean-risksets",
        s
        + "remeasure_current_frames.py; "
        + s
        + "summarize_current_frames.py; "
        + s
        + "finalize_supplement_temporal.py",
        d + "current_frames_v2/; " + d + "current_statistics_v2/analysis_rows.csv",
        d
        + "current_statistics_v2/statistics.json; "
        + d
        + "final_protein_frame_control/primary_temporal_statistics.json; "
        + t
        + "supp_table_S2.csv; "
        + t
        + "supp_table_S2.tex",
        "Expanded 81 to 46/45/31; primary 52 to 30/29/20. Alias tab:prefix-clean-auroc.",
    )
    add(
        "Table S3",
        "tab:time-to-event",
        s + "summarize_current_frames.py; " + s + "finalize_supplement_temporal.py",
        d + "current_frames_v2/; " + d + "current_statistics_v2/crossing_times.csv",
        d
        + "final_protein_frame_control/primary_temporal_statistics.json; "
        + t
        + "supp_table_S3.csv; "
        + t
        + "supp_table_S3.tex",
        "Retrospective sustained crossing, not demonstrated predictive lead time; denominator 39 primary non-retained trajectories.",
    )
    add(
        "Table S4",
        "tab:pocket-frame-control",
        s + "remeasure_current_local_pocket.py; " + s + "summarize_current_frames.py",
        d
        + "current_statistics_v2/local_pocket_rows.csv; "
        + d
        + "current_statistics_v2/analysis_rows.csv",
        d
        + "current_statistics_v2/local_pocket_statistics.json; "
        + t
        + "supp_table_S4.csv; "
        + t
        + "supp_table_S4.tex",
        "Current frozen 8 angstrom frame. Archived 6/8/10 angstrom radius sensitivity is listed separately below.",
    )
    add(
        "Table S5",
        "tab:protein-frame-stability",
        s
        + "finalize_protein_frame_control.py; "
        + s
        + "remeasure_current_local_pocket.py; "
        + s
        + "75_protein_side_stability_audit.py",
        d
        + "final_protein_frame_control/local_checkpoint_measurements.csv; "
        + d
        + "current_statistics_v2/analysis_rows.csv; "
        + d
        + "secondary/protein_side_stability/protein_side_stability_frames.csv",
        d
        + "final_protein_frame_control/protein_side_stability_summary.csv; "
        + d
        + "final_protein_frame_control/protein_side_stability_statistics.json; "
        + t
        + "supp_table_S5.csv; "
        + t
        + "supp_table_S5.tex",
        "C-alpha fits reuse archived per-frame traces with current exact-index aggregation; local ligand scores are recomputed. Half-up display rounds 0.8125 to 0.813.",
    )
    add(
        "Table S6",
        "tab:contact-retention-control",
        s
        + "48b_priority2_horizon_and_contact_analysis_n81.py; "
        + s
        + "43_priority1_reanalysis.py",
        d
        + "secondary/master/expanded_trajectory_source_of_truth.csv; "
        + d
        + "secondary/consolidated_100ns_extension_n81.csv",
        d
        + "secondary/contact_retention.csv; "
        + t
        + "supp_table_S6.csv; "
        + t
        + "supp_table_S6.tex",
        "Archived secondary contact measurement, not silently relabeled as current remeasurement. Available-pair denominators 80/81/80/80.",
    )
    add(
        "Table S7",
        "tab:boundary-sensitivity",
        s
        + "current_boundary_sensitivity.py; "
        + s
        + "reanalyze_current_primary.py; "
        + s
        + "77_analyze_confirmation_v1_fresh100ns_preliminary.py",
        d
        + "current_statistics_v2/analysis_rows.csv; "
        + d
        + "prospective_validation/measurements.csv",
        d
        + "current_boundary_sensitivity/boundary_sensitivity.json; "
        + t
        + "supp_table_S7.csv; "
        + t
        + "supp_table_S7.tex",
        "Post-hoc relabeling at 2.5/3/4 angstrom, not prospective threshold selection; original cohort RNGs preserved.",
    )
    add(
        "Figure 1",
        "fig:primary_idea",
        "Authored conceptual illustration; no numerical generator",
        p + "figure_index.json; src/posegate/geometry.py; src/posegate/trajectory.py",
        p + "figures/primary_idea.png",
        "Illustrates the measurement/evaluation definitions; not an empirical data plot.",
    )
    add(
        "Figure 2",
        "fig:early-to-late-checkpoints",
        s
        + "make_current_version_figures.py; "
        + s
        + "make_figure3_early_to_late_and_checkpoints.py",
        d
        + "current_statistics_v2/analysis_rows.csv; "
        + d
        + "current_statistics_v2/statistics.json",
        p
        + "figures/current_v2/figure3_early_to_late_and_checkpoints.png; "
        + p
        + "figures/current_v2/checkpoint_plot_data.csv",
        "Source filenames retain earlier figure numbers; manuscript number here is authoritative.",
    )
    add(
        "Figure 3",
        "fig:coord-attrib-scope",
        s
        + "make_current_version_figures.py; "
        + s
        + "make_figure3_physical_attribution_combined.py",
        d
        + "current_statistics_v2/statistics.json; "
        + d
        + "current_statistics_v2/analysis_rows.csv",
        p
        + "figures/current_v2/figure3_physical_attribution.png; "
        + p
        + "figures/current_v2/descriptor_plot_data.csv; "
        + p
        + "figures/current_v2/descriptor_paired_plot_data.csv",
        "Whole-protein geometric descriptors and paired AUROC differences; local controls in S4/S5.",
    )
    add(
        "Figure 4",
        "fig:two-fate-poses",
        p
        + "source_visualization/select_and_extract_frames.py; "
        + p
        + "source_visualization/v3/render_two_fate_pocket_poses_v3.py; "
        + p
        + "source_visualization/v3/compose_two_fate_figure_v3.py",
        p
        + "source_visualization/data/; "
        + p
        + "source_visualization/v3/renders/; "
        + d
        + "current_statistics_v2/analysis_rows.csv; "
        + d
        + "secondary/two_fate_example_prpd_windows.csv",
        p + "figures/two_fate_pose_overlay_figure.png",
        "Historical representative PDB snapshots and metadata; current numeric overlay uses analysis_rows.csv. PyMOL required to rerender, full DCDs required to repeat frame selection.",
    )
    add(
        "Figure 5",
        "fig:preboundary-link-qualification",
        s
        + "make_current_version_figures.py; "
        + s
        + "make_figure4_early_detection_redesign.py",
        d
        + "current_statistics_v2/statistics.json; "
        + d
        + "current_statistics_v2/crossing_times.csv",
        p
        + "figures/current_v2/figure4_early_clean_subsets_and_crossing_timing.png; "
        + d
        + "current_statistics_v2/crossing_summary.json",
        "Expanded-cohort early-clean discrimination and retrospective timing; no-crossing rows have a separate flag.",
    )
    write_csv(ROOT / "MANUSCRIPT_MAPPING.csv", entries)
    text = [
        "# Manuscript mapping\n",
        "Release v0.2.0: six main tables, seven supplementary tables, and five main figures. The supplement contains no figures.\n",
        "The corresponding machine-readable file is [MANUSCRIPT_MAPPING.csv](MANUSCRIPT_MAPPING.csv). Paths are relative to the repository root. Numerical CSVs retain full precision; TeX snapshots preserve published formatting. Definition/protocol tables have exact TeX exports, not invented numerical CSVs.\n",
        "Statistical generators preserve original research paths. Use `python tools/verify_paper_results.py` for portable auditing; raw DCDs are not bundled because of size. Full remeasurement and PyMOL rendering are separate from derived-data auditing.\n",
        "| Item | Label | Scripts | Inputs | Outputs | Audit scope |",
        "|---|---|---|---|---|---|",
    ]
    for entry in entries:
        text.append(
            "| "
            + " | ".join(
                str(entry[k]).replace("|", "\\|").replace("; ", "<br>")
                for k in (
                    "item",
                    "label",
                    "scripts",
                    "inputs",
                    "outputs",
                    "audit_scope",
                )
            )
            + " |"
        )
    text += [
        "\n## Additional controls\n",
        "- Separate-launch results: `source_analysis/summarize_current_separate_launch.py` consumes `data/current_separate_statistics_v2/current_separate_paired_rows.csv` and produces `separate_launch_statistics.json` in that directory. There are 52 target-outcome pairs but 51 complete own-short-horizon pairs. Prefixes here are under `paper/current/`.\n",
        "- Archived radius sensitivity: `source_analysis/35_local_pocket_alignment_control.py` and `36_local_pocket_control_evaluation.py`; inputs/outputs `data/secondary/local_pocket_control_measurements_n81_fixed.csv`, `local_pocket_control_results_n81_fixed.csv`, and `local_pocket_control_theta_agreement.csv`. These 6/8/10 angstrom results predate the current S4 treatment.\n",
        "- Table 5 resource accounting: 9 flagged x (100 - 5) = 855 ns; 8 successful x 95 = 760 ns; 1 false x 95 = 95 ns; 855/2400 = 35.625%, displayed 35.6%. All runs completed. The archived 6 ns summary is preserved, not substituted for these figures.\n",
        "- The primary corrected AUROC interval is the `reported_corrected_auroc` entry of `statistics.json`, tied to the original Python RNG. The other coordinate-control intervals use the recorded NumPy RNG. No interval is selected by its apparent favorability.\n",
        "- Frozen threshold provenance: `paper/current/policy/same_trajectory_5ns_rmsd_threshold_v2.json`, `threshold_selection_v2.json`, and `strict_remeasurement_v2/strict_development_rows_v2.csv`; source selection script `paper/current/source_analysis/55_select_shadow_threshold_v2.py`. The prospective cohort was fixed in `paper/current/cohorts/PROSPECTIVE_CONFIRMATION_V1_MANIFEST.json`.\n",
        "- `paper/current/manifest.json` records original source paths and SHA-256 hashes. `SHA256SUMS` additionally freezes generated CSV/TeX exports, mapping files, configs, fixtures, and code. Selected-coordinate hashes do not cover full DCD bytes.\n",
    ]
    (ROOT / "MANUSCRIPT_MAPPING.md").write_text("\n".join(text) + "\n")
    print(
        f"Exported {len(entries)} mapping entries, 13 TeX tables, and 10 numerical/composition CSVs"
    )


if __name__ == "__main__":
    export_tables()
    export_mapping()
