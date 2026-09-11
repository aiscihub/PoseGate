#!/usr/bin/env python3
"""Check the exported paper tables against current numerical inputs, without MD."""

from __future__ import annotations

import hashlib
from decimal import Decimal, ROUND_HALF_UP
import importlib.util
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper/current"


def module(name):
    path = PAPER / "source_analysis" / (name + ".py")
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def auc(y, x):
    labels, values = np.asarray(y, dtype=int), np.asarray(x, dtype=float)
    n1, n0 = labels.sum(), len(labels) - labels.sum()
    if n1 == 0 or n0 == 0:
        return None
    return float((rankdata(values)[labels == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def close(a, b):
    assert np.isclose(a, b, rtol=0, atol=1e-10), (a, b)


def display(value):
    """Match manuscript half-up rounding, including the exact AUROC 0.8125."""
    return str(Decimal(str(value)).quantize(Decimal("0.001"), rounding=ROUND_HALF_UP))


def main():
    manifest = json.loads((PAPER / "manifest.json").read_text())
    for name, entry in manifest["files"].items():
        assert (
            hashlib.sha256((PAPER / name).read_bytes()).hexdigest() == entry["sha256"]
        ), name
    tables = json.loads((PAPER / "table_index.json").read_text())
    assert len(tables) == 13

    def contains(label, values):
        text = tables[label]["latex"]
        for value in values:
            assert value in text, (label, value)

    base = PAPER / "data"
    rows = pd.read_csv(base / "current_statistics_v2/analysis_rows.csv")
    stats = json.loads((base / "current_statistics_v2/statistics.json").read_text())
    primary = rows.loc[rows.primary]
    assert (
        len(primary),
        primary.protein.nunique(),
        len(rows),
        rows.protein.nunique(),
    ) == (52, 15, 81, 16)
    cohort = json.loads((PAPER / "cohorts/final_expanded_81.json").read_text())["cases"]
    assert len(cohort) == 81 and sum(c["primary"] for c in cohort) == 52
    assert {c["complex_id"] for c in cohort} == set(rows.complex_id)
    indexed = rows.set_index("complex_id")
    for case in cohort:
        row = indexed.loc[case["complex_id"]]
        close(case["X5_A"], row.corrected_A)
        close(case["L100_A"], row.late_corrected_A)
        trace = pd.read_csv(PAPER / case["frames"])
        assert trace.frame_index.tolist() == list(range(500))
        close(trace.corrected_A.iloc[15:25].mean(), row.corrected_A)
        close(trace.corrected_A.iloc[350:500].median(), row.late_corrected_A)
    observed = {}
    for cohort, data in (
        ("Primary milbemycin", primary),
        ("Expanded milbemycin", rows),
    ):
        key = "primary" if cohort.startswith("Primary") else "expanded"
        observed[cohort] = {}
        for score in ("corrected_A", "pbc_naive_A", "ligand_self_A"):
            value = auc(data.nonretained, data[score])
            close(value, stats[key]["metrics"][score]["estimate"])
            observed[cohort][score] = value
    recomputed = module("summarize_current_frames").clustered_scores(
        primary, ["pbc_naive_A", "ligand_self_A"]
    )
    for score, result in recomputed["metrics"].items():
        for field in ("estimate", "ci_low", "ci_high"):
            close(result[field], stats["primary"]["metrics"][score][field])

    for score in ("corrected_A", "pbc_naive_A", "ligand_self_A"):
        r = stats["correlations"][score]
        close(
            spearmanr(primary[score], primary["late_" + score]).statistic,
            r["pooled"]["estimate"],
        )
        x, y = rankdata(primary[score]), rankdata(primary["late_" + score])
        for protein in primary.protein.unique():
            mask = primary.protein.to_numpy() == protein
            x[mask] -= x[mask].mean()
            y[mask] -= y[mask].mean()
        close(np.corrcoef(x, y)[0, 1], r["within"]["estimate"])
        contains(
            "tab:representation-matched",
            [
                f"{v:.3f}"
                for s in ("pooled", "within")
                for v in (r[s]["estimate"], r[s]["ci_low"], r[s]["ci_high"])
            ],
        )

    for name, expected in (
        ("early_low", 46),
        ("rolling_prefix_clean", 45),
        ("framewise_clean", 31),
    ):
        assert int(rows[name].sum()) == expected
    assert int(primary.framewise_clean.sum()) == 20
    temporal = json.loads(
        (
            base / "final_protein_frame_control/primary_temporal_statistics.json"
        ).read_text()
    )
    for data, subsets in ((rows, stats["subsets"]), (primary, temporal["subsets"])):
        for name, result in subsets.items():
            subset = data.loc[data[name]]
            for score in ("corrected_A", "reorientation_deg"):
                metric = result["metrics"][score]
                assert len(subset) == metric["n"]
                close(auc(subset.nonretained, subset[score]), metric["estimate"])
                contains(
                    "tab:prefix-clean-risksets",
                    [display(metric[k]) for k in ("estimate", "ci_low", "ci_high")],
                )
    events = pd.read_csv(base / "current_statistics_v2/crossing_times.csv")
    times = events.loc[events.complex_id.isin(primary.complex_id) & events.nonretained]
    assert len(times) == 39
    assert int((times.event_observed & (times.t_star_ns <= 5)).sum()) == 18
    after = times.loc[times.event_observed & (times.t_star_ns > 5), "t_star_ns"] - 5
    assert len(after) == 20 and int((~times.event_observed).sum()) == 1
    np.testing.assert_allclose(np.percentile(after, [25, 50, 75]), [11.5, 53, 68.5])
    contains("tab:time-to-event", ["18 (46", "20 (51", "53 [11.5, 68.5]", "1 (3"])
    for name, result in stats["subsets"].items():
        for score in ("corrected_A", "reorientation_deg"):
            close(
                auc(rows.loc[rows[name], "nonretained"], rows.loc[rows[name], score]),
                result["metrics"][score]["estimate"],
            )
    for checkpoint, values in stats["checkpoints"].items():
        for cohort, data in (("primary", primary), ("expanded", rows)):
            close(
                auc(data.nonretained, data["corrected_A_" + checkpoint]),
                values[cohort]["metrics"]["corrected_A"]["estimate"],
            )

    separate = pd.read_csv(
        base / "current_separate_statistics_v2/current_separate_paired_rows.csv"
    )
    separate_stats = json.loads(
        (
            base / "current_separate_statistics_v2/separate_launch_statistics.json"
        ).read_text()
    )["statistics"]
    assert len(separate) == 52 and separate.protein.nunique() == 15
    for score, key in (
        ("primary_x5", "same_trajectory_x5_vs_original_y100"),
        ("separate_x5", "separate_launch_x5_vs_original_y100"),
    ):
        close(
            auc(separate.primary_nonretained_100, separate[score]),
            separate_stats[key]["estimate"],
        )
    complete = separate.loc[separate.separate_late_complete]
    assert len(complete) == 51
    for score, label, key in (
        ("primary_x5", "primary_nonretained_20", "original_x5_vs_own_y20"),
        ("separate_x5", "separate_nonretained_20", "independent_x5_vs_own_y20"),
    ):
        close(auc(complete[label], complete[score]), separate_stats[key]["estimate"])

    validation = None
    for cohort, folder, late in (
        ("Prospective validation", "prospective_validation", "70_100"),
        ("Prospective pilot", "prospective_pilot", "14_20"),
    ):
        data = pd.read_csv(base / folder / "measurements.csv")
        metrics = pd.read_csv(base / folder / "metrics.csv")
        assert (
            len(data) == 24
            and data.protein.nunique() == 16
            and (data.status == "ok").all()
        )
        observed[cohort] = {}
        for result in metrics.itertuples():
            subset = data
            if result.metric_id.endswith("_strict_rolling_prefix_clean"):
                subset = data.loc[data.strict_rolling_prefix_clean_5ns]
            elif result.metric_id.endswith("_framewise_prefix_clean"):
                subset = data.loc[data.framewise_prefix_clean_5ns]
            assert len(subset) == result.n
            value = auc(subset["late_nonretained_" + late], subset[result.score_column])
            close(value, result.auroc)
            if result.metric_id in ("corrected_X5", "pbc_naive_X5", "ligand_self_X5"):
                observed[cohort][result.metric_id] = value
        if folder == "prospective_validation":
            validation = data
    y = validation.late_nonretained_70_100.to_numpy()
    flagged = validation.corrected_rmsd_X5_mean_3_5_A.to_numpy() > 2.7634173197224974
    assert (
        int((~y).sum()),
        int(y.sum()),
        int(flagged.sum()),
        int((flagged & y).sum()),
        int((flagged & ~y).sum()),
    ) == (9, 15, 9, 8, 1)
    assert int(flagged.sum()) * (100 - 5) == 855
    contains(
        "tab:prospective-validation", ["855", "760", "95", "35.6", "0.748", "0.658"]
    )
    vstats = module("77_analyze_confirmation_v1_fresh100ns_preliminary")
    interval = vstats.auc_summary(
        validation,
        "corrected_rmsd_X5_mean_3_5_A",
        n_bootstrap=5000,
        seed=vstats.BOOTSTRAP_SEED,
    )
    contains(
        "tab:prospective-validation",
        [f"{interval[v]:.3f}" for v in ("auroc", "ci_lo", "ci_hi")],
    )

    cross = pd.read_csv(base / "current_statistics_v2/crossligand_rows.csv")
    for name, count in (("beauvericin", 29), ("verapamil", 45)):
        data = cross.loc[cross.cohort == name]
        assert len(data) == count
        observed[name.title()] = {
            s: auc(data.nonretained, data[s])
            for s in ("corrected_A", "pbc_naive_A", "ligand_self_A")
        }
    for name, values in observed.items():
        if name != "Expanded milbemycin":
            contains(
                "tab:cross_ligand_results",
                [name] + [f"{v:.3f}" for v in values.values()],
            )

    final = pd.read_csv(
        base / "final_protein_frame_control/protein_side_stability_summary.csv"
    )
    stability = json.loads(
        (
            base / "final_protein_frame_control/protein_side_stability_statistics.json"
        ).read_text()
    )
    for key, result in stability["metrics"].items():
        close(final[key].median(), result["median"])
        close(auc(final.nonretained, final[key]), result["auroc_for_nonretained"])
        close(
            spearmanr(final[key], final.late_corrected_A).statistic,
            result["spearman_vs_late_corrected_A"],
        )
        if key != "protein_ca_rmsd_mean_8_10_A":
            contains(
                "tab:protein-frame-stability",
                [
                    display(result[k])
                    for k in (
                        "median",
                        "auroc_for_nonretained",
                        "spearman_vs_late_corrected_A",
                    )
                ]
                + [
                    display(result[k][v])
                    for k in ("spearman_ci", "auroc_ci")
                    for v in ("ci_low", "ci_high")
                ],
            )
    local = rows.merge(
        pd.read_csv(base / "current_statistics_v2/local_pocket_rows.csv"),
        on=["complex_id", "protein"],
        validate="one_to_one",
    )
    local_stats = json.loads(
        (base / "current_statistics_v2/local_pocket_statistics.json").read_text()
    )
    for checkpoint, result in local_stats["checkpoints"].items():
        whole, pocket = (
            local["reorientation_deg_" + checkpoint],
            local["theta_" + checkpoint],
        )
        close(auc(local.nonretained, whole), result["whole"]["estimate"])
        close(auc(local.nonretained, pocket), result["local"]["estimate"])
        close(spearmanr(whole, pocket).statistic, result["spearman"])
        contains(
            "tab:pocket-frame-control",
            [
                display(result[frame][field])
                for frame in ("whole", "local")
                for field in ("estimate", "ci_low", "ci_high")
            ]
            + [
                display(result[k])
                for k in ("delta_local_minus_whole", "spearman", "pearson")
            ],
        )
    boundary = pd.read_csv(
        base / "current_boundary_sensitivity/boundary_sensitivity.csv"
    )
    for row in boundary.itertuples():
        data = primary if row.cohort == "Primary milbemycin" else validation
        x = (
            data.corrected_A
            if row.cohort == "Primary milbemycin"
            else data.corrected_rmsd_X5_mean_3_5_A
        )
        late = (
            data.late_corrected_A
            if row.cohort == "Primary milbemycin"
            else data.late_corrected_rmsd_median_70_100_A
        )
        close(auc(late >= row.boundary_A, x), row.estimate)
        assert int((late < row.boundary_A).sum()) == row.retained
        contains(
            "tab:boundary-sensitivity",
            [f"{row.estimate:.3f}", f"{row.ci_low:.3f}", f"{row.ci_high:.3f}"],
        )
    contact = pd.read_csv(base / "secondary/contact_retention.csv")
    master = pd.read_csv(
        base / "secondary/master/expanded_trajectory_source_of_truth.csv"
    )
    extension = pd.read_csv(base / "secondary/consolidated_100ns_extension_n81.csv")
    contact_inputs = master.loc[master.eligible_primary_cohort].merge(
        extension, on="complex_id", suffixes=("", "_ext"), validate="one_to_one"
    )
    for row in contact.itertuples():
        pairs = contact_inputs.dropna(subset=[row.predictor, row.target])
        assert len(pairs) == row.n
        close(
            spearmanr(pairs[row.predictor], pairs[row.target]).statistic,
            row.spearman_rho,
        )
        contains(
            "tab:contact-retention-control",
            [
                f"{row.spearman_rho:.3f}",
                f"{row.spearman_ci_lo:.3f}",
                f"{row.spearman_ci_hi:.3f}",
            ],
        )
    print(
        json.dumps(
            {
                "coordinate_auroc": observed,
                "validation_flagged": 9,
                "potential_remaining_ns": 855,
            },
            indent=2,
        )
    )
    print(
        "PASS: file hashes, table snapshots, current point estimates, primary control and validation AUROC intervals."
    )
    print(
        "Other reported intervals are preserved source estimates; this check does not rerun raw trajectories."
    )


if __name__ == "__main__":
    main()
