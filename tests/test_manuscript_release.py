"""Release wiring and policy parity against distributed manuscript inputs."""

import csv
import json
from pathlib import Path
import re

import numpy as np

from posegate.config import load_config
from posegate.policy import evaluate_policy

ROOT = Path(__file__).resolve().parents[1]
PAPER = ROOT / "paper/current"


def test_table5_uses_the_bundled_one_sided_rule():
    config = load_config(ROOT / "configs/posegate_5ns_v2.yaml")
    frozen = json.loads((ROOT / config.provenance["source_policy"]).read_text())
    assert (
        frozen["decision"]["stop_threshold_A"] == config.policy.stop_threshold_angstrom
    )
    assert config.policy.stop_threshold_angstrom == 2.7634173197224974
    assert config.policy.stop_direction == "greater_than"
    assert config.policy.tie_rule == "no_stop"
    with (PAPER / "data/prospective_validation/measurements.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    flagged = []
    for row in rows:
        decision = evaluate_policy(
            config.policy,
            {
                "corrected_pose_rmsd_mean_angstrom": float(
                    row["corrected_rmsd_X5_mean_3_5_A"]
                )
            },
        )
        assert decision.uncalibrated_retention_score is None
        if decision.forecast == "POSE_NONRETAINED":
            flagged.append(row)
        else:
            assert decision.forecast == "DEFER"
    assert len(rows) == 24 and len(flagged) == 9
    assert sum(row["late_nonretained_70_100"].lower() == "true" for row in flagged) == 8
    threshold = config.policy.stop_threshold_angstrom
    for value, expected in (
        (np.nextafter(threshold, -np.inf), "DEFER"),
        (threshold, "DEFER"),
        (np.nextafter(threshold, np.inf), "POSE_NONRETAINED"),
    ):
        assert (
            evaluate_policy(
                config.policy, {"corrected_pose_rmsd_mean_angstrom": value}
            ).forecast
            == expected
        )


def test_mapping_covers_every_table_and_figure_and_links_exist():
    for stem in ("yau_aj", "yau_supplementary_materials"):
        assert not list(PAPER.glob(stem + ".*")), "Unsubmitted manuscript must not be distributed"
    with (ROOT / "MANUSCRIPT_MAPPING.csv").open() as handle:
        mapping = list(csv.DictReader(handle))
    tables = json.loads((PAPER / "table_index.json").read_text())
    figures = set(json.loads((PAPER / "figure_index.json").read_text()))
    assert len(figures) == 5
    assert len(mapping) == 18
    assert {row["label"] for row in mapping} == set(tables) | figures
    markdown = (ROOT / "MANUSCRIPT_MAPPING.md").read_text()
    for row in mapping:
        assert row["label"] in markdown
        for kind in ("scripts", "inputs", "outputs"):
            for value in row[kind].split("; "):
                if value.startswith(("paper/", "tools/", "src/", "configs/")):
                    assert (ROOT / value).exists(), value
    assert len(list((PAPER / "tables").glob("*.tex"))) == 13
    assert len(list((PAPER / "tables").glob("*.csv"))) == 10
    assert {p.name for p in (ROOT / "configs").glob("*.yaml")} == {
        "posegate_5ns_v2.yaml",
        "posegate_pilot20ns_v2.yaml",
    }
    for obsolete in (
        "assets",
        "posegate-trajectory-triage-source-of-truth",
        "shadow-test-archive",
    ):
        assert not (ROOT / obsolete).exists()


def test_extracted_latex_tables_match_the_frozen_sources():
    tables = json.loads((PAPER / "table_index.json").read_text())
    counters = {"yau_aj.tex": 0, "yau_supplementary_materials.tex": 0}
    for entry in tables.values():
        counters[entry["source"]] += 1
        prefix = "main_table_" if entry["source"] == "yau_aj.tex" else "supp_table_S"
        exported = PAPER / "tables" / f"{prefix}{counters[entry['source']]}.tex"
        assert exported.read_text() == entry["latex"] + "\n"
