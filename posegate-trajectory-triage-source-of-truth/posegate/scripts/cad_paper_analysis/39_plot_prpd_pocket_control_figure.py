#!/usr/bin/env python3
"""Main-text figure for the PRPD / local-pocket-alignment control section.

Three panels, values hardcoded from the verified results already reported
in this analysis chain (scripts 33-38 and their evaluation runs):
  A. Component AUROC by checkpoint: RMSD, T, Theta (whole-frame,
     frame-relative), I.
  B. Whole-frame vs pocket-frame (8 A) Theta AUROC by checkpoint.
  C. Compute saved at a matched 10% false-stop budget: RMSD vs
     whole-frame T+Theta vs pocket-frame T+Theta.

Values are not recomputed here -- this script only plots the numbers
already produced and reported by 34_prpd_frozen_lopo_joint_model.py and
36_local_pocket_control_evaluation.py, to keep the figure-generation step
separated from the (expensive, trajectory-reading) analysis step.

Usage:
    python 39_plot_prpd_pocket_control_figure.py --output ../../../doc/YAU_competition_2026/figures/prpd_pocket_control.pdf
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

NAVY = "#173F73"
BLUE = "#3D78B5"
GREEN = "#2E7D4A"
ORANGE = "#D97706"
RED = "#B6403A"
GRAY = "#6B7280"
LIGHT_GRAY = "#E5E7EB"

CHECKPOINTS = [5, 10, 20]

# Panel A: component AUROC (frozen LOPO logistic for RMSD/T/Theta; standalone raw AUROC for I)
RMSD_AUROC = [0.834, 0.805, 0.824]
T_AUROC = [0.779, 0.714, 0.689]
THETA_WHOLE_AUROC = [0.836, 0.800, 0.847]
I_AUROC = [0.449, 0.433, 0.462]  # frozen LOPO single-feature score, matching T/Theta/RMSD methodology

# Panel B: whole- vs pocket-frame Theta AUROC
THETA_POCKET_AUROC = [0.828, 0.814, 0.853]

# Panel C: compute saved at matched 10% false-stop budget
SAVED_RMSD = [0.468, 0.395, 0.384]
SAVED_T_THETA_WHOLE = [0.468, 0.419, 0.395]
SAVED_T_THETA_POCKET = [0.482, 0.382, 0.405]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    x = np.arange(len(CHECKPOINTS))

    ax = axes[0]
    ax.plot(x, RMSD_AUROC, "o-", color=NAVY, label="Corrected RMSD", linewidth=2)
    ax.plot(x, T_AUROC, "s-", color=BLUE, label="T (translation)", linewidth=2)
    ax.plot(x, THETA_WHOLE_AUROC, "^-", color=ORANGE, label=r"$\Theta$ (rotation)", linewidth=2)
    ax.plot(x, I_AUROC, "d-", color=GRAY, label="I (internal deformation)", linewidth=2)
    ax.axhline(0.5, color=LIGHT_GRAY, linewidth=1, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c} ns" for c in CHECKPOINTS])
    ax.set_ylabel("AUROC")
    ax.set_ylim(0.4, 0.95)
    ax.set_title("A. Component AUROC by checkpoint")
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[1]
    width = 0.35
    ax.bar(x - width / 2, THETA_WHOLE_AUROC, width, color=ORANGE, label=r"$\Theta_{\mathrm{whole}}$")
    ax.bar(x + width / 2, THETA_POCKET_AUROC, width, color=GREEN, label=r"$\Theta_{\mathrm{pocket}}$ (8 \AA)")
    ax.axhline(0.5, color=LIGHT_GRAY, linewidth=1, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c} ns" for c in CHECKPOINTS])
    ax.set_ylim(0.4, 0.95)
    ax.set_title(r"B. Whole-frame vs pocket-frame $\Theta$")
    ax.legend(fontsize=8, loc="lower left")

    ax = axes[2]
    width = 0.27
    ax.bar(x - width, [v * 100 for v in SAVED_RMSD], width, color=NAVY, label="RMSD")
    ax.bar(x, [v * 100 for v in SAVED_T_THETA_WHOLE], width, color=ORANGE, label=r"$T{+}\Theta_{\mathrm{whole}}$")
    ax.bar(x + width, [v * 100 for v in SAVED_T_THETA_POCKET], width, color=GREEN, label=r"$T{+}\Theta_{\mathrm{pocket}}$")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{c} ns" for c in CHECKPOINTS])
    ax.set_ylabel("Compute saved (%)")
    ax.set_title("C. Compute saved at 10% false-stop budget")
    ax.legend(fontsize=8, loc="upper right")

    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
