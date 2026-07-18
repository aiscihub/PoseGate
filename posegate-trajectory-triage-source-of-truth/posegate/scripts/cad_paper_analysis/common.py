"""Shared paths and constants for the Cage-Anchor Dynamics (CAD) paper analysis.

This is the single place to edit when reproducing this analysis on a new
dataset (new ligand, new protein panel, new simulation batch). Everything
downstream reads from these paths.
"""

from pathlib import Path

# Repo root (this file lives at yau_competition/scripts/cad_paper_analysis/common.py)
REPO_ROOT = Path(__file__).resolve().parents[3]

# Root of the raw 100 ns explicit-solvent trajectories, one subdir per protein:
#   BASE_100NS / <protein> / simulation_explicit / <pocket> / <replica> /
#       <protein>_prepared_<ligand>_<pocket>_complex_recombined_complex_explicit_initial_frame.pdb
#       <protein>_prepared_<ligand>_<pocket>_complex_recombined_complex_explicit_trajectory.dcd
#       plip_results/all_plip_interactions_summary.csv
# Point this at a new dataset's simulation directory to extend the analysis.
BASE_100NS = Path("/media/zhenli/datadrive/valleyfevermutation/simulation_100ns_md")

# Root of the genuinely independent 20 ns explicit-solvent trajectories (a separate
# simulation batch from BASE_100NS, not a truncation of it -- same directory layout).
# Used by 16_pocket_aware_state_20ns.py so the 20 ns pocket-aware state is computed
# without any dependence on the 100 ns run it will later be validated against.
BASE_20NS = Path("/media/zhenli/datadrive/valleyfevermutation/simulation_20ns_md")

# Cached feature/label tables produced by pipeline/training/yau_recapture_screen.py
SCREEN_ROOT = REPO_ROOT / "pipeline" / "training" / "outputs" / "yau_recapture_screen"
SCREENING_TABLE = SCREEN_ROOT / "yau_recapture_screening_table.csv"
RBE_TRACES = SCREEN_ROOT / "rbe_traces.csv"  # per-frame anchor_retention + (uncorrected) cage_distance
EARLY_WINDOW_TABLE = SCREEN_ROOT / "early_window_feature_summary_by_complex.csv"

# Output directory for this analysis (created if missing)
OUT_DIR = SCREEN_ROOT / "cad_states"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Corrected per-frame trace file produced by 01_recompute_pbc_corrected_traces.py
PBC_TRACES = OUT_DIR / "rbe_traces_pbc_corrected.csv"

# Figure output (referenced by doc/YAU_competition_2026/yau_aj.tex)
FIGURE_OUT = REPO_ROOT / "doc" / "YAU_competition_2026" / "figures" / "state_transition_result.png"

# Ligand / contact-detection constants (must match pipeline/training/yau_recapture_screen.py)
LIGAND_RESNAME = "UNK"
CONTACT_CUTOFF_A = 4.5
ANCHOR_BASELINE_FREQ = 0.20
CAGE_BUFFER_A = 2.0

# Cage-Anchor Dynamics state thresholds (Eqs. cage_indicator / anchor_indicator in the paper)
TAU_C = 1.0   # C(t) <= TAU_C  => geometrically caged
TAU_R = 0.5   # R(t) >= TAU_R  => anchored

TRANSITION_CLASSES = ["persistent_stable", "persistent_unstable", "delayed_failure", "recovery"]
