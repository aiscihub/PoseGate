"""Small saved-frame fixtures from the current reanalysis, not new simulations."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from posegate.geometry import (
    frame_relative_geometry,
    protein_relative_coordinates,
    whole_ligand_relative_to_protein,
    zero_geometry,
)

DATA = Path(__file__).parent / "data/current_reference_frames"
FIXTURES = json.loads((DATA / "manifest.json").read_text())


@pytest.mark.parametrize("entry", FIXTURES, ids=lambda e: e["complex_id"])
def test_current_saved_frame_measurements(entry):
    path = DATA / entry["file"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
    with np.load(path, allow_pickle=False) as data:
        reference_protein, reference_ligand = protein_relative_coordinates(
            data["protein"][0], data["ligand"][0], data["box"][0]
        )
        reference_raw = data["ligand"][0].astype(float) - data["protein"][0].astype(
            float
        ).mean(axis=0)
        measured = []
        for i, (protein, ligand, box) in enumerate(
            zip(data["protein"], data["ligand"], data["box"], strict=True)
        ):
            centered, relative = protein_relative_coordinates(protein, ligand, box)
            center = protein.astype(float).mean(axis=0)
            _, disagreement = whole_ligand_relative_to_protein(ligand, center, box)
            geometry = (
                zero_geometry()
                if i == 0
                else frame_relative_geometry(
                    centered,
                    reference_protein,
                    relative,
                    reference_ligand,
                    raw_ligand_relative=ligand.astype(float) - center,
                    reference_raw_ligand_relative=reference_raw,
                )
            )
            measured.append(
                [
                    geometry.corrected_pose_rmsd_angstrom,
                    geometry.pbc_naive_pose_rmsd_angstrom,
                    geometry.internal_deformation_angstrom,
                    geometry.centroid_displacement_angstrom,
                    np.degrees(geometry.frame_relative_reorientation_radians),
                    disagreement,
                ]
            )
        # Different supported MDAnalysis/NumPy versions need not be bitwise equal.
        np.testing.assert_allclose(measured, data["expected"], rtol=0, atol=2e-5)
        early = (data["frame_indices"] >= 15) & (data["frame_indices"] < 25)
        assert early.sum() == 10
        np.testing.assert_allclose(
            np.asarray(measured)[early, 0].mean(),
            data["expected_x5"],
            rtol=0,
            atol=2e-5,
        )
