"""The current manuscript's coordinate and nominal-window contract."""
from dataclasses import replace
from pathlib import Path

import MDAnalysis as mda
from MDAnalysis.coordinates.memory import MemoryReader
import numpy as np
import pytest

from posegate.config import load_config, config_from_record, SUPPORTED_PBC_METHOD
from posegate.exceptions import CheckpointNotReady, InputValidationError, ConfigurationError
from posegate.geometry import protein_relative_coordinates, validate_raw_ligand_wholeness
from posegate.policy import evaluate_policy
from posegate.trajectory import measure_prefix, measure_late_outcome, nominal_window_mask

ROOT = Path(__file__).resolve().parents[1]


def current_config():
    return load_config(ROOT / "configs/posegate_5ns_v2.yaml")


def test_current_config_matches_the_frozen_threshold_and_round_trips():
    config = current_config()
    assert config.timing.frame_interval_ns == 0.2
    assert config.policy.stop_threshold_angstrom == 2.7634173197224974
    assert config.policy.type == "threshold_rule"
    assert config_from_record(config.scientific_dict()).scientific_dict() == config.scientific_dict()


def test_common_lattice_shift_preserves_geometry_at_an_image_boundary():
    protein = np.array([[-1.,0,0], [1.,0,0], [0,1.,0], [0,-1.,0]])
    ligand = np.array([[4.8,0,0], [5.2,0,0], [5.,0.3,0]])
    box = np.array([10.,10,10,90,90,90])
    _, whole = protein_relative_coordinates(protein, ligand, box)
    _, legacy = protein_relative_coordinates(protein, ligand, box, method=SUPPORTED_PBC_METHOD)
    assert np.linalg.norm(whole[0]-whole[1]) == pytest.approx(0.4)
    assert np.linalg.norm(legacy[0]-legacy[1]) == pytest.approx(9.6)
    assert np.allclose(whole-ligand, (whole-ligand)[0])


def test_common_shift_in_a_triclinic_box():
    from MDAnalysis.lib.mdamath import triclinic_vectors
    protein = np.array([[-1.,0,0], [1.,0,0], [0,1.,0], [0,-1.,0]])
    ligand = np.array([[1.,1,1], [2.,1,1], [1.,2,1]])
    box = np.array([30.,32,34,75,80,65])
    lattice = triclinic_vectors(box).astype(float)
    _, baseline = protein_relative_coordinates(protein, ligand, box)
    _, shifted = protein_relative_coordinates(protein, ligand + lattice[0] - lattice[2], box)
    assert np.allclose(baseline, shifted, atol=1e-5)


def test_torn_raw_ligand_is_rejected_before_imaging():
    with pytest.raises(InputValidationError, match="whole coordinates"):
        validate_raw_ligand_wholeness(np.array([[1.,0,0], [99.,0,0]]),
                                     np.array([100.,100,100,90,90,90]))


@pytest.mark.parametrize("dt,start,end,first,last,count", [
    (0.2,3,5,15,24,10), (0.02,3,5,150,249,100),
    (0.2,70,100,350,499,150), (0.02,14,20,700,999,300),
])
def test_nominal_windows_match_supplementary_table_s1(dt,start,end,first,last,count):
    selected = np.flatnonzero(nominal_window_mask(round(end/dt),start,end,dt))
    assert (selected[0],selected[-1],len(selected)) == (first,last,count)


def test_incomplete_and_off_grid_windows_fail():
    with pytest.raises(CheckpointNotReady):
        nominal_window_mask(24,3,5,0.2)
    with pytest.raises(InputValidationError):
        nominal_window_mask(30,3.1,5,0.2)


def universe(frames=500, dt=200.0000006):
    protein = np.array([[40,40,40],[45,40,40],[40,45,40],[40,40,45]], dtype=float)
    ligand = np.array([[42,42,42],[43,42,42],[42,43,42],[42,42,43]], dtype=float)
    coords = np.repeat(np.vstack([protein,ligand])[None,:,:],frames,axis=0)
    coords[1:,4:,:] += np.array([3.,0,0])
    coords[1:,4:,0] += 100
    u = mda.Universe.empty(8,n_residues=2,atom_resindex=[0]*4+[1]*4)
    u.add_TopologyAttr("names", ["CA"]*4+["C1","C2","C3","C4"])
    u.add_TopologyAttr("resnames", ["ALA","UNK"])
    u.load_new(coords,format=MemoryReader,dt=dt,dimensions=np.tile([100,100,100,90,90,90],(frames,1)))
    return u


def test_current_prefix_controls_and_late_boundary(monkeypatch):
    u = universe()
    monkeypatch.setattr("posegate.trajectory._universe",lambda *_:u)
    config = current_config()
    prefix = measure_prefix("mock.pdb","mock.dcd",config)
    assert prefix.window_frame_count == 10
    assert prefix.corrected_pose_rmsd_mean_angstrom == pytest.approx(3.,abs=1e-5)
    assert prefix.pbc_naive_pose_rmsd_mean_angstrom == pytest.approx(103.,abs=1e-5)
    assert prefix.ligand_self_aligned_rmsd_mean_angstrom == pytest.approx(0.,abs=1e-5)
    outcome = measure_late_outcome("mock.pdb","mock.dcd",config)
    assert outcome.late_window_frame_count == 150
    assert not outcome.pose_retained


def test_prefix_never_reads_a_future_frame(monkeypatch):
    u = universe(frames=100)
    original = type(u.trajectory).__getitem__
    accessed = []
    def guarded(reader,index):
        assert index < 25
        accessed.append(index)
        return original(reader,index)
    monkeypatch.setattr(type(u.trajectory),"__getitem__",guarded)
    monkeypatch.setattr("posegate.trajectory._universe",lambda *_:u)
    measure_prefix("mock.pdb","mock.dcd",current_config())
    assert accessed == list(range(25))


def test_nominal_cadence_is_checked(monkeypatch):
    monkeypatch.setattr("posegate.trajectory._universe",lambda *_:universe(dt=20))
    with pytest.raises(InputValidationError,match="cadence"):
        measure_prefix("mock.pdb","mock.dcd",current_config())


def test_nonfinite_values_cannot_trigger_or_pass_a_policy():
    for value in [np.nan, np.inf, -np.inf]:
        with pytest.raises(ConfigurationError,match="finite"):
            evaluate_policy(current_config().policy,{"corrected_pose_rmsd_mean_angstrom":value})
