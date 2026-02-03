"""Golden tests for XYZ reader."""

from pathlib import Path

import pytest

from xmoltoppm.io.xyz_reader import read_xyz


def test_traj_xyz_first_frame():
    """First frame of traj.xyz has expected nat, molname, iteration, energy, cell, first atom."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found in project root")
    frames = list(read_xyz(traj_path))
    assert len(frames) >= 1
    f = frames[0]
    assert f.nat == 6
    assert f.molname.strip() == "Acetylene_O2"
    assert f.iteration == 0
    assert abs(f.energy - (-605.26)) < 0.01
    assert f.symbols[0] == "H"
    assert abs(f.coords[0, 0] - 38.45539) < 1e-5
    assert abs(f.coords[0, 1] - 41.24817) < 1e-5
    assert abs(f.coords[0, 2] - 38.89843) < 1e-5
    assert f.cell is not None
    assert f.cell.shape == (3, 3)
    # Orthorhombic 80,80,80 -> diagonal-ish
    assert abs(f.cell[0, 0] - 80.0) < 1.0
    assert abs(f.cell[2, 2] - 80.0) < 1.0


def test_minimal_xyz_two_frames():
    """Minimal 2-atom 2-frame XYZ parses correctly."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "minimal.xyz"
    frames = list(read_xyz(path))
    assert len(frames) == 2
    assert frames[0].nat == 2
    assert frames[0].symbols == ["H", "H"]
    assert frames[0].iteration == 0
    assert frames[1].iteration == 1
    assert abs(frames[0].coords[1, 0] - 1.0) < 1e-9
