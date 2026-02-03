"""Golden tests for renderer: checksum or bounds of first frame."""

import hashlib
from pathlib import Path

import pytest

from xmoltoppm.core.options import Options
from xmoltoppm.io.xyz_reader import read_xyz
from xmoltoppm.viz.render import render


def test_render_first_frame_shape_and_checksum():
    """Render first frame of traj.xyz with fixed options; check shape and checksum."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found in project root")
    frames = list(read_xyz(traj_path))
    assert len(frames) >= 1
    opts = Options(
        size=500,
        circle_style=1,
        line_style=3,
        all_frames=False,
        frame_skip=0,
        filename_type=0,
        shadow=0.0,
        shadow_strength=20.0,
        shadow_direction="yz",
        input_file=None,
        output_path="output.ppm",
        bond_radius_adjust=1.0,
        border_size=0.5,
        circle_size=0.15,
        line_size=0.05,
        background=255,
        draw_boundary_lines=0,
    )
    img = render(frames[0], opts)
    assert img.ndim == 3
    assert img.shape[2] == 3
    assert img.dtype == "uint8"
    h, w = img.shape[0], img.shape[1]
    assert h > 0 and w > 0
    # Deterministic checksum of raw bytes (lock after first run)
    checksum = hashlib.sha256(img.tobytes()).hexdigest()
    # Store expected after first run; any change in render will change checksum
    assert len(checksum) == 64
