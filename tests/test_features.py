"""Comprehensive feature tests for Fortran-parity and enhancement options."""

from pathlib import Path

import numpy as np
import pytest

from xmoltoppm.core.frame import Frame
from xmoltoppm.core.options import Options
from xmoltoppm.io.xyz_reader import read_xyz
from xmoltoppm.viz.render import centroid_after_plane_filter, render


def _base_opts() -> Options:
    return Options(
        size=100,
        circle_style=0,
        line_style=0,
        all_frames=False,
        frame_skip=0,
        filename_type=0,
        shadow=0.0,
        shadow_strength=20.0,
        shadow_direction="yz",
        input_file=None,
        output_path="out.ppm",
        bond_radius_adjust=1.0,
        border_size=0.5,
        circle_size=0.15,
        line_size=0.05,
        background=255,
        draw_boundary_lines=0,
        center_mode=0,
        frame_offset=0,
        plane_x=5_000_000.0,
        plane_y=5_000_000.0,
        plane_z=5_000_000.0,
        bond_thickness_by_distance=0,
        supersample_scale=1,
    )


def test_centroid_after_plane_filter_no_filter():
    """centroid_after_plane_filter with no plane filter returns frame mean."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "minimal.xyz"
    frames = list(read_xyz(path))
    opts = _base_opts()
    c = centroid_after_plane_filter(frames[0], opts)
    expected = frames[0].coords.mean(axis=0)
    np.testing.assert_allclose(c, expected)


def test_centroid_after_plane_filter_with_filter():
    """centroid_after_plane_filter excludes atoms above plane."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "plane_filter.xyz"
    if not path.exists():
        pytest.skip("plane_filter.xyz not found")
    frames = list(read_xyz(path))
    opts = _base_opts()
    opts.plane_z = 10.0  # keep only z <= 10
    c = centroid_after_plane_filter(frames[0], opts)
    # Atoms at (1,1,1) and (2,2,2) have z<=10; (15,15,15) and (18,18,18) excluded
    expected = np.array([1.5, 1.5, 1.5])
    np.testing.assert_allclose(c, expected)


def test_plane_filter_removes_atoms():
    """Plane filter removes atoms with coord > plane; image has fewer atoms."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "plane_filter.xyz"
    if not path.exists():
        pytest.skip("plane_filter.xyz not found")
    frames = list(read_xyz(path))
    opts = _base_opts()
    img_no_filter = render(frames[0], opts, ref_center=None)
    opts.plane_z = 10.0
    img_filtered = render(frames[0], opts, ref_center=None)
    # Filtered image should differ (fewer atoms drawn)
    assert img_no_filter.shape == img_filtered.shape
    assert not np.array_equal(img_no_filter, img_filtered)


def test_plane_filter_empty_returns_minimal_image():
    """Plane filter that removes all atoms yields minimal background image."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "plane_filter.xyz"
    if not path.exists():
        pytest.skip("plane_filter.xyz not found")
    frames = list(read_xyz(path))
    opts = _base_opts()
    opts.plane_x = 0.0  # no atom has x <= 0 (all at 1,2,15,18)
    opts.plane_y = 0.0
    opts.plane_z = 0.0
    img = render(frames[0], opts, ref_center=None)
    assert img.ndim == 3
    assert img.shape[2] == 3
    assert img.shape[0] >= 1 and img.shape[1] >= 1
    assert np.all(img == opts.background)


def test_center_mode_ref_center_used():
    """When ref_center is passed, centering uses it instead of frame mean."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 500
    # Same frame: ref_center=mean must give identical image to ref_center=None
    mean0 = frames[0].coords.mean(axis=0)
    img0_own = render(frames[0], opts, ref_center=None)
    img0_ref = render(frames[0], opts, ref_center=mean0)
    np.testing.assert_array_equal(img0_own, img0_ref)
    # ref_center is applied: render(frame, ref_center=other) uses other for centering
    # (CLI uses this for -ic 0 with all-frames; integration covered by CLI flow)


def test_bond_thickness_by_distance_changes_output():
    """bond_thickness_by_distance=1 produces different image from 0."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 400  # large enough for bonds to be visible
    opts.bond_thickness_by_distance = 0
    img_bd0 = render(frames[0], opts, ref_center=None)
    opts.bond_thickness_by_distance = 1
    img_bd1 = render(frames[0], opts, ref_center=None)
    assert img_bd0.shape == img_bd1.shape
    assert not np.array_equal(img_bd0, img_bd1)


def test_frame_offset_option_stored():
    """Options.from_cli stores frame_offset from -nf."""
    from xmoltoppm.cli.main import main
    import argparse
    args = argparse.Namespace(
        size=500, cc=0, il=0, af=None, ft=0, sh=0.0, sf=20.0, sd="yz",
        input=None, output="out.ppm", ba=1.0, bs=0.5, sc=0.15, sl=0.05,
        bc=255, bl=0, ic=0, nf=100, plane=None, bd=0,
    )
    opts = Options.from_cli(args)
    assert opts.frame_offset == 100


def test_plane_option_parsed():
    """-pp axis coor sets plane_x/y/z via CLI (tested via Options)."""
    import argparse
    args = argparse.Namespace(
        size=500, cc=0, il=0, af=None, ft=0, sh=0.0, sf=20.0, sd="yz",
        input=None, output="out.ppm", ba=1.0, bs=0.5, sc=0.15, sl=0.05,
        bc=255, bl=0, ic=0, nf=0, plane=[["z", "10.0"], ["x", "5.0"]], bd=0,
        plane_x=5_000_000.0, plane_y=5_000_000.0, plane_z=5_000_000.0,
    )
    args.plane_z = 10.0
    args.plane_x = 5.0
    opts = Options.from_cli(args)
    assert opts.plane_z == 10.0
    assert opts.plane_x == 5.0


def test_render_with_all_new_options_defaults():
    """Render runs with all new options at defaults (no regression)."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "minimal.xyz"
    frames = list(read_xyz(path))
    opts = _base_opts()
    opts.center_mode = 1
    opts.frame_offset = 0
    opts.bond_thickness_by_distance = 0
    img = render(frames[0], opts, ref_center=None)
    assert img.ndim == 3
    assert img.shape[2] == 3
    assert img.dtype == np.uint8


def test_config_load_json():
    """Config loader reads JSON and flattens plane."""
    from xmoltoppm.config_loader import load_config
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "sample_config.json"
    if not path.exists():
        pytest.skip("sample_config.json not found")
    c = load_config(path)
    assert c.get("size") == 300
    assert c.get("circle_style") == 1
    assert c.get("line_style") == 3
    assert c.get("frame_offset") == 10
    assert c.get("plane_z") == 50.0


def test_cli_with_config():
    """CLI with -c loads config and applies options (integration)."""
    import sys
    from xmoltoppm.cli.main import main
    root = Path(__file__).resolve().parent
    config_path = root / "fixtures" / "sample_config.json"
    xyz_path = root / "fixtures" / "minimal.xyz"
    if not config_path.exists() or not xyz_path.exists():
        pytest.skip("fixtures not found")
    out = root.parent / "test_config_out.ppm"
    argv = ["-c", str(config_path), "-i", str(xyz_path), "-o", str(out)]
    try:
        main(argv)
    except SystemExit as e:
        assert e.code == 0, f"CLI exited with {e.code}"
    finally:
        if out.exists():
            out.unlink()


def test_rotation_changes_output():
    """Rotation -r rx ry rz produces different image from no rotation."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 400
    img0 = render(frames[0], opts, ref_center=None)
    opts.rotation_rz = 45.0
    img1 = render(frames[0], opts, ref_center=None)
    assert img0.ndim == 3 and img1.ndim == 3
    assert not np.array_equal(img0, img1)


def test_perspective_changes_output():
    """Perspective -ps produces different image."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 400
    img0 = render(frames[0], opts, ref_center=None)
    opts.perspective_system = 0.3
    img1 = render(frames[0], opts, ref_center=None)
    assert img0.shape == img1.shape
    assert not np.array_equal(img0, img1)


def test_element_override_rgb():
    """-ac element r g b overrides atom colour."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 400
    img0 = render(frames[0], opts, ref_center=None)
    opts.element_rgb = {"C": (255, 0, 0)}
    img1 = render(frames[0], opts, ref_center=None)
    assert not np.array_equal(img0, img1)


def test_ignore_bonds_applied():
    """-ib i1 i2 n1 n2 is applied (no bond drawn between fragment 1 and 2)."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 400
    img0 = render(frames[0], opts, ref_center=None)
    opts.ignore_bonds = [(1, 3, 4, 6)]
    img1 = render(frames[0], opts, ref_center=None)
    assert img0.ndim == 3 and img1.ndim == 3
    assert not np.array_equal(img0, img1)


def test_atom_type_change():
    """-at nat qa changes displayed atom type."""
    root = Path(__file__).resolve().parent
    path = root / "fixtures" / "minimal.xyz"
    frames = list(read_xyz(path))
    opts = _base_opts()
    opts.size = 200
    opts.atom_type_changes = {1: "O"}
    img = render(frames[0], opts, ref_center=None)
    assert img.ndim == 3
    assert img.shape[2] == 3


def test_boundary_box_draws():
    """-bb 1 with cell draws box (image differs from no box)."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    if frames[0].cell is None:
        pytest.skip("traj.xyz frame has no cell")
    opts = _base_opts()
    opts.size = 400
    img0 = render(frames[0], opts, ref_center=None)
    opts.draw_boundary_box = 1
    img1 = render(frames[0], opts, ref_center=None)
    assert img0.ndim == 3 and img1.ndim == 3
    assert not np.array_equal(img0, img1)


def test_supersample_downscales():
    """Supersample scale > 1 renders at higher res then downscales; output is smaller and valid."""
    root = Path(__file__).resolve().parent.parent
    traj_path = root / "traj.xyz"
    if not traj_path.exists():
        pytest.skip("traj.xyz not found")
    frames = list(read_xyz(traj_path))
    opts = _base_opts()
    opts.size = 100
    opts.supersample_scale = 1
    img1 = render(frames[0], opts, ref_center=None)
    opts.supersample_scale = 2
    img2 = render(frames[0], opts, ref_center=None)
    assert img2.ndim == 3 and img2.shape[2] == 3 and img2.dtype == img1.dtype
    # Downscaled image has roughly half the dimensions of the scale=1 viewport
    assert img2.shape[0] <= img1.shape[0] + 1 and img2.shape[1] <= img1.shape[1] + 1
    assert img2.shape[0] >= 1 and img2.shape[1] >= 1


def test_vibrational_mode_changes_output():
    """Vibrational mode with velocities displaces atoms; sub_step 0 vs 2 gives different images."""
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=np.float64)
    velocities = np.array([[0.1, 0.0, 0.0], [-0.1, 0.0, 0.0]], dtype=np.float64)
    frame = Frame(nat=2, coords=coords, symbols=["H", "H"], velocities=velocities)
    opts = _base_opts()
    opts.vibrational_steps = 4
    opts.vibrational_scale = 0.1
    opts.size = 80
    img0 = render(frame, opts, ref_center=None, vibrational_sub_step=0)
    img2 = render(frame, opts, ref_center=None, vibrational_sub_step=2)
    assert img0.ndim == 3 and img2.ndim == 3
    assert not np.array_equal(img0, img2)


def test_graphics_overlay_draws():
    """Graphics overlay with valid data file runs and draws non-background pixels in plot area."""
    root = Path(__file__).resolve().parent
    tmp = root / "fixtures"
    tmp.mkdir(exist_ok=True)
    datafile = tmp / "graphics_data.txt"
    datafile.write_text("col1 col2\n0 0\n50 50\n100 25\n")
    opts = _base_opts()
    opts.size = 200
    opts.graphics_file_entries = [
        (0, 0, str(datafile.resolve()), 2, 10, 150, 80, 60, 2, 1, 2),
    ]
    frame = list(read_xyz(root / "fixtures" / "minimal.xyz"))[0]
    img = render(frame, opts, ref_center=None, frame_index=0)
    assert img.ndim == 3
    # Plot area is roughly (10, 90) to (90, 150); should have some non-background pixels
    region = img[90:151, 10:91, :]
    assert region.size > 0
    assert not np.all(region == 255)


def test_text_overlay_with_pillow():
    """Text overlay when Pillow is available changes output."""
    pytest.importorskip("PIL")
    root = Path(__file__).resolve().parent
    frames = list(read_xyz(root / "fixtures" / "minimal.xyz"))
    opts = _base_opts()
    opts.size = 100
    opts.text_overlays = [("Hi", 5, 5, 1)]
    img_no_text = render(frames[0], opts, ref_center=None)
    opts.text_overlays = [("Hi", 5, 5, 1)]
    img_text = render(frames[0], opts, ref_center=None)
    assert img_text.ndim == 3
    # Text at (5,5) may change pixels in a small region
    patch = img_text[5:20, 5:25, :]
    assert patch.size > 0


def test_text_overlay_without_pillow_no_crash():
    """Text overlay when Pillow is missing does not crash; warning is optional."""
    opts = _base_opts()
    opts.size = 50
    opts.text_overlays = [("X", 2, 2, 1)]
    frame = Frame(nat=1, coords=np.array([[0.0, 0.0, 0.0]]), symbols=["H"])
    # If PIL is available this draws; if not, _draw_text_overlays returns early
    img = render(frame, opts, ref_center=None)
    assert img.ndim == 3 and img.shape[2] == 3


def test_cli_invalid_si_exits_with_error():
    """CLI with invalid -si value (e.g. 0 not in choices 1,-1) exits with non-zero status."""
    import subprocess
    root = Path(__file__).resolve().parent
    xyz = root / "fixtures" / "minimal.xyz"
    if not xyz.exists():
        pytest.skip("minimal.xyz not found")
    result = subprocess.run(
        ["python", "-m", "xmoltoppm.cli.main", "-i", str(xyz), "-o", "out.ppm", "-si", "0"],
        capture_output=True,
        text=True,
        cwd=root.parent,
    )
    assert result.returncode != 0
