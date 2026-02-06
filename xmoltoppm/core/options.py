"""Options: render and CLI options (mirrors Fortran -s, -cc, -il, etc.)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Default plane values: no filtering (Fortran 5000000.0)
_PLANE_DEFAULT = 5_000_000.0


@dataclass
class Options:
    """Render and output options."""

    size: int = 500
    circle_style: int = 0  # 0 flat, 1 shaded
    line_style: int = 0  # 0-4
    all_frames: bool = False
    frame_skip: int = 0  # iskip: output every (frame_skip+1)-th frame
    filename_type: int = 0  # 0: 0001.ppm, 1: 1.ppm, 2: molname.ppm
    shadow: float = 0.0  # 0 no, 1 max
    shadow_strength: float = 20.0  # directional shadow (il=4)
    shadow_direction: str = "yz"  # xy, xz, yz
    shadow_invert: int = 1  # -si: 1 normal, -1 invert directional shadow (il=4)
    input_file: Optional[str] = None
    output_path: str = "output.ppm"
    # Extra Fortran-equivalent defaults
    bond_radius_adjust: float = 1.0  # ba
    border_size: float = 0.5  # brdersize
    circle_size: float = 0.15  # vcirc
    line_size: float = 0.05  # vline
    background: int = 255  # ibagr
    draw_boundary_lines: int = 0  # ibline: 0 black lines, 1 no
    # -ic: 0 keep center from first frame, 1 re-center every frame
    center_mode: int = 0
    # -nf: add to frame counter in output filenames/text
    frame_offset: int = 0
    # -pp: remove atoms with coord > plane (per-axis; default no filter)
    plane_x: float = _PLANE_DEFAULT
    plane_y: float = _PLANE_DEFAULT
    plane_z: float = _PLANE_DEFAULT
    # -bd: 0 fixed line thickness, 1 scale by bond distance
    bond_thickness_by_distance: int = 0
    # Supersample: render at size*scale then downscale (antialiasing)
    supersample_scale: int = 1
    # Per-element overrides: -ar (circle radius), -ab (bond radius), -ac (rgb), -pb (pair bond), -cb (pair bond color)
    element_radius_circle: Dict[str, float] = field(default_factory=dict)
    element_radius_bond: Dict[str, float] = field(default_factory=dict)
    element_rgb: Dict[str, Tuple[int, int, int]] = field(default_factory=dict)
    pair_bond_cutoff: Dict[Tuple[str, str], float] = field(default_factory=dict)
    pair_bond_rgb: Dict[Tuple[str, str], Tuple[int, int, int]] = field(default_factory=dict)
    # -na iatr vradius: only show atoms within vradius of atom iatr (1-based); -ir refresh every frame
    radius_center_atom: int = -1  # 1-based, <=0 = off
    radius_radius: float = 1e9
    radius_refresh: int = 0  # 1 = recompute every frame
    # -ib i1 i2 n1 n2: ignore bonds between fragment [i1..i2] and [n1..n2] (1-based); list of (i1,i2,n1,n2)
    ignore_bonds: List[Tuple[int, int, int, int]] = field(default_factory=list)
    # -at nat qa: change atom index nat (1-based) to symbol qa
    atom_type_changes: Dict[int, str] = field(default_factory=dict)
    # -r rx ry rz: rotation degrees around x, y, z
    rotation_rx: float = 0.0
    rotation_ry: float = 0.0
    rotation_rz: float = 0.0
    # -ps, -pc, -pl: perspective ratio (0 = none)
    perspective_system: float = 0.0
    perspective_circle: float = 0.0
    perspective_line: float = 0.0
    # -bb: draw boundary box (1 = yes)
    draw_boundary_box: int = 0
    # -is, -sr: strain coloring (1 = use estrain), scale
    strain_coloring: int = 0
    strain_scale: float = 1.0
    # -iv, -vv: velocity coloring (1=temp from estrain, 2=from velocity vector), reference temp K
    velocity_coloring: int = 0
    velocity_reference_temp: float = 300.0
    # -ut: use temp range from input (estrain as temp), min/max for gradient
    use_temp_range: int = 0
    temp_min: float = 0.0
    temp_max: float = 1000.0
    # -ia, -as, -ah: draw velocity arrows, scale, head size (pixels)
    draw_arrows: int = 0
    arrow_scale: float = 1.0
    arrow_head_size: int = 10
    # -t file, -ti i1 i2: transparent atoms (1-based indices)
    transparent_atoms: frozenset = field(default_factory=frozenset)
    # -tx TEXT ix iy scale: text overlay (repeat for multiple)
    text_overlays: List[Tuple[str, int, int, int]] = field(default_factory=list)
    # -tf FILE: text file entries (start_frame, end_frame, special, text, ix, iy, scale, transp); special 0=literal 1=ENERGY 2=ITERATION 3=VOLUME
    text_file_entries: List[Tuple[int, int, int, str, int, int, int, int]] = field(default_factory=list)
    # -ri axis irend irstep: rotation animation (axis x/y/z, step degrees per frame)
    rotation_animation_axis: str = ""
    rotation_animation_step: float = 0.0
    # -gf FILE: graphics overlay entries (start, end, datafile, ncols, ix, iy, width, height, point_size, col_x, col_y)
    graphics_file_entries: List[Tuple[int, int, str, int, int, int, int, int, int, int, int]] = field(default_factory=list)
    # -vi nstep, -vs vscale, -vr irepeat: vibrational mode (displace coords along velocity/mode vector)
    vibrational_steps: int = 0
    vibrational_scale: float = 1.0
    vibrational_repeat: int = 0
    # --max-frames N: limit to first N frames (None = no limit)
    max_frames: Optional[int] = None
    # Pixel bbox (px_min, py_min, px_max, py_max) to avoid when placing graphics overlay (e.g. from -gfe)
    graphics_avoid_bbox: Optional[Tuple[int, int, int, int]] = None
    # When True, use side-by-side layout: molecule left, graphics overlay right (no overlay on top of molecule)
    graphics_side_by_side: bool = False

    @classmethod
    def from_cli(cls, args) -> Options:
        """Build Options from argparse namespace (see cli.main)."""
        return cls(
            size=getattr(args, "size", 500),
            circle_style=getattr(args, "cc", 0),
            line_style=getattr(args, "il", 0),
            all_frames=getattr(args, "af", None) is not None,
            frame_skip=int(args.af) if getattr(args, "af", None) is not None else 0,
            filename_type=getattr(args, "ft", 0),
            shadow=getattr(args, "sh", 0.0),
            shadow_strength=getattr(args, "sf", 20.0),
            shadow_direction=getattr(args, "sd", "yz"),
            shadow_invert=getattr(args, "shadow_invert", 1),
            input_file=getattr(args, "input", None),
            output_path=getattr(args, "output", "output.ppm"),
            bond_radius_adjust=getattr(args, "ba", 1.0),
            border_size=getattr(args, "bs", 0.5),
            circle_size=getattr(args, "sc", 0.15),
            line_size=getattr(args, "sl", 0.05),
            background=getattr(args, "bc", 255),
            draw_boundary_lines=getattr(args, "bl", 0),
            center_mode=getattr(args, "ic", 0),
            frame_offset=getattr(args, "nf", 0),
            plane_x=getattr(args, "plane_x", _PLANE_DEFAULT),
            plane_y=getattr(args, "plane_y", _PLANE_DEFAULT),
            plane_z=getattr(args, "plane_z", _PLANE_DEFAULT),
            bond_thickness_by_distance=getattr(args, "bd", 0),
            supersample_scale=getattr(args, "supersample_scale", 1),
            element_radius_circle=getattr(args, "element_radius_circle", {}),
            element_radius_bond=getattr(args, "element_radius_bond", {}),
            element_rgb=getattr(args, "element_rgb", {}),
            pair_bond_cutoff=getattr(args, "pair_bond_cutoff", {}),
            pair_bond_rgb=getattr(args, "pair_bond_rgb", {}),
            radius_center_atom=getattr(args, "radius_center_atom", -1),
            radius_radius=getattr(args, "radius_radius", 1e9),
            radius_refresh=getattr(args, "radius_refresh", 0),
            ignore_bonds=getattr(args, "ignore_bonds", []),
            atom_type_changes=getattr(args, "atom_type_changes", {}),
            rotation_rx=getattr(args, "rotation_rx", 0.0),
            rotation_ry=getattr(args, "rotation_ry", 0.0),
            rotation_rz=getattr(args, "rotation_rz", 0.0),
            perspective_system=getattr(args, "perspective_system", 0.0),
            perspective_circle=getattr(args, "perspective_circle", 0.0),
            perspective_line=getattr(args, "perspective_line", 0.0),
            draw_boundary_box=getattr(args, "draw_boundary_box", 0),
            strain_coloring=getattr(args, "strain_coloring", 0),
            strain_scale=getattr(args, "strain_scale", 1.0),
            velocity_coloring=getattr(args, "velocity_coloring", 0),
            velocity_reference_temp=getattr(args, "velocity_reference_temp", 300.0),
            use_temp_range=getattr(args, "use_temp_range", 0),
            temp_min=getattr(args, "temp_min", 0.0),
            temp_max=getattr(args, "temp_max", 1000.0),
            draw_arrows=getattr(args, "draw_arrows", 0),
            arrow_scale=getattr(args, "arrow_scale", 1.0),
            arrow_head_size=getattr(args, "arrow_head_size", 10),
            transparent_atoms=getattr(args, "transparent_atoms", frozenset()),
            text_overlays=getattr(args, "text_overlays", []),
            text_file_entries=getattr(args, "text_file_entries", []),
            rotation_animation_axis=getattr(args, "rotation_animation_axis", ""),
            rotation_animation_step=getattr(args, "rotation_animation_step", 0.0),
            graphics_file_entries=getattr(args, "graphics_file_entries", []),
            vibrational_steps=getattr(args, "vibrational_steps", 0),
            vibrational_scale=getattr(args, "vibrational_scale", 1.0),
            vibrational_repeat=getattr(args, "vibrational_repeat", 0),
            max_frames=getattr(args, "max_frames", None),
        )
