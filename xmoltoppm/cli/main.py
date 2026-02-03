"""CLI entrypoint: argparse, single-frame or all-frames."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from xmoltoppm import __version__
from xmoltoppm.config_loader import load_config
from xmoltoppm.core.options import Options
from xmoltoppm.io.ppm_writer import write_ppm, write_png
from xmoltoppm.io.xyz_reader import read_xyz
from xmoltoppm.viz.render import centroid_after_plane_filter, render

# Config file key -> argparse dest (for applying config as parser defaults)
_CONFIG_TO_ARGS = {
    "size": "size", "cc": "cc", "circle_style": "cc",
    "il": "il", "line_style": "il", "ft": "ft", "filename_type": "ft",
    "sh": "sh", "shadow": "sh", "sf": "sf", "sd": "sd",
    "ba": "ba", "bs": "bs", "sc": "sc", "sl": "sl", "bc": "bc", "bl": "bl",
    "ic": "ic", "center_mode": "ic", "nf": "nf", "frame_offset": "nf",
    "bd": "bd", "bond_thickness_by_distance": "bd",
    "supersample_scale": "supersample_scale", "supersample": "supersample_scale",
    "is_strain": "is_strain", "sr": "sr", "iv": "iv", "vv": "vv",
    "ia": "ia", "as_arrow": "as_arrow", "ah": "ah",
}

EPILOG = """
Examples:
  Single frame (first frame of trajectory):
    xmoltoppm -i traj.xyz -o output.ppm
    xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -o output.png --png

  All frames (movie sequence; -af 0 = every frame):
    xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -af 0 -ft 0 -o frames/
"""


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="XYZ trajectory to PPM/PNG molecular visualization (xmoltoppm). Use -h for options and examples.",
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Input XYZ file")
    parser.add_argument("-o", "--output", default="output.ppm", help="Output path (single frame) or base for numbered files")
    parser.add_argument("-s", "--size", type=int, default=500, help="Canvas size in pixels (default 500)")
    parser.add_argument("-cc", "--circle-style", dest="cc", type=int, default=0, choices=(0, 1),
                        help="Circle style: 0 flat, 1 shaded")
    parser.add_argument("-il", "--line-style", dest="il", type=int, default=0, choices=(0, 1, 2, 3, 4),
                        help="Line style: 0 single colour, 1 two-colour, 2 single+shadow, 3 two+shadow, 4 two+shadow+directional")
    parser.add_argument("-af", "--all-frames", dest="af", type=int, default=None, metavar="iskip",
                        help="Output all frames; iskip = frames to skip between outputs (0=every frame)")
    parser.add_argument("-ft", "--filename-type", dest="ft", type=int, default=0, choices=(0, 1, 2),
                        help="Filename type: 0=0001.ppm, 1=1.ppm, 2=molname.ppm")
    parser.add_argument("-sh", "--shadow", dest="sh", type=float, default=0.0,
                        help="Shadow effect 0-1")
    parser.add_argument("-sf", "--shadow-strength", dest="sf", type=float, default=20.0,
                        help="Directional shadow strength (for -il 4)")
    parser.add_argument("-sd", "--shadow-dir", dest="sd", type=str, default="yz", choices=("xy", "xz", "yz"),
                        help="Shadow direction")
    parser.add_argument("-si", dest="si", type=int, default=1, choices=(1, -1),
                        help="Shadow inversion: 1 normal, -1 invert directional shadow (for -il 4)")
    parser.add_argument("-ba", "--bond-adjust", dest="ba", type=float, default=1.0,
                        help="Bond radius adjustment")
    parser.add_argument("-bs", "--border-size", dest="bs", type=float, default=0.5,
                        help="Border size in Angstrom")
    parser.add_argument("-sc", "--circle-size", dest="sc", type=float, default=0.15,
                        help="Circle size factor")
    parser.add_argument("-sl", "--line-size", dest="sl", type=float, default=0.05,
                        help="Line size factor")
    parser.add_argument("-bc", "--background", dest="bc", type=int, default=255,
                        help="Background colour 0-255")
    parser.add_argument("-bl", "--no-boundary-lines", dest="bl", type=int, default=0, choices=(0, 1),
                        help="0: black lines around atoms+bonds, 1: no")
    parser.add_argument("-ic", "--center", dest="ic", type=int, default=0, choices=(0, 1),
                        help="0: keep center from first frame, 1: re-center every frame")
    parser.add_argument("-nf", "--frame-offset", dest="nf", type=int, default=0,
                        help="Add to frame counter in output filenames")
    parser.add_argument("-pp", "--plane", dest="plane", action="append", metavar=("AXIS", "COOR"), nargs=2,
                        help="Remove atoms with axis > coor (e.g. -pp z 10). Repeat for x,y,z.")
    parser.add_argument("-bd", "--bond-thickness", dest="bd", type=int, default=0, choices=(0, 1),
                        help="1: scale bond line thickness by bond distance")
    parser.add_argument("--png", action="store_true", help="Write PNG instead of PPM (requires Pillow)")
    parser.add_argument("--write-run-json", action="store_true",
                        help="Write sidecar JSON with input path, options, and version")
    parser.add_argument("-p", "--preset", dest="preset", type=str, default=None, choices=("movie",),
                        help="Preset: 'movie' = pretty defaults (size 800, shaded circles, line style 4, shadow)")
    parser.add_argument("-c", "--config", dest="config", type=str, default=None,
                        help="Load options from JSON or YAML file (CLI overrides config)")
    parser.add_argument("--supersample", dest="supersample_scale", type=int, default=1, metavar="N",
                        help="Render at Nx resolution then downscale (antialiasing)")
    parser.add_argument("-ar", dest="ar", action="append", nargs=2, metavar=("ELEMENT", "RADIUS"),
                        help="Override circle radius for element (repeat for multiple)")
    parser.add_argument("-ab", dest="ab", action="append", nargs=2, metavar=("ELEMENT", "RADIUS"),
                        help="Override bond radius for element")
    parser.add_argument("-pb", dest="pb", action="append", nargs=3, metavar=("EL1", "EL2", "RADIUS"),
                        help="Override bond cutoff for element pair")
    parser.add_argument("-ac", dest="ac", action="append", nargs=4, metavar=("ELEMENT", "R", "G", "B"),
                        help="Override element colour (0-255)")
    parser.add_argument("-cb", dest="cb", action="append", nargs=5, metavar=("EL1", "EL2", "R", "G", "B"),
                        help="Override bond colour for element pair")
    parser.add_argument("-na", dest="na", nargs=2, metavar=("IATR", "VRADIUS"),
                        help="Only show atoms within VRADIUS of atom IATR (1-based)")
    parser.add_argument("-ir", dest="ir", type=int, default=0, choices=(0, 1),
                        help="1: refresh radius selection every frame")
    parser.add_argument("-ib", dest="ib", action="append", nargs=4, type=int, metavar=("I1", "I2", "N1", "N2"),
                        help="Ignore bonds between fragment [I1..I2] and [N1..N2] (1-based)")
    parser.add_argument("-at", dest="at", action="append", nargs=2, metavar=("NAT", "QA"),
                        help="Change atom index NAT (1-based) to symbol QA")
    parser.add_argument("-r", "--rotation", dest="r", nargs=3, type=float, metavar=("RX", "RY", "RZ"),
                        help="Rotation degrees around x, y, z")
    parser.add_argument("-ps", dest="ps", type=float, default=0.0, help="Perspective ratio (system)")
    parser.add_argument("-pc", dest="pc", type=float, default=0.0, help="Perspective ratio (circles)")
    parser.add_argument("-pl", dest="pl", type=float, default=0.0, help="Perspective ratio (lines)")
    parser.add_argument("-bb", dest="bb", type=int, default=0, choices=(0, 1),
                        help="1: draw boundary box from cell")
    parser.add_argument("-is", dest="is_strain", type=int, default=0, choices=(0, 1),
                        help="1: color atoms by strain (estrain column)")
    parser.add_argument("-sr", dest="sr", type=float, default=1.0, help="Strain color scale (default 1.0)")
    parser.add_argument("-iv", dest="iv", type=int, default=0, choices=(0, 1, 2),
                        help="Velocity coloring: 0 off, 1 temp from estrain as speed, 2 from velocity vector")
    parser.add_argument("-vv", dest="vv", type=float, default=300.0,
                        help="Reference temperature in K for velocity coloring (default 300)")
    parser.add_argument("-ut", dest="ut", nargs=2, type=float, metavar=("TMIN", "TMAX"), default=None,
                        help="Use estrain as temp; gradient from TMIN to TMAX")
    parser.add_argument("-ia", dest="ia", type=int, default=0, choices=(0, 1),
                        help="1: draw velocity arrows")
    parser.add_argument("-as", dest="as_arrow", type=float, default=1.0,
                        help="Arrow scale (default 1.0)")
    parser.add_argument("-ah", dest="ah", type=int, default=10, help="Arrow head size in pixels (default 10)")
    parser.add_argument("-t", dest="t_file", type=str, default=None, metavar="FILE",
                        help="File listing 1-based atom indices for transparent atoms (one per line)")
    parser.add_argument("-ti", dest="ti", action="append", nargs=2, type=int, metavar=("I1", "I2"), default=None,
                        help="Transparent atoms: 1-based index range I1..I2 (repeat for multiple)")
    parser.add_argument("-ri", dest="ri", nargs=3, metavar=("AXIS", "IREND", "IRSTEP"), default=None,
                        help="Rotation animation: axis (x/y/z), total degrees, step per frame")
    parser.add_argument("-tx", dest="tx", action="append", nargs=4, metavar=("TEXT", "IX", "IY", "SCALE"), default=None,
                        help="Text overlay: TEXT at pixel (IX, IY) with SCALE (repeat for multiple)")
    parser.add_argument("-tf", dest="tf", type=str, default=None, metavar="FILE",
                        help="Text overlay file: lines 'start end text ix iy scale [transp]'; text ENERGY/ITERATION/VOLUME/FRAME = frame data")
    parser.add_argument("-gf", dest="gf", type=str, default=None, metavar="FILE",
                        help="Graphics overlay config: lines 'start end datafile ncols ix iy width height point_size col_x col_y'")
    parser.add_argument("-vi", dest="vi", type=int, default=0, metavar="NSTEP",
                        help="Vibrational mode: NSTEP sub-frames per trajectory frame (displace along velocity/mode vector); 0=off")
    parser.add_argument("-vs", dest="vs", type=float, default=1.0, help="Vibrational mode scale (default 1.0)")
    parser.add_argument("-vr", dest="vr", type=int, default=0, metavar="IREPEAT",
                        help="Vibrational mode: repeat cycle IREPEAT times per trajectory frame (default 0)")
    parser.add_argument("--max-frames", "-mf", dest="max_frames", type=int, default=None, metavar="N",
                        help="Limit to first N frames (default: no limit)")
    # Two-pass: get config path, load config, set parser defaults, then parse
    args_pre, _ = parser.parse_known_args(argv)
    config_dict = {}
    if getattr(args_pre, "config", None):
        try:
            config_dict = load_config(args_pre.config)
        except Exception as e:
            print(f"Error: failed to load config: {e}", file=sys.stderr)
            sys.exit(1)
        defaults = {}
        for ckey, arg_attr in _CONFIG_TO_ARGS.items():
            if ckey in config_dict:
                defaults[arg_attr] = config_dict[ckey]
        if defaults:
            parser.set_defaults(**defaults)
    args = parser.parse_args(argv)

    # Apply -pp axis coor -> plane_x/y/z; then config plane_* if no -pp
    plane_default = 5_000_000.0
    args.plane_x = config_dict.get("plane_x", plane_default)
    args.plane_y = config_dict.get("plane_y", plane_default)
    args.plane_z = config_dict.get("plane_z", plane_default)
    for item in (args.plane or []):
        ax, val = item[0].strip().lower(), float(item[1])
        if ax == "x":
            args.plane_x = val
        elif ax == "y":
            args.plane_y = val
        elif ax == "z":
            args.plane_z = val

    if args.preset == "movie":
        args.s = 800
        args.cc = 1
        args.il = 4
        args.sh = 0.6
        args.sf = 20.0
        args.sd = "yz"

    # Build override dicts from -ar, -ab, -pb, -ac, -cb
    args.element_radius_circle = {}
    for item in (getattr(args, "ar", None) or []):
        args.element_radius_circle[item[0].strip()] = float(item[1])
    args.element_radius_bond = {}
    for item in (getattr(args, "ab", None) or []):
        args.element_radius_bond[item[0].strip()] = float(item[1])
    args.pair_bond_cutoff = {}
    for item in (getattr(args, "pb", None) or []):
        k = (item[0].strip(), item[1].strip())
        args.pair_bond_cutoff[k] = float(item[2])
    args.element_rgb = {}
    for item in (getattr(args, "ac", None) or []):
        args.element_rgb[item[0].strip()] = (int(item[1]), int(item[2]), int(item[3]))
    args.pair_bond_rgb = {}
    for item in (getattr(args, "cb", None) or []):
        k = (item[0].strip(), item[1].strip())
        args.pair_bond_rgb[k] = (int(item[2]), int(item[3]), int(item[4]))
    # -na, -ir, -ib, -at
    if getattr(args, "na", None) is not None and len(args.na) == 2:
        args.radius_center_atom = int(args.na[0])
        args.radius_radius = float(args.na[1])
    else:
        args.radius_center_atom = -1
        args.radius_radius = 1e9
    args.radius_refresh = getattr(args, "ir", 0)
    args.ignore_bonds = [(a, b, c, d) for a, b, c, d in (getattr(args, "ib", None) or [])]
    args.atom_type_changes = {}
    for item in (getattr(args, "at", None) or []):
        args.atom_type_changes[int(item[0])] = item[1].strip()
    # -r, -ps, -pc, -pl, -bb
    if getattr(args, "r", None) is not None and len(args.r) == 3:
        args.rotation_rx, args.rotation_ry, args.rotation_rz = float(args.r[0]), float(args.r[1]), float(args.r[2])
    else:
        args.rotation_rx = args.rotation_ry = args.rotation_rz = 0.0
    args.perspective_system = getattr(args, "ps", 0.0)
    args.perspective_circle = getattr(args, "pc", 0.0)
    args.perspective_line = getattr(args, "pl", 0.0)
    args.draw_boundary_box = getattr(args, "bb", 0)
    args.shadow_invert = getattr(args, "si", 1)
    # Strain/velocity coloring and arrows
    args.strain_coloring = getattr(args, "is_strain", 0)
    args.strain_scale = getattr(args, "sr", 1.0)
    args.velocity_coloring = getattr(args, "iv", 0)
    args.velocity_reference_temp = getattr(args, "vv", 300.0)
    ut_arg = getattr(args, "ut", None)
    args.use_temp_range = 1 if (ut_arg is not None and len(ut_arg) == 2) else 0
    args.temp_min = float(ut_arg[0]) if args.use_temp_range else 0.0
    args.temp_max = float(ut_arg[1]) if args.use_temp_range else 1000.0
    args.draw_arrows = getattr(args, "ia", 0)
    args.arrow_scale = getattr(args, "as_arrow", 1.0)
    args.arrow_head_size = getattr(args, "ah", 10)
    # Transparent atoms: -t file and -ti i1 i2
    transparent_set = set()
    t_file = getattr(args, "t_file", None)
    if t_file:
        try:
            with open(t_file) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#"):
                        for part in line.split():
                            try:
                                transparent_set.add(int(part))
                            except ValueError:
                                pass
        except FileNotFoundError:
            print(f"Warning: transparent atoms file not found: {t_file}", file=sys.stderr)
    for item in (getattr(args, "ti", None) or []):
        if len(item) == 2:
            i1, i2 = int(item[0]), int(item[1])
            for idx in range(min(i1, i2), max(i1, i2) + 1):
                transparent_set.add(idx)
    args.transparent_atoms = frozenset(transparent_set)
    # -tx TEXT ix iy scale
    args.text_overlays = []
    for item in (getattr(args, "tx", None) or []):
        if len(item) == 4:
            try:
                args.text_overlays.append((str(item[0]), int(item[1]), int(item[2]), max(1, int(item[3]))))
            except (ValueError, TypeError):
                pass
    # -tf FILE: lines start end text ix iy scale [transp]; special from text (ENERGY=1, ITERATION=2, VOLUME=3)
    args.text_file_entries = []
    tf_path = getattr(args, "tf", None)
    if tf_path:
        try:
            with open(tf_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            start_f, end_f = int(parts[0]), int(parts[1])
                            text_s = parts[2]
                            ix_f, iy_f = int(parts[3]), int(parts[4])
                            scale_f = int(parts[5])
                            transp_f = int(parts[6]) if len(parts) > 6 else 0
                            special = 0
                            if text_s.upper().startswith("ENERGY"):
                                special = 1
                            elif text_s.upper().startswith("ITERATION"):
                                special = 2
                            elif text_s.upper().startswith("VOLUME"):
                                special = 3
                            elif text_s.upper().startswith("FRAME"):
                                special = 4
                            args.text_file_entries.append((start_f, end_f, special, text_s, ix_f, iy_f, scale_f, transp_f))
                        except (ValueError, IndexError):
                            pass
        except FileNotFoundError:
            print(f"Warning: text file not found: {tf_path}", file=sys.stderr)
    # -gf FILE: graphics overlay config (one line per window); resolve datafile relative to config dir
    args.graphics_file_entries = []
    gf_path = getattr(args, "gf", None)
    if gf_path:
        try:
            gf_dir = Path(gf_path).resolve().parent
            with open(gf_path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    if len(parts) >= 11:
                        try:
                            start_g = int(parts[0])
                            end_g = int(parts[1])
                            datafile = parts[2]
                            if not Path(datafile).is_absolute():
                                datafile = str((gf_dir / datafile).resolve())
                            ncols_g = int(parts[3])
                            ix_g = int(parts[4])
                            iy_g = int(parts[5])
                            width_g = int(parts[6])
                            height_g = int(parts[7])
                            point_size_g = int(parts[8])
                            col_x_g = int(parts[9])
                            col_y_g = int(parts[10])
                            args.graphics_file_entries.append(
                                (start_g, end_g, datafile, ncols_g, ix_g, iy_g, width_g, height_g, point_size_g, col_x_g, col_y_g)
                            )
                        except (ValueError, IndexError):
                            pass
        except FileNotFoundError:
            print(f"Warning: graphics config file not found: {gf_path}", file=sys.stderr)
    # -ri axis irend irstep: rotation animation
    ri_arg = getattr(args, "ri", None)
    if ri_arg is not None and len(ri_arg) == 3:
        args.rotation_animation_axis = str(ri_arg[0]).strip().lower()
        args.rotation_animation_step = float(ri_arg[2])  # step per frame (irend unused for now)
    else:
        args.rotation_animation_axis = ""
        args.rotation_animation_step = 0.0

    args.vibrational_steps = getattr(args, "vi", 0)
    args.vibrational_scale = getattr(args, "vs", 1.0)
    args.vibrational_repeat = getattr(args, "vr", 0)

    opts = Options.from_cli(args)
    opts.input_file = args.input
    opts.output_path = args.output

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    try:
        frames = list(read_xyz(input_path))
    except Exception as e:
        print(f"Error: failed to read XYZ file: {e}", file=sys.stderr)
        print("Check XYZ format: line 1 = atom count (integer), line 2 = comment, then N lines of symbol x y z.", file=sys.stderr)
        sys.exit(1)
    if not frames:
        print("Error: no frames read from XYZ file (file may be empty or invalid).", file=sys.stderr)
        print("Check XYZ format: line 1 = atom count (integer), line 2 = comment, then N lines of symbol x y z.", file=sys.stderr)
        sys.exit(1)
    if opts.max_frames is not None and opts.max_frames > 0:
        frames = frames[: opts.max_frames]

    def write_image(path: Path, img) -> None:
        if args.png or (str(path).lower().endswith(".png")):
            try:
                write_png(path, img)
            except ImportError:
                print("PNG output requires Pillow. Install with: pip install Pillow", file=sys.stderr)
                sys.exit(1)
        else:
            write_ppm(path, img)

    def write_run_json(out_path: Path, num_frames: int) -> None:
        sidecar = out_path.with_suffix(out_path.suffix + ".run.json")
        data = {
            "input_file": str(input_path.resolve()),
            "version": __version__,
            "options": {
                "size": opts.size,
                "circle_style": opts.circle_style,
                "line_style": opts.line_style,
                "all_frames": opts.all_frames,
                "frame_skip": opts.frame_skip,
                "filename_type": opts.filename_type,
                "shadow": opts.shadow,
                "shadow_strength": opts.shadow_strength,
                "shadow_direction": opts.shadow_direction,
                "output_path": opts.output_path,
            },
            "num_frames_rendered": num_frames,
        }
        with open(sidecar, "w") as f:
            json.dump(data, f, indent=2)
        print(f"Wrote run metadata to {sidecar}")

    if args.af is not None:
        # All-frames: output every (args.af + 1)-th frame
        iskip = args.af
        out_index = 0
        ext = ".png" if args.png else ".ppm"
        first_out_path = None
        to_render = [(i, f) for i, f in enumerate(frames) if (i % (iskip + 1)) == 0]
        total = len(to_render)
        ref_center = None
        if opts.center_mode == 0 and to_render:
            ref_center = centroid_after_plane_filter(to_render[0][1], opts)
        n_vib = opts.vibrational_steps * (opts.vibrational_repeat + 1) if opts.vibrational_steps > 0 else 0
        for idx, (i, frame) in enumerate(to_render):
            use_vib = n_vib > 0 and frame.velocities is not None
            if use_vib:
                sub_total = n_vib
            else:
                sub_total = 1
            for vib_idx in range(sub_total):
                out_index += 1
                if (total > 1 or sub_total > 1) and vib_idx == 0:
                    print(f"Rendering frame {idx + 1}/{total} ...", file=sys.stderr)
                n = out_index + opts.frame_offset
                if args.ft == 0:
                    out_name = f"{n:04d}{ext}"
                elif args.ft == 1:
                    out_name = f"{n}{ext}"
                else:
                    molname = (frame.molname or "frame").strip()
                    out_name = f"{molname}{ext}"
                out_path = Path(args.output).parent / out_name
                if Path(args.output).parent == Path("."):
                    out_path = Path(out_name)
                if first_out_path is None:
                    first_out_path = out_path
                sub_step = (vib_idx % opts.vibrational_steps) if use_vib else None
                img = render(
                    frame, opts, ref_center=ref_center, frame_index=out_index - 1,
                    vibrational_sub_step=sub_step,
                )
                write_image(out_path, img)
                print(f"Wrote {out_path}")
        if getattr(args, "write_run_json", False) and first_out_path is not None:
            write_run_json(first_out_path, out_index)
        print(f"Normal end; wrote {out_index} frame(s)")
    else:
        # Single frame: first frame (ref_center only used when all_frames and -ic 0)
        frame = frames[0]
        n_vib = opts.vibrational_steps * (opts.vibrational_repeat + 1) if opts.vibrational_steps > 0 else 0
        use_vib = n_vib > 0 and frame.velocities is not None
        if use_vib:
            base_path = Path(args.output)
            if args.png:
                base_path = base_path.with_suffix(".png")
            for vib_idx in range(n_vib):
                sub_step = vib_idx % opts.vibrational_steps
                n = vib_idx + 1 + opts.frame_offset
                out_name = f"{n:04d}{base_path.suffix}" if opts.filename_type == 0 else f"{n}{base_path.suffix}"
                out_path = base_path.parent / out_name
                if base_path.parent == Path("."):
                    out_path = Path(out_name)
                img = render(frame, opts, ref_center=None, frame_index=vib_idx, vibrational_sub_step=sub_step)
                write_image(out_path, img)
                print(f"Wrote {out_path}")
            if getattr(args, "write_run_json", False):
                write_run_json(out_path, n_vib)
            print(f"Normal end of program; wrote {n_vib} frame(s)")
        else:
            img = render(frame, opts, ref_center=None, frame_index=0)
            out_path = Path(args.output)
            if args.png:
                out_path = out_path.with_suffix(".png")
            write_image(out_path, img)
            if getattr(args, "write_run_json", False):
                write_run_json(out_path, 1)
            print(f"Normal end of program; output in {out_path}")


if __name__ == "__main__":
    main()
