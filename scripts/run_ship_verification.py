#!/usr/bin/env python3
"""
Run ship verification: for each feature group, render full video from traj.xyz
(-af 0), measure wall time, and write outputs + timing report to verification_review/.
Run from project root: python scripts/run_ship_verification.py
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TRAJ = ROOT / "traj.xyz"
OUT_BASE = ROOT / "verification_review"
REPORT = OUT_BASE / "timing_report.txt"

# (group_name, extra_argv, description)
# -af 0 and -o OUT_BASE/group/frame.ppm are added by the runner
GROUPS = [
    ("baseline", [], "Minimal: default size, circle and line style"),
    ("circle_line_styles_il0", ["-cc", "0", "-il", "0"], "Circle flat, line single colour"),
    ("circle_line_styles_il3", ["-cc", "1", "-il", "3"], "Circle shaded, line two+shadow"),
    ("circle_line_styles_il4", ["-cc", "1", "-il", "4"], "Line two+shadow+directional"),
    ("shadow", ["-cc", "1", "-il", "4", "-sh", "0.6", "-sf", "20", "-sd", "yz", "-si", "1"], "Shadow normal"),
    ("shadow_invert", ["-cc", "1", "-il", "4", "-sh", "0.6", "-si", "-1"], "Shadow inverted"),
    ("png_output", ["--png"], "PNG instead of PPM (needs Pillow)"),
    ("overlays_text", ["-tf", str(ROOT / "tests" / "fixtures" / "overlay_frame.txt")], "Text overlay: Frame 1, Frame 2, ... (needs Pillow)"),
    ("overlays_graphics", ["-gf", str(ROOT / "tests" / "fixtures" / "gf_config.txt")], "Graphics overlay"),
    ("energy_overlay", ["-gfe"], "Energy vs iteration overlay from XYZ comment line"),
    ("strain_coloring", ["-is", "1", "-sr", "1.0"], "Strain coloring (estrain column)"),
    ("transparent", ["-ti", "1", "3"], "Transparent atoms range"),
    ("rotation_anim", ["-ri", "z", "360", "2"], "Rotation animation"),
    ("max_frames", ["-mf", "2"], "Limit to first 2 frames (only group using -mf)"),
    ("config_preset", ["-c", str(ROOT / "tests" / "fixtures" / "sample_config.json")], "Config file"),
    ("preset_movie", ["-p", "movie"], "Preset movie"),
    ("recenter", ["-ic", "1"], "Re-center every frame"),
    ("supersample", ["--supersample", "2"], "2x supersample then downscale"),
    ("boundary_box", ["-bb", "1"], "Draw boundary box from cell"),
    ("write_run_json", ["--write-run-json"], "Sidecar run JSON"),
]


def main() -> None:
    if not TRAJ.exists():
        print(f"Error: {TRAJ} not found. Run from project root.", file=sys.stderr)
        sys.exit(1)

    OUT_BASE.mkdir(parents=True, exist_ok=True)
    with open(REPORT, "w") as f:
        f.write("group\tframes\twall_seconds\n")
    total_seconds = 0.0
    skipped: list[str] = []

    for group_name, extra_argv, _desc in GROUPS:
        out_dir = OUT_BASE / group_name
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / "frame.ppm"

        # max_frames: -af 0 -mf 2 so we get exactly 2 frames for speed
        if group_name == "max_frames":
            base_argv = [
                sys.executable, "-m", "xmoltoppm.cli.main",
                "-i", str(TRAJ),
                "-af", "0",
                "-o", str(out_path),
            ]
            cmd = base_argv + extra_argv  # extra_argv = ["-mf", "2"]
        else:
            base_argv = [
                sys.executable, "-m", "xmoltoppm.cli.main",
                "-i", str(TRAJ),
                "-af", "0",
                "-o", str(out_path),
            ]
            cmd = base_argv + extra_argv

        print(f"Running: {group_name} ...", file=sys.stderr)
        t0 = time.perf_counter()
        try:
            result = subprocess.run(
                cmd,
                cwd=str(ROOT),
                capture_output=True,
                text=True,
                timeout=600,
            )
            elapsed = time.perf_counter() - t0
        except subprocess.TimeoutExpired:
            with open(REPORT, "a") as f:
                f.write(f"{group_name}\tTIMEOUT\t-\n")
            skipped.append(group_name)
            print(f"  Timeout.", file=sys.stderr)
            continue
        except Exception as e:
            with open(REPORT, "a") as f:
                f.write(f"{group_name}\tERROR\t{e!s}\n")
            skipped.append(group_name)
            print(f"  Error: {e}", file=sys.stderr)
            continue

        if result.returncode != 0:
            with open(REPORT, "a") as f:
                f.write(f"{group_name}\tFAIL\t{result.returncode}\n")
            skipped.append(group_name)
            print(f"  stderr: {result.stderr[:200]}", file=sys.stderr)
            continue

        # Count frames: number of .ppm or .png in out_dir
        count = sum(1 for _ in out_dir.glob("*.ppm")) + sum(1 for _ in out_dir.glob("*.png"))
        with open(REPORT, "a") as f:
            f.write(f"{group_name}\t{count}\t{elapsed:.2f}\n")
        total_seconds += elapsed
        print(f"  {count} frames, {elapsed:.2f}s", file=sys.stderr)

    # Append summary to report
    with open(REPORT, "a") as f:
        f.write("\n")
        f.write(f"# Total wall time (completed groups): {total_seconds:.2f}s\n")
        if skipped:
            f.write(f"# Skipped/failed: {', '.join(skipped)}\n")
    print(f"\nTiming report: {REPORT}", file=sys.stderr)


if __name__ == "__main__":
    main()
