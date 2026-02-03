# xmoltoppm

**Turn XYZ molecular dynamics trajectories into images or movies in one command.**

For the **full capability set** (all options and flags), run `xmoltoppm -h` or open **`docs/guide.md`** (reference-style; detailed options and workflows).

[![Tests](https://github.com/asepehri93/xmoltoppm/actions/workflows/tests.yml/badge.svg)](https://github.com/asepehri93/xmoltoppm/actions/workflows/tests.yml)
[![PyPI version](https://img.shields.io/pypi/v/xmoltoppm)](https://pypi.org/project/xmoltoppm/) — *[Publish to PyPI](CONTRIBUTING.md#publishing-to-pypi) once so the badge and `pip install xmoltoppm` work.*
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

This package is a Python port of the original Fortran **xmoltoppm** by **Adri C.T. van Duin**. See [What's new in this version](#whats-new-in-this-version) below for updates in this port.

---

## What it does

xmoltoppm reads XYZ-format trajectory files and draws 2D pictures: atoms as coloured circles, bonds by distance. You get a single snapshot or a numbered sequence of frames ready for video encoding. No interactive prompts—everything is controlled from the command line.

---

## ✨ Features

- **Single command to movie**: Render every frame (or every Nth) from a trajectory into numbered PPM or PNG files; feed them to ffmpeg and you have a movie.
- **One snapshot in seconds**: Export the first frame (or any frame) with size, style, and format options in one line.
- **Flexible XYZ input**: Works with standard 4-column XYZ and with extended formats (e.g. ReaxFF comment line, optional estrain/velocity columns).
- **PNG and run metadata**: Use `--png` for PNG output (with Pillow); use `--write-run-json` to record input path, options, and version for reproducibility.
- **Rich styling**: Shaded circles, multiple line styles (including directional shadow), strain or velocity coloring, velocity arrows, text and graphics overlays, rotation animation, transparent atoms, config files and presets.

---

## 📦 Installation

Get up and running with pip:

```bash
pip install xmoltoppm
```

Optional PNG support (recommended if you want PNG output):

```bash
pip install xmoltoppm[png]
# or
pip install Pillow
```

Requires Python 3.9+ and NumPy.

---

## Quick Start

The easiest way to use xmoltoppm:

```bash
xmoltoppm -i traj.xyz -o output.ppm
xmoltoppm -i traj.xyz -af 0 -ft 0 -o frames/
xmoltoppm -i traj.xyz --preset movie -af 0 -o movie_frames/
```

- Line 1: first frame only → one image.
- Line 2: all frames → `0001.ppm`, `0002.ppm`, … in `frames/`.
- Line 3: same with pretty defaults (size 800, shaded circles, shadow).

---

## 📖 Practical Examples

### Single snapshot

Get one image from the first frame:

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -o snapshot.ppm
```

For PNG:

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -o snapshot.png --png
```

### Full trajectory to movie frames

Render every frame into numbered files, then encode with ffmpeg:

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -af 0 -ft 0 -o frames/
# Produces frames/0001.ppm, frames/0002.ppm, ...
ffmpeg -framerate 25 -i frames/%04d.ppm -pix_fmt yuv420p movie.mp4
```

### Pretty movie preset and run metadata

Use the built-in “movie” preset and record options for reproducibility:

```bash
xmoltoppm -i traj.xyz --preset movie -af 0 -ft 0 -o movie_frames/ --write-run-json
```

This uses size 800, shaded circles, line style 4 with shadow, and writes a `.run.json` next to the first output with input path, options, and version.

---

## Key options

| Option | Description |
|--------|-------------|
| `-i`, `--input` | Input XYZ file (required) |
| `-o`, `--output` | Output path or directory base (default: `output.ppm`) |
| `-s`, `--size` | Canvas size in pixels (default: 500) |
| `-cc` | Circle style: 0 flat, 1 shaded |
| `-il` | Line style: 0–4 (single/two-colour, shadow, directional) |
| `-af iskip` | Output all frames; iskip = frames to skip (0 = every frame) |
| `-ft` | Filename type: 0 `0001.ppm`, 1 `1.ppm`, 2 `molname.ppm` |
| `--png` | Write PNG instead of PPM (requires Pillow) |
| `--preset movie` | Pretty defaults: size 800, shaded, line style 4, shadow |
| `-c`, `--config` | Load options from JSON or YAML file |
| `--max-frames N` | Limit to first N frames |
| `--write-run-json` | Write sidecar JSON with input, options, version |

Run `xmoltoppm -h` for the full list of options.

---

## XYZ format

- **Line 1**: Number of atoms (integer).
- **Line 2**: Comment. If it has 9 space-separated fields, they are parsed as molname, iteration, energy, cell lengths/angles; otherwise defaults are used.
- **Next N lines**: One per atom. Minimum: `symbol x y z`. Optional extra columns: type, estrain, or velocity components.

So any standard XYZ with at least 4 columns per atom works.

---

## FAQ / Troubleshooting

- **“PNG output requires Pillow”**  
  Install Pillow: `pip install Pillow` or `pip install xmoltoppm[png]`.

- **Does it work with my XYZ file?**  
  Yes, if each atom line has at least `symbol x y z`. Extra columns (e.g. type, estrain, velocities) are used when present for coloring, arrows, or vibrational mode.

- **No frames written / “no frames read”**  
  Check that the first line is the atom count (integer) and the next N lines are atom lines. See `docs/guide.md` for detailed XYZ format.

- **Why does shadow seem to move or flip between frames?**  
  Shadow is **view-depth (z) based**: “lit” atoms follow which side faces the viewer. As the molecule rotates, the lit side can shift (e.g. top vs bottom). Line style 4 adds a z-based component on bonds; for a fixed world light direction, a future version may add an option.

- **Do output images have a fixed size for movies?**  
  Yes. Every frame is rendered at the same pixel dimensions (square canvas from `-s`), so ffmpeg and other encoders get a consistent sequence.

---

## What's new in this version

This port adds the following relative to the original Fortran xmoltoppm:

- **Full CLI** — No interactive prompts or temp files; all options via arguments.
- **PNG and run metadata** — `--png` (with Pillow) and `--write-run-json` for reproducibility.
- **Robust XYZ** — Standard 4-column plus optional ReaxFF header and extra columns (estrain, velocities).
- **Fixed canvas for movies** — Output dimensions stay constant across frames so ffmpeg and other encoders get a consistent sequence.
- **Modular code** — Separate IO, core, viz, and CLI; tests and extensions are straightforward.

---

## 🤝 Contributing and support

Found a bug or have a question? Open an issue on the project repository. See [CONTRIBUTING.md](CONTRIBUTING.md) for development setup, tests, and how to add options. For detailed documentation (XYZ format, all options, workflows), see **`docs/guide.md`**. See [CHANGELOG.md](CHANGELOG.md) for version history.

---

## Appendix A: verification — where to look and what to check

After running `python scripts/run_ship_verification.py` (requires `traj.xyz` in the project root), outputs go to `verification_review/<group>/`. Use this to evaluate accuracy and quality:

| Group | Where | What to check |
|-------|--------|----------------|
| **baseline** | `verification_review/baseline/` | 102 PPMs; default style; molecule centered, bonds visible; **fixed image size** every frame. |
| **circle_line_styles_il0 / il3 / il4** | `…/circle_line_styles_il0/`, `il3/`, `il4/` | il0: flat circles, single-colour bonds. il3: shaded circles, two-colour + shadow. il4: adds z-based directional component on bonds (may look similar to il3). |
| **shadow** / **shadow_invert** | `…/shadow/`, `…/shadow_invert/` | Shadow is **view-depth (z) based**; as the molecule rotates, lit atoms can shift (top/bottom). shadow_invert flips the bond-direction component. |
| **png_output** | `…/png_output/` | 102 PNGs; same content as baseline, PNG format. |
| **overlays_text** | `…/overlays_text/` | Text shows **Frame 1**, **Frame 2**, … (output frame number) at fixed position. |
| **overlays_graphics** | `…/overlays_graphics/` | Small inset plot in a **fixed corner**, not over the atoms. |
| **strain_coloring** | `…/strain_coloring/` | Atoms coloured by estrain; colour gradient vs baseline. |
| **transparent** | `…/transparent/` | Atoms 1–3 outline-only (transparent fill). |
| **rotation_anim** | `…/rotation_anim/` | Frames 0001→0102: molecule rotates 360° around z. |
| **max_frames** | `…/max_frames/` | Exactly 2 frames (0001, 0002). |
| **config_preset** | `…/config_preset/` | Options from `sample_config.json`. |
| **preset_movie** | `…/preset_movie/` | Size 800, shaded, line style 4, shadow. |
| **recenter** | `…/recenter/` | Molecule re-centred each frame. |
| **supersample** | `…/supersample/` | Sharper (2× render then downscale). |
| **boundary_box** | `…/boundary_box/` | **2D frame** (0.95× inset rectangle only); view scaled to molecule; no 3D box atoms. |
| **write_run_json** | `…/write_run_json/` | 102 PPMs plus `0001.ppm.run.json` with input, options, version. |

---

## Appendix B: Example run times

Typical wall times for rendering **all frames** of a 102-frame trajectory (6 atoms per frame) on a single run:

| Scenario | Frames | Wall time (s) |
|---------|--------|----------------|
| Baseline (default style) | 102 | ~63 |
| PNG output | 102 | ~14 |
| Text overlay | 102 | ~97 |
| Strain coloring | 102 | ~66 |

Times are from `verification_review/timing_report.txt` (run `python scripts/run_ship_verification.py` to regenerate). Your times will depend on trajectory size and hardware.
