# xmoltoppm User Guide

Detailed documentation for the Python xmoltoppm application: XYZ format, options, workflow, and programmatic use.

---

## 1. Overview

xmoltoppm reads XYZ-format molecular dynamics trajectories and produces 2D raster images (PPM or PNG). Each frame is projected onto the xy-plane (z used for depth ordering and optional shading). Atoms are drawn as coloured circles; bonds are inferred from interatomic distances and drawn as lines. Output can be a single image or a numbered sequence for video encoding.

---

## 2. XYZ Input Format

### 2.1 Structure

- **Line 1**: Integer — number of atoms in the frame.
- **Line 2**: Comment line. Parsed if it contains 9 space-separated fields (ReaxFF-style).
- **Next N lines**: One line per atom.

Frames repeat: after N atom lines, the next line is the atom count of the next frame.

### 2.2 Comment Line (Optional 9-Field)

If the second line has 9 fields, they are interpreted as:

| Field   | Meaning        | Example   |
|--------|-----------------|-----------|
| 1      | Molecule name   | Acetylene_O2 |
| 2      | Iteration/step  | 0         |
| 3      | Energy          | -605.26   |
| 4–6    | Cell lengths a, b, c (Å) | 80 80 80 |
| 7–9    | Cell angles α, β, γ (degrees) | 90 90 90 |

If parsing fails or there are fewer than 9 fields, defaults are used (e.g. empty name, zero iteration/energy, no cell) and **coordinates are still read**. This allows standard XYZ files from other tools to work.

### 2.3 Atom Lines

Minimum: **4 columns** — `symbol x y z`.

Optional extra columns (used when present):

- Columns 5–6: e.g. type index, estrain (or other scalar).
- Columns 7–9: velocity components vx, vy, vz.

Example:

```
H   38.45539  41.24817  38.89843     1          2.80642
```

Symbols are matched to the built-in element table (CPK-style colours and radii); unknown symbols default to carbon.

---

## 3. Command-Line Reference

### 3.1 Input and output

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-i`, `--input` | path | (required) | Input XYZ file. |
| `-o`, `--output` | path | output.ppm | Output path (single frame) or base for numbered files; with `-af`, filenames from `-ft` and frame index. |
| `--png` | flag | off | Write PNG instead of PPM (requires Pillow). |
| `--write-run-json` | flag | off | Write sidecar `.run.json` with input path, options, and version. |

### 3.2 Canvas and layout

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-s`, `--size` | int | 500 | Canvas size in pixels; output is square (size × size). |
| `-bs` | float | 0.5 | Border size in Å around the molecule. |
| `-ic` | 0 or 1 | 0 | Center: 0 keep from first frame, 1 re-center every frame. |
| `-nf` | int | 0 | Add to frame counter in output filenames. |

### 3.4 Atom and Bond Appearance

- **`-cc`** — Circle style: `0` flat, `1` shaded (z-dependent).
- **`-sc`** — Circle size factor. Default: 0.15.
- **`-il` / `--line-style`** — Line (bond) style: `0` single colour, `1` two-colour (gradient), `2` single + shadow, `3` two-colour + shadow, `4` two-colour + directional shadow. Use `3` or `4` for a more 3D look.
- **`-sl`** — Line thickness factor. Default: 0.05.
- **`-ba`** — Bond radius adjustment (distance cutoff scale). Default: 1.0.
- **`-bl`** — `0` draw black outline around atoms/bonds, `1` no black outline.
- **`-bc`** — Background colour (0–255). Default: 255 (white).

### 3.4 Shadow (line styles 2–4)

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-sh` | float | 0.0 | Shadow strength (0–1). |
| `-sf` | float | 20.0 | Directional shadow strength (for `-il 4`). |
| `-sd` | xy \| xz \| yz | yz | Shadow direction. |
| `-si` | 1 or -1 | 1 | Shadow inversion for directional shadow (`-il 4`): 1 normal, -1 invert. |

### 3.5 All-frames and preset

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-af iskip` | int | — | Output all frames; iskip = frames to skip (0 = every frame). |
| `-ft` | 0–2 | 0 | Filename type: 0 0001.ppm, 1 1.ppm, 2 molname.ppm. |
| `--max-frames`, `-mf` | int | — | Limit to first N frames. |
| `-p`, `--preset` | movie | — | Preset "movie": size 800, cc 1, il 4, shadow 0.6, sf 20, sd yz. |

### 3.6 Filtering and selection

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-pp` axis coor | append | — | Remove atoms with coord > coor on axis (e.g. `-pp z 10`); repeat for x, y, z. Defaults (no filter) are large. |
| `-na` IATR VRADIUS | 2 args | — | Only show atoms within VRADIUS of atom IATR (1-based). |
| `-ir` | 0 or 1 | 0 | 1: refresh radius selection every frame (with `-na`). |
| `-ib` I1 I2 N1 N2 | 4 ints | — | Ignore bonds between fragment [I1..I2] and [N1..N2] (1-based); repeat for multiple. |
| `-at` NAT QA | 2 args | — | Change atom index NAT (1-based) to symbol QA; repeat for multiple. |

### 3.7 Rotation, perspective, and boundary box

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-r` RX RY RZ | 3 floats | — | Rotation in degrees around x, then y, then z. |
| `-ri` AXIS IREND IRSTEP | 3 args | — | Rotation animation: axis (x/y/z), total degrees, step per frame. |
| `-ps`, `-pc`, `-pl` | float | 0.0 | Perspective ratio (system, circles, lines); 0 = none. |
| `-bb` | 0 or 1 | 0 | 1: draw boundary box. Draws a **2D frame only** — a black rectangle at 0.95× image size (inset from edges); no 3D box atoms or bonds. Requires `frame.cell` in XYZ. Viewport is scaled to molecule. |

### 3.8 Coloring and velocity

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-is` | 0 or 1 | 0 | 1: colour atoms by strain (estrain column). |
| `-sr` | float | 1.0 | Strain colour scale. |
| `-iv` | 0–2 | 0 | Velocity colouring: 0 off, 1 temp from estrain as speed, 2 from velocity vector. |
| `-vv` | float | 300.0 | Reference temperature (K) for velocity colouring. |
| `-ut` TMIN TMAX | 2 floats | — | Use estrain as temp; gradient from TMIN to TMAX. |
| `-ia` | 0 or 1 | 0 | 1: draw velocity arrows. |
| `-as` | float | 1.0 | Arrow scale. |
| `-ah` | int | 10 | Arrow head size in pixels. |

### 3.9 Overlays

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-tx` TEXT IX IY SCALE | 4 args | — | Text overlay at pixel (IX, IY) with SCALE; repeat for multiple. Requires Pillow. |
| `-tf` FILE | path | — | Text overlay file; see §3.15 Text overlay file format. |
| `-gf` FILE | path | — | Graphics overlay config; see §3.16 Graphics overlay config. |

### 3.10 Transparent atoms and vibrational mode

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-t` FILE | path | — | File listing 1-based atom indices for transparent (outline-only) atoms, one per line. |
| `-ti` I1 I2 | 2 ints | — | Transparent atoms in range I1..I2 (1-based); repeat for multiple. |
| `-vi` NSTEP | int | 0 | Vibrational mode: NSTEP sub-frames per trajectory frame (displace along velocity/mode); 0 = off. |
| `-vs` | float | 1.0 | Vibrational mode scale. |
| `-vr` IREPEAT | int | 0 | Vibrational mode: repeat cycle IREPEAT times per frame. |

### 3.11 Config and element overrides

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `-c`, `--config` FILE | path | — | Load options from JSON or YAML; CLI overrides config. See §3.17 Config file. |
| `--supersample` N | int | 1 | Render at N× resolution then downscale (antialiasing). |
| `-ar` ELEMENT RADIUS | repeat | — | Override circle radius for element. |
| `-ab` ELEMENT RADIUS | repeat | — | Override bond radius for element. |
| `-pb` EL1 EL2 RADIUS | repeat | — | Override bond cutoff for element pair. |
| `-ac` ELEMENT R G B | repeat | — | Override element colour (0–255). |
| `-cb` EL1 EL2 R G B | repeat | — | Override bond colour for element pair. |

### 3.15 Text overlay file (`-tf`)

One line per overlay. Fields: **start end text ix iy scale [transp]** (space-separated). Frame range is 0-based (start, end inclusive). Optional 7th field `transp` (0/1). If **text** is a special word, it is replaced by frame data:

| Text | Replaced by |
|------|-------------|
| ENERGY | Frame energy (e.g. from ReaxFF comment line). |
| ITERATION | Frame iteration/step. |
| VOLUME | Cell volume (requires cell in XYZ). |
| FRAME | "Frame <n>" (output frame number, 1-based). |

Example: `0 99 FRAME 10 20 2` shows "Frame 1", "Frame 2", … at pixel (10, 20) with scale 2 for frames 0–99.

### 3.16 Graphics overlay config (`-gf`)

One line per overlay window. Fields: **start end datafile ncols ix iy width height point_size col_x col_y**. Frame range start/end (0-based). **datafile** path is relative to the config file's directory unless absolute. **ncols** = number of columns in the data file; **ix, iy** = position of the window; **width, height** = window size in pixels; **point_size** = point size for plotting; **col_x, col_y** = column indices (0-based) for x and y data.

### 3.17 Config file (`-c`)

JSON or YAML file. Keys match CLI option names used in config (e.g. `size`, `cc`, `il`, `sh`, `bs`, `sc`, `sl`, `bc`, `bl`, `ic`, `nf`, `bd`, `supersample_scale`, `is_strain`, `sr`, `iv`, `vv`, `ia`, `as_arrow`, `ah`). Command-line arguments override config values. Example: `{"size": 800, "cc": 1, "il": 4}`.

---

## 4. Workflows

### 4.1 Single Snapshot

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -o snapshot.ppm
```

### 4.2 PNG Snapshot with Run Metadata

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -o snapshot.png --png --write-run-json
```

### 4.3 Full Trajectory as Numbered Frames

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -af 0 -ft 0 -o frames/frame.ppm
# Writes frames/0001.ppm, frames/0002.ppm, ...
```

### 4.4 Every 5th Frame as PNG

```bash
xmoltoppm -i traj.xyz -s 500 -cc 1 -il 3 -af 4 --png -o movie/m.png
# Writes movie/0001.png, movie/0002.png, ... (every 5th frame)
```

### 4.5 Encode Movie with ffmpeg

From PPM sequence:

```bash
ffmpeg -framerate 25 -i %04d.ppm -pix_fmt yuv420p movie.mp4
```

From PNG sequence:

```bash
ffmpeg -framerate 25 -i %04d.png -pix_fmt yuv420p movie.mp4
```

---

## 5. Programmatic Use

You can use the package without the CLI:

```python
from xmoltoppm.io.xyz_reader import read_xyz
from xmoltoppm.io.ppm_writer import write_ppm, write_png
from xmoltoppm.core.options import Options
from xmoltoppm.viz.render import render

# Read first frame
frames = list(read_xyz("traj.xyz"))
frame = frames[0]

# Default options
opts = Options(size=500, circle_style=1, line_style=3)

# Render to RGB array (H, W, 3) uint8
img = render(frame, opts)

# Write
write_ppm("out.ppm", img)
# or
write_png("out.png", img)
```

`Frame` holds: `nat`, `coords` (N×3), `symbols`, `cell` (3×3 or None), `molname`, `iteration`, `energy`, optional `velocities`, `estrain`.  
`Options` holds all rendering parameters. Field definitions are in `xmoltoppm.core.frame` and `xmoltoppm.core.options`.

---

## 6. Run Metadata JSON

With `--write-run-json`, a sidecar file is written (e.g. `output.ppm.run.json`) containing:

- `input_file` — Absolute path to the XYZ file.
- `version` — xmoltoppm version.
- `options` — Key options used (size, circle_style, line_style, etc.).
- `num_frames_rendered` — Number of frames written.

Use this to reproduce or document a run.

---

## 7. Troubleshooting

- **Comment line not 9 fields** — xmoltoppm uses defaults (empty name, zero iteration/energy, no cell) and still loads coordinates. Standard XYZ from other tools should work.
- **Image too small or clipped** — Increase `-s` (canvas size) or reduce `-bs` (border size).
- **Wrong or missing bonds** — Adjust `-ba` (bond radius adjustment). Bonds are drawn when distance < cutoff; cutoff depends on element radii and `-ba`.
- **All-frames output location** — With `-af`, `-o` is the path/directory for numbered files. Example: `-o frames/` writes `frames/0001.ppm`, `frames/0002.ppm`, … (with `-ft 0`). The base filename in `-o` is only used to derive the directory; output filenames are 0001, 0002, etc., or 1, 2, or molname depending on `-ft`.
- **PNG fails** — Install Pillow: `pip install Pillow`.
- **Boundary box (`-bb`)** — This port draws a 2D frame only (0.95× inset rectangle); no 3D box atoms or center dot. Requires cell data in the XYZ comment line.

---

## 8. Comparison with the Fortran Version

- **CLI only**: No interactive prompts; no `xmol2ppmtemp.dat`.
- **Robust XYZ**: Standard XYZ (any comment line) and 4-column atom lines work; 9-field header and extra columns are optional.
- **PNG**: Optional PNG output via Pillow.
- **Reproducibility**: Optional run JSON.
- **Implemented**: Strain/velocity colouring (`-is`, `-sr`, `-iv`, `-vv`, `-ut`), velocity arrows (`-ia`, `-as`, `-ah`), transparent atoms (`-t`, `-ti`), rotation animation (`-ri axis irend irstep`), text overlays (`-tx`, `-tf`; requires Pillow; `-tf` supports ENERGY/ITERATION/VOLUME/FRAME for frame data), graphics overlays (`-gf`), vibrational mode (`-vi NSTEP`, `-vs`, `-vr`). Core rendering (bonds + circles, styles 0–4) is implemented.

