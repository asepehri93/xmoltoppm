# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Graphics from energy (`-gfe`)** — Build a graphics overlay from per-frame energies in the XYZ comment line (iteration, energy). Writes `<output_stem>.energy.txt` and plots it in a corner; use with `-af` so the overlay cursor matches the current frame. Optional `-gfe IX IY WIDTH HEIGHT` for overlay position/size.
- **Test** — `test_graphics_overlay_highlight_matches_frame_index` ensures the overlay highlight (cursor) changes with frame index.

### Fixed

- **Graphics overlay cursor** — In all-frames mode, the overlay now highlights the row corresponding to the **trajectory frame index** (not the output image index), so the cursor matches the frame being shown.

## [1.0.0] - 2025-01-30

### Added

- **Shadow inversion (`-si`)** — Fortran parity: option to invert directional shadow (line style 4). Use `-si -1` for inverted shadow (e.g. for projectors). Default: `1`.
- **Pillow warning for text overlays** — When `-tx` or `-tf` is used but Pillow is not installed, a one-time warning is printed to stderr and text is skipped.
- **Graphics overlay path resolution** — Data file paths in `-gf` config files are resolved relative to the config file’s directory when not absolute.
- **`--max-frames` / `-mf N`** — Limit rendering to the first N frames. Useful for large trajectories.
- **New tests** — Vibrational mode (displacement changes output), graphics overlay (-gf), text overlay (with/without Pillow), invalid CLI option (-si 0).

### Changed

- **Version** — Bumped from 0.1.0 to 1.0.0 for first stable release.

See [README.md](README.md) and [docs/guide.md](docs/guide.md) for full option list and usage.
