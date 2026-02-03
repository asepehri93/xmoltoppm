# Contributing to xmoltoppm

Thanks for your interest in contributing. This document covers development setup, code layout, and how to add options or report issues.

## Development setup

1. Clone the repository.
2. Install in development mode (editable install):
   ```bash
   pip install -e .
   ```
   Optional PNG support: `pip install -e ".[png]"`.
3. Run tests:
   ```bash
   pytest tests/ -v
   ```

## Code layout

| Directory / module | Role |
|--------------------|------|
| `xmoltoppm/cli/main.py` | CLI entrypoint: argparse, single-frame vs all-frames, I/O wiring. |
| `xmoltoppm/core/options.py` | `Options` dataclass (all rendering parameters). |
| `xmoltoppm/core/frame.py` | `Frame` dataclass (atoms, coords, cell, etc.). |
| `xmoltoppm/core/elements.py` | Element table (radii, colours). |
| `xmoltoppm/io/xyz_reader.py` | Read XYZ trajectories. |
| `xmoltoppm/io/ppm_writer.py` | Write PPM and PNG. |
| `xmoltoppm/viz/render.py` | Render frame to RGB image. |
| `xmoltoppm/config_loader.py` | Load JSON/YAML config for CLI defaults. |

## Adding a new option

1. **`xmoltoppm/core/options.py`** — Add the field to the `Options` dataclass with a default.
2. **`xmoltoppm/cli/main.py`** — Add an argument (e.g. `parser.add_argument(...)`) and map it into `Options` (e.g. in `Options.from_cli(args)` or the post-parse block). If the option should be loadable from config, add the config key and dest to `_CONFIG_TO_ARGS`.
3. **`xmoltoppm/viz/render.py`** — Use the option where needed (e.g. `getattr(options, "option_name", default)`).
4. **`docs/guide.md`** — Document the option in the appropriate CLI reference table (§3) and in workflows or file-format sections if relevant.

## Publishing to PyPI

The package is **not** on PyPI until you publish it. After the first publish, `pip install xmoltoppm` and the PyPI badge in the README will work.

### Option A: Publish via GitHub Release (recommended)

1. Create a [PyPI account](https://pypi.org/account/register/) and a [PyPI API token](https://pypi.org/manage/account/token/) (scope: entire account or just this project).
2. In your GitHub repo: **Settings → Secrets and variables → Actions** → **New repository secret**. Name: `PYPI_API_TOKEN`, value: your token (prefix with `pypi-` if PyPI shows it that way).
3. Create a release: **Releases → Create a new release** → choose tag (e.g. `v1.0.0`, must match `version` in `pyproject.toml` for that release) → publish. The workflow [`.github/workflows/publish.yml`](.github/workflows/publish.yml) will build and upload to PyPI.

### Option B: Publish manually

From the project root:

```bash
pip install build twine
python -m build
twine check dist/*
twine upload dist/*
```

Use your PyPI username and the API token as the password. To test first: `twine upload --repository testpypi dist/*`.

## Reporting issues

Open an issue on the project repository. For bugs, include the command you ran, the XYZ input (or a minimal example), and the expected vs actual behaviour. For feature requests, a short description and use case helps.
