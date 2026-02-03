"""Load options from YAML or JSON config file."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict


def load_config(path: str | Path) -> Dict[str, Any]:
    """Load a config file (JSON or YAML). Returns a flat dict of option names -> values."""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")
    suffix = path.suffix.lower()
    with open(path, "r") as f:
        if suffix == ".json":
            data = json.load(f)
        elif suffix in (".yaml", ".yml"):
            try:
                import yaml
            except ImportError:
                raise ImportError("YAML config requires PyYAML: pip install PyYAML") from None
            data = yaml.safe_load(f)
            if data is None:
                data = {}
        else:
            raise ValueError(f"Config file must be .json or .yaml/.yml, got {suffix}")
    out: Dict[str, Any] = {}
    # Top-level "plane": {x: 10, z: 5} -> plane_x, plane_z
    if "plane" in data and isinstance(data["plane"], dict):
        for ax in ("x", "y", "z"):
            if ax in data["plane"]:
                out[f"plane_{ax}"] = float(data["plane"][ax])
        rest = {k: v for k, v in data.items() if k != "plane"}
    else:
        rest = data
    out.update(_flatten_config(rest))
    return out


def _flatten_config(data: Dict[str, Any], prefix: str = "") -> Dict[str, Any]:
    """Flatten nested dict to one level; keys match CLI/Options attribute names."""
    out: Dict[str, Any] = {}
    for key, value in data.items():
        if isinstance(value, dict) and not _is_option_value(value):
            out.update(_flatten_config(value, prefix=f"{prefix}{key}_"))
        else:
            name = key if not prefix else prefix + key
            out[name] = value
    return out


def _is_option_value(v: Any) -> bool:
    if not isinstance(v, dict):
        return True
    return False
