"""Element colors, radii, and bond cutoffs (from Fortran tables)."""

from __future__ import annotations

from typing import Dict, Optional, Tuple

# Fortran: ratyc = circle radius, raty = bond radius, ntyred/green/blue, amas
# Unknown -> C (index 6)
_DEFAULT = {
    "radius_circle": 1.80,
    "radius_bond": 1.80,
    "rgb": (100, 100, 100),
    "mass": 12.011,
}

_ELEMENTS: Dict[str, Dict[str, object]] = {
    "H": {"radius_circle": 1.00, "radius_bond": 1.00, "rgb": (250, 250, 250), "mass": 1.008},
    "He": {"radius_circle": 1.20, "radius_bond": 1.20, "rgb": (180, 180, 255), "mass": 4.002},
    "Li": {"radius_circle": 2.20, "radius_bond": 2.20, "rgb": (180, 100, 240), "mass": 6.941},
    "B": {"radius_circle": 1.60, "radius_bond": 1.60, "rgb": (160, 255, 160), "mass": 10.811},
    "C": {"radius_circle": 1.80, "radius_bond": 1.80, "rgb": (100, 100, 100), "mass": 12.011},
    "N": {"radius_circle": 1.65, "radius_bond": 1.65, "rgb": (0, 191, 255), "mass": 14.007},
    "O": {"radius_circle": 1.50, "radius_bond": 1.50, "rgb": (255, 50, 50), "mass": 15.9994},
    "F": {"radius_circle": 1.70, "radius_bond": 1.70, "rgb": (200, 50, 200), "mass": 18.998},
    "Na": {"radius_circle": 2.00, "radius_bond": 2.00, "rgb": (90, 120, 220), "mass": 22.990},
    "Mg": {"radius_circle": 1.90, "radius_bond": 1.90, "rgb": (100, 100, 250), "mass": 24.305},
    "Al": {"radius_circle": 1.90, "radius_bond": 1.90, "rgb": (200, 200, 250), "mass": 26.982},
    "Si": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (200, 200, 100), "mass": 28.086},
    "P": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (200, 0, 200), "mass": 30.974},
    "S": {"radius_circle": 2.40, "radius_bond": 2.40, "rgb": (255, 255, 0), "mass": 32.066},
    "Cl": {"radius_circle": 2.40, "radius_bond": 2.40, "rgb": (0, 255, 0), "mass": 35.453},
    "Ar": {"radius_circle": 2.40, "radius_bond": 2.40, "rgb": (50, 220, 120), "mass": 39.948},
    "K": {"radius_circle": 2.00, "radius_bond": 2.00, "rgb": (220, 220, 150), "mass": 39.098},
    "Ca": {"radius_circle": 2.40, "radius_bond": 2.40, "rgb": (100, 150, 200), "mass": 40.078},
    "Sc": {"radius_circle": 2.75, "radius_bond": 2.75, "rgb": (140, 200, 200), "mass": 44.956},
    "Ti": {"radius_circle": 2.75, "radius_bond": 2.75, "rgb": (100, 180, 220), "mass": 47.88},
    "V": {"radius_circle": 2.80, "radius_bond": 2.80, "rgb": (120, 150, 200), "mass": 50.942},
    "Cr": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (200, 125, 230), "mass": 51.996},
    "Mn": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (150, 225, 100), "mass": 54.939},
    "Fe": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (255, 125, 50), "mass": 55.847},
    "Co": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (150, 150, 255), "mass": 58.933},
    "Ni": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (120, 180, 180), "mass": 58.69},
    "Cu": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (200, 120, 100), "mass": 63.55},
    "Zn": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (110, 140, 210), "mass": 65.39},
    "Se": {"radius_circle": 2.50, "radius_bond": 2.50, "rgb": (230, 230, 40), "mass": 78.96},
    "Kr": {"radius_circle": 2.60, "radius_bond": 2.60, "rgb": (150, 200, 220), "mass": 83.80},
    "Rb": {"radius_circle": 2.10, "radius_bond": 2.10, "rgb": (180, 180, 130), "mass": 85.468},
    "Y": {"radius_circle": 3.00, "radius_bond": 3.00, "rgb": (180, 120, 120), "mass": 88.906},
    "Zr": {"radius_circle": 3.00, "radius_bond": 3.00, "rgb": (180, 180, 120), "mass": 91.224},
    "Mo": {"radius_circle": 3.00, "radius_bond": 3.00, "rgb": (180, 120, 180), "mass": 95.94},
    "Ru": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (80, 130, 210), "mass": 101.07},
    "Rh": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (180, 30, 250), "mass": 102.91},
    "Pd": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (210, 220, 250), "mass": 106.42},
    "Sb": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (130, 180, 50), "mass": 121.75},
    "Te": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (180, 130, 120), "mass": 127.60},
    "Xe": {"radius_circle": 3.00, "radius_bond": 1.00, "rgb": (140, 230, 150), "mass": 131.29},
    "Ba": {"radius_circle": 3.00, "radius_bond": 3.00, "rgb": (120, 180, 120), "mass": 137.327},
    "W": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (150, 20, 220), "mass": 195.08},
    "Pt": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (120, 120, 220), "mass": 195.08},
    "Au": {"radius_circle": 2.70, "radius_bond": 2.70, "rgb": (255, 220, 0), "mass": 196.97},
    "Bi": {"radius_circle": 2.80, "radius_bond": 2.80, "rgb": (100, 220, 150), "mass": 208.98},
    "Bo": {"radius_circle": 0.00, "radius_bond": 1000000.0, "rgb": (0, 0, 0), "mass": 0.10},
    "D": {"radius_circle": 1.00, "radius_bond": 1.00, "rgb": (0, 250, 250), "mass": 2.00},
}


def _norm(symbol: str) -> str:
    return symbol.strip()


def get_element(symbol: str) -> Dict[str, object]:
    """Return radius_circle, radius_bond, rgb (tuple of 3 ints), mass. Unknown -> C."""
    s = _norm(symbol)
    return _ELEMENTS.get(s, _DEFAULT).copy()


def get_radius_circle(
    symbol: str,
    override: Optional[Dict[str, float]] = None,
) -> float:
    """Circle radius for element; override from -ar if provided."""
    s = _norm(symbol)
    if override is not None and s in override:
        return override[s]
    return float(get_element(symbol)["radius_circle"])


def get_radius_bond(
    symbol: str,
    override: Optional[Dict[str, float]] = None,
) -> float:
    """Bond radius for element; override from -ab if provided."""
    s = _norm(symbol)
    if override is not None and s in override:
        return override[s]
    return float(get_element(symbol)["radius_bond"])


def get_bond_cutoff(
    sym1: str,
    sym2: str,
    pair_override: Optional[Dict[Tuple[str, str], float]] = None,
    element_bond_override: Optional[Dict[str, float]] = None,
) -> float:
    """Bond distance cutoff = 0.5 * (r1 + r2) or pair override (-pb). Bo-other -> 0."""
    s1, s2 = _norm(sym1), _norm(sym2)
    if pair_override:
        key = (s1, s2)
        if key in pair_override:
            return pair_override[key]
        key = (s2, s1)
        if key in pair_override:
            return pair_override[key]
    r1 = get_radius_bond(sym1, element_bond_override)
    r2 = get_radius_bond(sym2, element_bond_override)
    if (s1 == "Bo" and s2 != "Bo") or (s2 == "Bo" and s1 != "Bo"):
        return 0.0
    if s1 == "Bo" and s2 == "Bo":
        return 0.0
    return 0.5 * (r1 + r2)


def get_mass(symbol: str) -> float:
    """Return atomic mass (g/mol) for velocity/temperature coloring. Unknown -> C."""
    return float(get_element(symbol)["mass"])


def get_rgb(
    symbol: str,
    override: Optional[Dict[str, Tuple[int, int, int]]] = None,
) -> Tuple[int, int, int]:
    """Return (r, g, b) 0-255 for element; override from -ac if provided."""
    s = _norm(symbol)
    if override is not None and s in override:
        return override[s]
    d = get_element(symbol)
    return d["rgb"]  # type: ignore


def get_bond_rgb(
    sym1: str,
    sym2: str,
    override: Optional[Dict[Tuple[str, str], Tuple[int, int, int]]] = None,
) -> Optional[Tuple[int, int, int]]:
    """Bond color override (-cb); returns None if no override."""
    if not override:
        return None
    s1, s2 = _norm(sym1), _norm(sym2)
    return override.get((s1, s2)) or override.get((s2, s1))
