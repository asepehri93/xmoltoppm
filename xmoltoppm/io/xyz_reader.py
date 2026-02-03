"""XYZ trajectory reader (9-field comment + flexible atom lines).

Standard XYZ: line 1 = nat, line 2 = optional comment, then nat lines of
symbol x y z [optional columns]. If the second line has 9 fields
(molname, iteration, energy, a, b, c, alpha, beta, gamma), they are parsed
as ReaxFF-style header and cell; otherwise defaults are used and coordinates
are still loaded. Atom lines require at least 4 columns (symbol, x, y, z);
optional columns 5-6 (e.g. type, estrain) and 7-9 (vx, vy, vz) are read when present.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterator, Optional, Tuple

import numpy as np

from xmoltoppm.core.frame import Frame


def _parse_comment_line(line: str) -> Tuple[str, int, float, Optional[np.ndarray]]:
    """Parse second XYZ line: molname, iteration, energy, a,b,c, alpha,beta,gamma.
    Returns (molname, iteration, energy, cell_3x3 or None).
    """
    line = line.strip()
    parts = line.split()
    # Need at least 9 fields: molname, iteration, energy, a, b, c, alpha, beta, gamma
    if len(parts) >= 9:
        try:
            molname = parts[0]
            iteration = int(parts[1])
            energy = float(parts[2])
            a, b, c = float(parts[3]), float(parts[4]), float(parts[5])
            alpha, beta, gamma = float(parts[6]), float(parts[7]), float(parts[8])
            cell = _cell_matrix(a, b, c, alpha, beta, gamma)
            return molname, iteration, energy, cell
        except (ValueError, IndexError):
            pass
    return "" if not parts else parts[0], 0, 0.0, None


def _cell_matrix(a: float, b: float, c: float, alpha: float, beta: float, gamma: float) -> np.ndarray:
    """Build 3x3 cell matrix from lengths and angles (degrees). Same as Fortran tm11..tm33."""
    d2r = np.pi / 180.0
    ha, hb, hg = alpha * d2r, beta * d2r, gamma * d2r
    sa, ca = np.sin(ha), np.cos(ha)
    sb, cb = np.sin(hb), np.cos(hb)
    cphi = (np.cos(hg) - ca * cb) / (sa * sb) if sa * sb != 0 else 0.0
    cphi = max(-1.0, min(1.0, cphi))
    sphi = np.sqrt(1.0 - cphi * cphi)
    # Fortran: tm11=axis1*sinbet*sinphi, tm21=axis1*sinbet*cosphi, tm31=axis1*cosbet
    #          tm22=axis2*sinalf, tm32=axis2*cosalf, tm33=axis3
    cell = np.zeros((3, 3))
    cell[0, 0] = a * sb * sphi
    cell[0, 1] = a * sb * cphi
    cell[0, 2] = a * cb
    cell[1, 1] = b * sa
    cell[1, 2] = b * ca
    cell[2, 2] = c
    return cell


def read_xyz(path: str | Path) -> Iterator[Frame]:
    """Read XYZ file; yield one Frame per trajectory frame."""
    path = Path(path)
    with open(path, "r") as f:
        while True:
            nat_line = f.readline()
            if not nat_line.strip():
                break
            try:
                nat = int(nat_line.strip())
            except ValueError:
                break
            comment = f.readline()
            molname, iteration, energy, cell = _parse_comment_line(comment)
            symbols: list[str] = []
            coords_list: list[list[float]] = []
            estrain_list: list[float] = []
            vel_list: list[list[float]] = []
            for _ in range(nat):
                line = f.readline()
                if not line:
                    raise ValueError("Unexpected end of file in XYZ")
                tokens = line.split()
                if len(tokens) < 4:
                    raise ValueError("Atom line must have at least symbol x y z")
                symbols.append(tokens[0].strip())
                coords_list.append([float(tokens[1]), float(tokens[2]), float(tokens[3])])
                e = None
                if len(tokens) >= 6:
                    try:
                        e = float(tokens[5])
                    except ValueError:
                        pass
                estrain_list.append(e if e is not None else 0.0)
                v = None
                if len(tokens) >= 9:
                    try:
                        v = [float(tokens[6]), float(tokens[7]), float(tokens[8])]
                    except ValueError:
                        pass
                vel_list.append(v)

            coords = np.array(coords_list, dtype=np.float64)
            estrain = np.array(estrain_list, dtype=np.float64)
            velocities = None
            if all(v is not None for v in vel_list):
                velocities = np.array(vel_list, dtype=np.float64)
            yield Frame(
                nat=nat,
                coords=coords,
                symbols=symbols,
                cell=cell,
                molname=molname,
                iteration=iteration,
                energy=energy,
                velocities=velocities,
                estrain=estrain,
            )
