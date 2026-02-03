"""2D projection renderer: Frame + Options -> RGB array."""

from __future__ import annotations

import math
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

from xmoltoppm.core.elements import (
    get_bond_cutoff,
    get_bond_rgb,
    get_element,
    get_radius_bond,
    get_radius_circle,
    get_rgb,
)
from xmoltoppm.core.frame import Frame
from xmoltoppm.core.options import Options

_DEG2RAD = math.pi / 180.0
_text_overlay_warned = False


def centroid_after_plane_filter(frame: Frame, options: Options) -> np.ndarray:
    """Return centroid of frame after applying plane filter (for ref_center when -ic 0)."""
    coords = frame.coords.copy()
    plane_x, plane_y, plane_z = options.plane_x, options.plane_y, options.plane_z
    if plane_x < 1e9 or plane_y < 1e9 or plane_z < 1e9:
        keep = (
            (coords[:, 0] <= plane_x)
            & (coords[:, 1] <= plane_y)
            & (coords[:, 2] <= plane_z)
        )
        coords = coords[keep]
    if coords.shape[0] == 0:
        return np.zeros(3, dtype=np.float64)
    return coords.mean(axis=0)


def render(
    frame: Frame,
    options: Options,
    ref_center: Optional[np.ndarray] = None,
    frame_index: int = 0,
    vibrational_sub_step: Optional[int] = None,
) -> np.ndarray:
    """Render one frame to RGB array (H, W, 3) uint8.

    ref_center: when center_mode==0 (keep first-frame center), pass first frame
    centroid (3,) so all frames use the same center; else None.
    frame_index: 0-based output frame index; used with -ri to add frame_index * step
    to the rotation around the animation axis.
    vibrational_sub_step: when -vi nstep is set, 0..nstep-1; displace coords along
    velocity/mode vector (forward then back) for animation.
    """
    coords = frame.coords.copy()
    symbols = list(frame.symbols)
    nat = frame.nat
    velocities = frame.velocities
    estrain = frame.estrain

    # Vibrational mode (-vi): displace coords along velocity/mode vector for this sub-step
    nfreq = getattr(options, "vibrational_steps", 0)
    vscale = getattr(options, "vibrational_scale", 1.0)
    if nfreq > 0 and velocities is not None and vibrational_sub_step is not None and vibrational_sub_step > 0:
        k = vibrational_sub_step  # 0-based
        if k < nfreq // 2:
            factor = k
        else:
            factor = nfreq - 1 - k
        coords = coords + factor * vscale * velocities
    # 1-based original indices for transparent_atoms (-t, -ti)
    original_1based = np.arange(1, nat + 1, dtype=np.intp)

    # Plane filter: remove atoms with coord > plane (Fortran -pp)
    plane_x, plane_y, plane_z = options.plane_x, options.plane_y, options.plane_z
    if plane_x < 1e9 or plane_y < 1e9 or plane_z < 1e9:
        keep = (
            (coords[:, 0] <= plane_x)
            & (coords[:, 1] <= plane_y)
            & (coords[:, 2] <= plane_z)
        )
        coords = coords[keep]
        symbols = [s for s, k in zip(symbols, keep) if k]
        original_1based = original_1based[keep]
        nat = len(symbols)
        if velocities is not None:
            velocities = velocities[keep]
        if estrain is not None:
            estrain = estrain[keep]
    # -at: change atom types (1-based index -> symbol)
    for idx, sym in (getattr(options, "atom_type_changes", None) or {}).items():
        i0 = int(idx) - 1
        if 0 <= i0 < len(symbols):
            symbols[i0] = sym.strip()
    nat = len(symbols)

    # -na: only show atoms within radius_radius of atom radius_center_atom (1-based)
    rca = getattr(options, "radius_center_atom", -1)
    if rca > 0 and nat > 0:
        rrad = getattr(options, "radius_radius", 1e9)
        i0 = rca - 1
        if 0 <= i0 < nat:
            cx, cy, cz = coords[i0, 0], coords[i0, 1], coords[i0, 2]
            keep = np.zeros(nat, dtype=bool)
            for i in range(nat):
                d = math.sqrt(
                    (coords[i, 0] - cx) ** 2
                    + (coords[i, 1] - cy) ** 2
                    + (coords[i, 2] - cz) ** 2
                )
                keep[i] = d <= rrad
            coords = coords[keep]
            symbols = [s for s, k in zip(symbols, keep) if k]
            original_1based = original_1based[keep]
            nat = len(symbols)
            if velocities is not None:
                velocities = velocities[keep]
            if estrain is not None:
                estrain = estrain[keep]
    # -bb: boundary box is drawn later as a 2D frame only (_draw_image_frame); no 3D box atoms
    nat_prev = nat
    draw_bb = getattr(options, "draw_boundary_box", 0)
    if nat == 0:
        size = max(1, options.size)
        return np.full((size, size, 3), options.background, dtype=np.uint8)

    # Supersample
    scale = max(1, getattr(options, "supersample_scale", 1))
    size = options.size * scale

    # Center: use ref_center when -ic 0, else centroid of all atoms
    if ref_center is not None:
        center = ref_center
    else:
        center = coords.mean(axis=0)
    coords = coords - center

    # -r: rotation (degrees) around x, then y, then z; -ri: add frame_index * step to one axis
    rx = getattr(options, "rotation_rx", 0.0)
    ry = getattr(options, "rotation_ry", 0.0)
    rz = getattr(options, "rotation_rz", 0.0)
    raxis = (getattr(options, "rotation_animation_axis", "") or "").strip().lower()
    rstep = getattr(options, "rotation_animation_step", 0.0)
    if raxis and rstep != 0:
        extra = frame_index * rstep
        if raxis == "x":
            rx += extra
        elif raxis == "y":
            ry += extra
        elif raxis == "z":
            rz += extra
    rx, ry, rz = rx * _DEG2RAD, ry * _DEG2RAD, rz * _DEG2RAD
    if abs(rx) > 1e-10:
        c, s = math.cos(rx), math.sin(rx)
        y, z = coords[:, 1].copy(), coords[:, 2].copy()
        coords[:, 1] = y * c - z * s
        coords[:, 2] = y * s + z * c
    if abs(ry) > 1e-10:
        c, s = math.cos(ry), math.sin(ry)
        x, z = coords[:, 0].copy(), coords[:, 2].copy()
        coords[:, 0] = x * c + z * s
        coords[:, 2] = -x * s + z * c
    if abs(rz) > 1e-10:
        c, s = math.cos(rz), math.sin(rz)
        x, y = coords[:, 0].copy(), coords[:, 1].copy()
        coords[:, 0] = x * c - y * s
        coords[:, 1] = x * s + y * c

    # Apply same rotation to velocities for arrow drawing
    if velocities is not None:
        v = velocities.copy()
        if abs(rx) > 1e-10:
            c, s = math.cos(rx), math.sin(rx)
            y, z = v[:, 1].copy(), v[:, 2].copy()
            v[:, 1] = y * c - z * s
            v[:, 2] = y * s + z * c
        if abs(ry) > 1e-10:
            c, s = math.cos(ry), math.sin(ry)
            x, z = v[:, 0].copy(), v[:, 2].copy()
            v[:, 0] = x * c + z * s
            v[:, 2] = -x * s + z * c
        if abs(rz) > 1e-10:
            c, s = math.cos(rz), math.sin(rz)
            x, y = v[:, 0].copy(), v[:, 1].copy()
            v[:, 0] = x * c - y * s
            v[:, 1] = x * s + y * c
        velocities = v

    # Sort by z descending
    order = np.argsort(-coords[:, 2])
    coords = coords[order]
    symbols_sorted = [symbols[i] for i in order]
    original_1based_sorted = original_1based[order]
    estrain_sorted = estrain[order] if estrain is not None else None
    velocities_sorted = velocities[order] if velocities is not None else None

    # Viewport extent: when drawing boundary box, scale to molecule only so atoms are not tiny
    draw_bb = getattr(options, "draw_boundary_box", 0)
    if draw_bb and nat_prev < nat:
        mol_mask = original_1based_sorted != 0  # use sorted: coords are in z-order
        if np.any(mol_mask):
            xmin = coords[mol_mask, 0].min()
            xmax = coords[mol_mask, 0].max()
            ymin = coords[mol_mask, 1].min()
            ymax = coords[mol_mask, 1].max()
            zmin = coords[mol_mask, 2].min()
            zmax = coords[mol_mask, 2].max()
        else:
            xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
            ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
            zmin, zmax = coords[:, 2].min(), coords[:, 2].max()
    else:
        xmin, xmax = coords[:, 0].min(), coords[:, 0].max()
        ymin, ymax = coords[:, 1].min(), coords[:, 1].max()
        zmin, zmax = coords[:, 2].min(), coords[:, 2].max()
    zdif2 = max(zmax - zmin, 1e-4)

    brdersize = options.border_size
    # Fixed canvas (square) so every frame has the same dimensions for movies
    dxmax = 2.0 * brdersize + (xmax - xmin)
    dymax = 2.0 * brdersize + (ymax - ymin)
    extent_max = max(dxmax, dymax, 1e-10)
    nsizex = size
    nsizey = size
    vconv = size / extent_max

    xmiddle = 0.5 * (xmax + xmin)
    ymiddle = 0.5 * (ymax + ymin)
    vps = getattr(options, "perspective_system", 0.0)

    def to_pixel(x: float, y: float, z: float) -> Tuple[int, int]:
        pxh = 0.5 * nsizex + vconv * (x - xmiddle)
        pyh = 0.5 * nsizey + vconv * (y - ymiddle)
        if vps != 0:
            rper = 1.0 + vps * (-1.0 + (z - zmin) / zdif2)
            pxh = 0.5 * nsizex * (1.0 - rper) + pxh * rper
            pyh = 0.5 * nsizey * (1.0 - rper) + pyh * rper
        return int(pxh), int(pyh)

    rgb = np.full((nsizey, nsizex, 3), options.background, dtype=np.uint8)

    px_list: List[int] = []
    py_list: List[int] = []
    for i in range(nat):
        px, py = to_pixel(coords[i, 0], coords[i, 1], coords[i, 2])
        px_list.append(px)
        py_list.append(py)

    # Per-atom "temperature" for velocity coloring (-iv 1: from estrain as speed, -iv 2: from |v|, -ut: estrain as T)
    velocity_coloring = getattr(options, "velocity_coloring", 0)
    velocity_reference_temp = getattr(options, "velocity_reference_temp", 300.0)
    use_temp_range = getattr(options, "use_temp_range", 0)
    tempat_sorted: Optional[np.ndarray] = None
    if velocity_coloring > 0 or use_temp_range:
        tempat_sorted = np.zeros(nat, dtype=np.float64)
        for i in range(nat):
            sym = symbols_sorted[i]
            mass = get_mass(sym)
            if velocity_coloring == 1 and estrain_sorted is not None:
                # T from 0.5*m*v^2; estrain used as speed (Fortran: estrain(i1) as v)
                v = float(estrain_sorted[i])
                tempat_sorted[i] = 2.0 * (0.5 * mass * v * v / _CONVMD) / (3.0 * _RGASC * _XJOUCA / 1e3)
            elif velocity_coloring == 2 and velocities_sorted is not None:
                vx, vy, vz = velocities_sorted[i, 0], velocities_sorted[i, 1], velocities_sorted[i, 2]
                v = math.sqrt(vx * vx + vy * vy + vz * vz)
                tempat_sorted[i] = 2.0 * (0.5 * mass * v * v / _CONVMD) / (3.0 * _RGASC * _XJOUCA / 1e3)
            elif use_temp_range and estrain_sorted is not None:
                tempat_sorted[i] = float(estrain_sorted[i])
            else:
                tempat_sorted[i] = velocity_reference_temp

    vline = options.line_size * vconv
    ba = options.bond_radius_adjust
    vcirc = options.circle_size
    bd = options.bond_thickness_by_distance
    er_circle = getattr(options, "element_radius_circle", None) or {}
    er_bond = getattr(options, "element_radius_bond", None) or {}
    er_rgb = getattr(options, "element_rgb", None) or {}
    pair_cutoff = getattr(options, "pair_bond_cutoff", None) or {}
    pair_rgb = getattr(options, "pair_bond_rgb", None) or {}
    ignore_bonds = getattr(options, "ignore_bonds", None) or []

    def bond_ignored(oi: int, oj: int) -> bool:
        a1, a2 = oi + 1, oj + 1
        for (i1, i2, n1, n2) in ignore_bonds:
            if (i1 <= a1 <= i2 and n1 <= a2 <= n2) or (n1 <= a1 <= n2 and i1 <= a2 <= i2):
                return True
        return False

    # Box edges (pairs of corner indices 0..7 for 12 edges)
    _BOX_EDGES = {(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)}

    def is_box_edge(oi: int, oj: int) -> bool:
        if symbols_sorted[oi] != "Bo" or symbols_sorted[oj] != "Bo":
            return False
        bi = order[oi] - nat_prev if order[oi] >= nat_prev else None
        bj = order[oj] - nat_prev if order[oj] >= nat_prev else None
        if bi is None or bj is None:
            return False
        return (min(bi, bj), max(bi, bj)) in _BOX_EDGES

    # Draw bonds
    for i in range(nat):
        for j in range(i + 1, nat):
            if bond_ignored(i, j):
                continue
            if is_box_edge(i, j):
                # Box edges not drawn here when draw_bb (drawn as image frame later)
                continue
            cutoff = ba * get_bond_cutoff(
                symbols_sorted[i], symbols_sorted[j],
                pair_override=pair_cutoff,
                element_bond_override=er_bond,
            )
            if cutoff <= 0:
                continue
            d = math.sqrt(
                (coords[i, 0] - coords[j, 0]) ** 2
                + (coords[i, 1] - coords[j, 1]) ** 2
                + (coords[i, 2] - coords[j, 2]) ** 2
            )
            if d >= cutoff:
                continue
            if bd:
                r1 = get_radius_bond(symbols_sorted[i], er_bond)
                r2 = get_radius_bond(symbols_sorted[j], er_bond)
                qdis = 0.5 * (r1 + r2) / max(d, 1e-10)
                line_width = int(max(1, vline * qdis))
            else:
                line_width = int(max(1, vline))
            x0, y0 = px_list[i], py_list[i]
            x1, y1 = px_list[j], py_list[j]
            cb = get_bond_rgb(symbols_sorted[i], symbols_sorted[j], pair_rgb)
            if cb is not None:
                c0, c1 = cb, cb
            else:
                c0 = get_rgb(symbols_sorted[i], er_rgb)
                c1 = get_rgb(symbols_sorted[j], er_rgb)
            _draw_line(
                rgb, x0, y0, x1, y1,
                c0, c1,
                options.line_style, options.shadow, options.shadow_strength,
                options.shadow_direction, coords[i, 2], coords[j, 2],
                zmin, zdif2, line_width, options.draw_boundary_lines,
                getattr(options, "shadow_invert", 1),
            )

    # Draw atoms (circles)
    strain_coloring = getattr(options, "strain_coloring", 0)
    strain_scale = getattr(options, "strain_scale", 1.0)
    temp_min = getattr(options, "temp_min", 0.0)
    temp_max = getattr(options, "temp_max", 1000.0)
    transparent_atoms = getattr(options, "transparent_atoms", frozenset()) or frozenset()

    for i in range(nat):
        px, py = px_list[i], py_list[i]
        r, g, b = get_rgb(symbols_sorted[i], er_rgb)

        # Strain coloring (-is): red += scale*estrain, green -= scale*estrain
        if strain_coloring and estrain_sorted is not None:
            delta = int(strain_scale * float(estrain_sorted[i]))
            r = min(255, max(0, r + delta))
            g = min(255, max(0, g - delta))

        # Velocity coloring (-iv): r,g,b shift by (T - Tref)/4
        if velocity_coloring > 0 and tempat_sorted is not None:
            delta = int((float(tempat_sorted[i]) - velocity_reference_temp) / 4.0)
            r = min(255, max(0, r + delta))
            g = min(255, max(0, g + delta))
            b = min(255, max(0, b + delta))

        # Use temp range (-ut): gradient red (cold) -> blue (hot), green=0
        if use_temp_range and tempat_sorted is not None:
            temp_diff = max(temp_max - temp_min, 1e-10)
            hulp1 = (temp_max - float(tempat_sorted[i])) / temp_diff
            hulp2 = (float(tempat_sorted[i]) - temp_min) / temp_diff
            r = int(255 - 255 * hulp1)
            b = int(255 - 255 * hulp2)
            g = 0
            r = min(255, max(0, r))
            b = min(255, max(0, b))

        radius_ang = vcirc * get_radius_circle(symbols_sorted[i], er_circle)
        radius_pix = max(1, int(radius_ang * vconv))
        z = coords[i, 2]
        vshad = 1.0 + options.shadow * (-1.0 + (z - zmin) / zdif2)
        if options.shadow > 0:
            vshad = max(0.3, min(1.0, vshad))
        is_transparent = int(original_1based_sorted[i]) in transparent_atoms
        _draw_circle(
            rgb, px, py, radius_pix,
            (r, g, b), options.circle_style, z, zmin, zdif2, vshad,
            options.draw_boundary_lines,
            transparent=is_transparent,
        )

    # Draw velocity arrows (-ia)
    draw_arrows = getattr(options, "draw_arrows", 0)
    if draw_arrows and velocities_sorted is not None:
        arrow_scale = getattr(options, "arrow_scale", 1.0)
        arrow_head = max(2, getattr(options, "arrow_head_size", 10))
        for i in range(nat):
            if symbols_sorted[i] == "Bo":
                continue
            vx, vy = float(velocities_sorted[i, 0]), float(velocities_sorted[i, 1])
            if vx == 0 and vy == 0:
                continue
            px, py = px_list[i], py_list[i]
            ex = px + vconv * arrow_scale * vx
            ey = py + vconv * arrow_scale * vy
            r, g, b = get_rgb(symbols_sorted[i], er_rgb)
            _draw_arrow(rgb, px, py, int(round(ex)), int(round(ey)), (r, g, b), arrow_head)

    # Text overlays (-tx, -tf): draw after arrows, before downscale
    text_overlays_list = getattr(options, "text_overlays", None) or []
    text_file_entries_list = getattr(options, "text_file_entries", None) or []
    to_draw: List[Tuple[str, int, int, int]] = list(text_overlays_list)
    for (start, end, special, text, ix, iy, scl, _transp) in text_file_entries_list:
        if start <= frame_index <= end:
            if special == 1:
                text = f"ENERGY:{frame.energy:.2f} KCAL/MOL"
            elif special == 2:
                text = f"ITERATION:{frame.iteration}"
            elif special == 3 and frame.cell is not None:
                vol = abs(float(np.linalg.det(frame.cell)))
                text = f"VOLUME:{vol:.2f} ANGSTROM"
            elif special == 4:
                text = f"Frame {frame_index + 1}"
            to_draw.append((text, ix, iy, max(1, scl)))
    if to_draw:
        _draw_text_overlays(rgb, to_draw)

    # Graphics overlays (-gf): plot data in windows
    graphics_entries = getattr(options, "graphics_file_entries", None) or []
    for entry in graphics_entries:
        start, end, _path, _ncols, _ix, _iy, _w, _h, _psz, _cx, _cy = entry
        if start <= frame_index <= end:
            _draw_graphics_overlay(rgb, entry, frame_index)

    # Boundary box: draw as exact image frame (four edges)
    if draw_bb and frame.cell is not None:
        _draw_image_frame(rgb, nsizex, nsizey)

    if scale > 1:
        # Downsample: (H*scale, W*scale, 3) -> (H, W, 3) by block average
        h, w = rgb.shape[0], rgb.shape[1]
        h_out, w_out = h // scale, w // scale
        if h_out < 1:
            h_out = 1
        if w_out < 1:
            w_out = 1
        # Trim to multiple of scale
        rgb = rgb[: h_out * scale, : w_out * scale, :]
        rgb = (
            rgb.reshape(h_out, scale, w_out, scale, 3).mean(axis=(1, 3)).astype(np.uint8)
        )

    return rgb


def _draw_line(
    rgb: np.ndarray,
    x0: int, y0: int, x1: int, y1: int,
    color0: Tuple[int, int, int], color1: Tuple[int, int, int],
    line_style: int, shadow: float, shadow_strength: float,
    shadow_dir: str, z0: float, z1: float,
    zmin: float, zdif2: float, width: int,
    draw_boundary_lines: int,
    shadow_invert: int = 1,
) -> None:
    """Draw a thick line from (x0,y0) to (x1,y1)."""
    if x0 == x1 and y0 == y1:
        return  # avoid drawing a single dot for zero-length lines
    h, w = rgb.shape[0], rgb.shape[1]
    length = max(1, int(math.hypot(x1 - x0, y1 - y0)))
    dx = (x1 - x0) / length
    dy = (y1 - y0) / length
    perp_dx = -dy
    perp_dy = dx
    vwadj = math.sqrt(dx * dx + dy * dy)
    if vwadj < 1e-10:
        vwadj = 1.0
    half = max(0, width // 2)

    for step in range(length + 1):
        cx = int(x0 + step * dx)
        cy = int(y0 + step * dy)
        t = step / length if length > 0 else 0.5
        if line_style in (1, 3, 4):
            r = int(color0[0] * (1 - t) + color1[0] * t)
            g = int(color0[1] * (1 - t) + color1[1] * t)
            b = int(color0[2] * (1 - t) + color1[2] * t)
        else:
            r = (color0[0] + color1[0]) // 2
            g = (color0[1] + color1[1]) // 2
            b = (color0[2] + color1[2]) // 2
        if line_style >= 2 and shadow > 0:
            zmid = z0 * (1 - t) + z1 * t
            vshad = 1.0 + shadow * (-1.0 + (zmid - zmin) / zdif2)
            vshad = max(0.3, min(1.0, vshad))
            r = int(r * vshad)
            g = int(g * vshad)
            b = int(b * vshad)
        if line_style == 4 and shadow_strength != 0:
            dz = z1 - z0
            dlen = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2) + 1e-10
            if shadow_dir == "yz":
                vshadhu = (dz / dlen) * shadow_strength
            elif shadow_dir == "xz":
                vshadhu = (dz / dlen) * shadow_strength
            else:
                vshadhu = (dz / dlen) * shadow_strength
            vshadhu = shadow_invert * vshadhu
            r = int(max(0, min(255, r + vshadhu)))
            g = int(max(0, min(255, g + vshadhu)))
            b = int(max(0, min(255, b + vshadhu)))
        for dw in range(-half, half + 1):
            px = int(cx + perp_dx * dw)
            py = int(cy + perp_dy * dw)
            if 0 <= px < w and 0 <= py < h:
                if draw_boundary_lines == 0:
                    rgb[py, px, 0] = min(255, max(0, r))
                    rgb[py, px, 1] = min(255, max(0, g))
                    rgb[py, px, 2] = min(255, max(0, b))
                else:
                    rgb[py, px, 0] = min(255, max(0, r))
                    rgb[py, px, 1] = min(255, max(0, g))
                    rgb[py, px, 2] = min(255, max(0, b))


_FRAME_RATIO = 0.95  # boundary frame size as fraction of image (inset from edges)


def _draw_image_frame(rgb: np.ndarray, w: int, h: int) -> None:
    """Draw a black rectangle at _FRAME_RATIO of the image size (inset from edges)."""
    black = (0, 0, 0)
    inset = (1.0 - _FRAME_RATIO) / 2.0
    left = max(0, int(inset * w))
    right = min(w - 1, int((1.0 - inset) * w))
    top = max(0, int(inset * h))
    bottom = min(h - 1, int((1.0 - inset) * h))
    if left >= right or top >= bottom:
        return
    # top, right, bottom, left
    for (x0, y0, x1, y1) in [
        (left, top, right, top),
        (right, top, right, bottom),
        (right, bottom, left, bottom),
        (left, bottom, left, top),
    ]:
        _draw_line(
            rgb, x0, y0, x1, y1,
            black, black, 0, 0.0, 0.0, "yz", 0.0, 0.0, 1.0, 1, 1, 1,
        )


def _draw_circle(
    rgb: np.ndarray, cx: int, cy: int, radius_pix: int,
    color: Tuple[int, int, int], circle_style: int,
    z: float, zmin: float, zdif2: float, vshad: float,
    draw_boundary_lines: int,
    transparent: bool = False,
) -> None:
    """Draw a filled circle; circle_style 1 = shaded by distance from center.
    If transparent=True, draw only the outline (no fill).
    """
    h, w = rgb.shape[0], rgb.shape[1]
    r, g, b = color
    R_ang = radius_pix  # in pixel units
    outline_color = (0, 0, 0) if draw_boundary_lines == 0 else (r, g, b)
    for dy in range(-radius_pix, radius_pix + 1):
        for dx in range(-radius_pix, radius_pix + 1):
            dist = math.sqrt(dx * dx + dy * dy)
            if dist > R_ang:
                if dist <= R_ang + 1:
                    px, py = cx + dx, cy + dy
                    if 0 <= px < w and 0 <= py < h:
                        rgb[py, px, 0], rgb[py, px, 1], rgb[py, px, 2] = outline_color
                continue
            if transparent:
                continue  # no fill for transparent atoms
            px, py = cx + dx, cy + dy
            if px < 0 or px >= w or py < 0 or py >= h:
                continue
            ratio = 1.0
            if circle_style == 1 and R_ang > 0:
                ratioh = dist / R_ang
                ratio = 1.0 - ratioh * ratioh * ratioh
            nr = int(r * ratio * vshad)
            ng = int(g * ratio * vshad)
            nb = int(b * ratio * vshad)
            rgb[py, px, 0] = min(255, max(0, nr))
            rgb[py, px, 1] = min(255, max(0, ng))
            rgb[py, px, 2] = min(255, max(0, nb))


def _draw_arrow(
    rgb: np.ndarray,
    x0: int, y0: int, x1: int, y1: int,
    color: Tuple[int, int, int],
    head_size: int,
) -> None:
    """Draw an arrow from (x0,y0) to (x1,y1) with a triangular head."""
    h, w = rgb.shape[0], rgb.shape[1]
    r, g, b = color
    length = math.hypot(x1 - x0, y1 - y0)
    if length < 1e-6:
        return
    dx = (x1 - x0) / length
    dy = (y1 - y0) / length
    base_x = x1 - head_size * dx
    base_y = y1 - head_size * dy
    shaft_len = math.hypot(base_x - x0, base_y - y0)
    n_steps = max(1, int(shaft_len))
    for k in range(n_steps + 1):
        t = k / n_steps
        cx = int(x0 + t * (base_x - x0))
        cy = int(y0 + t * (base_y - y0))
        if 0 <= cx < w and 0 <= cy < h:
            rgb[cy, cx, 0], rgb[cy, cx, 1], rgb[cy, cx, 2] = r, g, b
    # Arrowhead: triangle at tip; tip at (x1,y1), base corners perpendicular
    perp_x = -dy
    perp_y = dx
    base_x = x1 - head_size * dx
    base_y = y1 - head_size * dy
    c1_x = int(base_x + perp_x * head_size * 0.6)
    c1_y = int(base_y + perp_y * head_size * 0.6)
    c2_x = int(base_x - perp_x * head_size * 0.6)
    c2_y = int(base_y - perp_y * head_size * 0.6)
    # Rasterize triangle (tip, c1, c2) with scanline or point-in-triangle
    for dy_ in range(-head_size - 1, head_size + 2):
        for dx_ in range(-head_size - 1, head_size + 2):
            px, py = x1 + dx_, y1 + dy_
            if px < 0 or px >= w or py < 0 or py >= h:
                continue
            # Barycentric: is (px,py) inside triangle (x1,y1), (c1_x,c1_y), (c2_x,c2_y)?
            v0x, v0y = c1_x - x1, c1_y - y1
            v1x, v1y = c2_x - x1, c2_y - y1
            v2x, v2y = px - x1, py - y1
            dot00 = v0x * v0x + v0y * v0y
            dot01 = v0x * v1x + v0y * v1y
            dot02 = v0x * v2x + v0y * v2y
            dot11 = v1x * v1x + v1y * v1y
            dot12 = v1x * v2x + v1y * v2y
            inv = dot00 * dot11 - dot01 * dot01
            if abs(inv) < 1e-10:
                continue
            inv = 1.0 / inv
            u = (dot11 * dot02 - dot01 * dot12) * inv
            v = (dot00 * dot12 - dot01 * dot02) * inv
            if u >= 0 and v >= 0 and u + v <= 1:
                rgb[py, px, 0], rgb[py, px, 1], rgb[py, px, 2] = r, g, b


def _draw_text_overlays(
    rgb: np.ndarray,
    to_draw: List[Tuple[str, int, int, int]],
) -> None:
    """Draw text overlays at (ix, iy) with scale. Requires Pillow."""
    global _text_overlay_warned
    if not to_draw:
        return
    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        if not _text_overlay_warned:
            print(
                "Warning: text overlays require Pillow (pip install Pillow); skipping text.",
                file=sys.stderr,
            )
            _text_overlay_warned = True
        return
    img = Image.fromarray(rgb, mode="RGB")
    draw = ImageDraw.Draw(img)
    default_font = ImageFont.load_default()
    font_size_candidates = [  # (path, size_key); try truetype then fall back to default
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/System/Library/Fonts/Helvetica.ttc",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ]
    for (text, ix, iy, scale) in to_draw:
        if not text:
            continue
        font_size = max(8, min(72, 8 * scale))
        font = default_font
        for path in font_size_candidates:
            try:
                font = ImageFont.truetype(path, font_size)
                break
            except (OSError, IOError):
                continue
        draw.text((ix, iy), text, fill=(0, 0, 0), font=font)
    out = np.array(img)
    rgb[:, :, :] = out[:, :, :]


def _load_graphics_data(path: str, ncols: int) -> Tuple[List[str], List[List[float]]]:
    """Load graphics data file: first line = headers (ncols tokens), rest = numeric rows. Returns (headers, rows)."""
    path_obj = Path(path)
    if not path_obj.exists():
        return [], []
    rows: List[List[float]] = []
    headers: List[str] = []
    with open(path_obj) as f:
        got_headers = False
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < ncols:
                continue
            if not got_headers:
                headers = parts[:ncols]
                got_headers = True
            else:
                try:
                    rows.append([float(parts[i]) for i in range(ncols)])
                except ValueError:
                    pass
    return headers, rows


def _draw_graphics_overlay(
    rgb: np.ndarray,
    entry: Tuple[int, int, str, int, int, int, int, int, int, int, int],
    frame_index: int,
) -> None:
    """Draw one graphics overlay: white rect, axes, data points and lines; highlight point at frame_index."""
    (_start, _end, path, ncols, ix, iy, width, height, point_size, col_x, col_y) = entry
    col_x = max(1, min(col_x, ncols))  # 1-based
    col_y = max(1, min(col_y, ncols))
    headers, rows = _load_graphics_data(path, ncols)
    if not rows or len(rows) == 0:
        return
    h, w = rgb.shape[0], rgb.shape[1]
    # Origin (ix, iy) = bottom-left of plot; plot area x in [ix, ix+width], y in [iy-height, iy]
    ix, iy = int(ix), int(iy)
    width, height = max(1, int(width)), max(1, int(height))
    point_size = max(0, int(point_size))
    # Clip to image
    x0, x1 = max(0, ix), min(w, ix + width)
    y0_plot_top = iy - height
    y1_plot_bottom = iy
    y0, y1 = max(0, y0_plot_top), min(h, y1_plot_bottom + 1)
    # Clear plot area to white
    for py in range(y0, y1):
        for px in range(x0, x1):
            if 0 <= px < w and 0 <= py < h:
                rgb[py, px, 0] = rgb[py, px, 1] = rgb[py, px, 2] = 255
    # Data columns (1-based)
    x_vals = [row[col_x - 1] for row in rows]
    y_vals = [row[col_y - 1] for row in rows]
    vmin1, vmax1 = min(x_vals), max(x_vals)
    vmin2, vmax2 = min(y_vals), max(y_vals)
    scale1 = (vmax1 - vmin1) / width if (vmax1 - vmin1) > 1e-10 else 1.0
    scale2 = (vmax2 - vmin2) / height if (vmax2 - vmin2) > 1e-10 else 1.0
    # Draw axes (bottom and left)
    for px in range(ix, min(w, ix + width + 1)):
        if 0 <= px < w and 0 <= iy < h:
            rgb[iy, px, 0] = rgb[iy, px, 1] = rgb[iy, px, 2] = 0
    for py in range(max(0, iy - height), min(h, iy + 1)):
        if 0 <= ix < w and 0 <= py < h:
            rgb[py, ix, 0] = rgb[py, ix, 1] = rgb[py, ix, 2] = 0
    # Points and lines (red); highlight frame_index in green with crosshairs
    nx_prev, ny_prev = None, None
    for i2, (xv, yv) in enumerate(zip(x_vals, y_vals)):
        nx = ix + int((xv - vmin1) / scale1)
        ny = iy - int((yv - vmin2) / scale2)  # y up in data -> py down in image
        nx = max(ix, min(ix + width, nx))
        ny = max(iy - height, min(iy, ny))
        is_highlight = i2 == frame_index
        psz = point_size + (1 if is_highlight else 0)
        r, g, b = (0, 255, 0) if is_highlight else (255, 0, 0)
        for dy in range(-psz, psz + 1):
            for dx in range(-psz, psz + 1):
                px_, py_ = nx + dx, ny + dy
                if 0 <= px_ < w and 0 <= py_ < h:
                    rgb[py_, px_, 0], rgb[py_, px_, 1], rgb[py_, px_, 2] = r, g, b
        if is_highlight:
            for px_ in range(ix, min(w, ix + width + 1)):
                if 0 <= ny < h:
                    rgb[ny, px_, 0] = rgb[ny, px_, 1] = rgb[ny, px_, 2] = 0
            for py_ in range(max(0, iy - height), min(h, iy + 1)):
                if 0 <= nx < w:
                    rgb[py_, nx, 0] = rgb[py_, nx, 1] = rgb[py_, nx, 2] = 0
        if nx_prev is not None and ny_prev is not None:
            nsteps = max(1, int(math.hypot(nx - nx_prev, ny - ny_prev)))
            for k in range(nsteps + 1):
                t = k / nsteps
                px_ = int(nx_prev + t * (nx - nx_prev))
                py_ = int(ny_prev + t * (ny - ny_prev))
                if 0 <= px_ < w and 0 <= py_ < h:
                    rgb[py_, px_, 0], rgb[py_, px_, 1], rgb[py_, px_, 2] = 255, 0, 0
        nx_prev, ny_prev = nx, ny
