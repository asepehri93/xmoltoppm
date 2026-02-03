"""PPM (and optional PNG) writer for RGB arrays."""

from __future__ import annotations

from pathlib import Path
from typing import Union

import numpy as np


def write_ppm(path: Union[str, Path], rgb: np.ndarray) -> None:
    """Write RGB array (H, W, 3) uint8 to PPM P3 format. maxval 255."""
    path = Path(path)
    if rgb.ndim != 3 or rgb.shape[2] != 3:
        raise ValueError("rgb must be (H, W, 3)")
    h, w = rgb.shape[0], rgb.shape[1]
    with open(path, "w") as f:
        f.write("P3\n")
        f.write(f"{w} {h}\n")
        f.write("255\n")
        for y in range(h):
            row = []
            for x in range(w):
                r, g, b = int(rgb[y, x, 0]), int(rgb[y, x, 1]), int(rgb[y, x, 2])
                row.extend([str(r), str(g), str(b)])
            f.write(" ".join(row) + "\n")


def write_png(path: Union[str, Path], rgb: np.ndarray) -> None:
    """Write RGB array (H, W, 3) uint8 to PNG. Requires Pillow."""
    try:
        from PIL import Image
    except ImportError:
        raise ImportError("PNG output requires Pillow: pip install Pillow") from None
    path = Path(path)
    if rgb.ndim != 3 or rgb.shape[2] != 3:
        raise ValueError("rgb must be (H, W, 3)")
    img = Image.fromarray(rgb, mode="RGB")
    img.save(path)
