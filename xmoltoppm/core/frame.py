"""Frame: one XYZ frame (coordinates, cell, metadata)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class Frame:
    """One trajectory frame: atoms, cell, optional velocities/strain."""

    nat: int
    coords: np.ndarray  # (N, 3)
    symbols: list[str]
    cell: Optional[np.ndarray] = None  # (3, 3) or None
    molname: str = ""
    iteration: int = 0
    energy: float = 0.0
    velocities: Optional[np.ndarray] = None  # (N, 3)
    estrain: Optional[np.ndarray] = None  # (N,)

    def __post_init__(self) -> None:
        if self.coords.shape[0] != self.nat or self.coords.shape[1] != 3:
            raise ValueError("coords must be (nat, 3)")
        if len(self.symbols) != self.nat:
            raise ValueError("symbols length must equal nat")
        if self.velocities is not None and self.velocities.shape != (self.nat, 3):
            raise ValueError("velocities must be (nat, 3)")
        if self.estrain is not None and self.estrain.shape != (self.nat,):
            raise ValueError("estrain must be (nat,)")
