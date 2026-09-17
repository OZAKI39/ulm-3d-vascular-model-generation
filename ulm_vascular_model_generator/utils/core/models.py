"""
The basic data structure of the vascular tree.
This file is intentionally small: all "vessel segments" in the algorithm are represented by `Vessel`.
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from ..geometry.geometry import _distance


@dataclass
class Vessel:
    """
    Define a single vessel segment. 
    In the current model, a `Vessel` is not a whole tree or a curved line, 
    but a straight cylinder from `x_p` to `x_d`. 
    Complex vascular networks are formed by connecting many such small cylinders end-to-end.
    """

    # segment ID
    vid: int
    # parent vessel ID; root vessels have no parent, represented by -1.
    parent_id: int
    # list of child vessel IDs; used in topological traversal and flow accumulation.
    children: list[int] = field(default_factory=list)

    # proximal point, representing the blood flow inlet side coordinate.
    x_p: np.ndarray = field(default_factory=lambda: np.zeros(3))
    # distal point, representing the blood flow outlet side coordinate.
    x_d: np.ndarray = field(default_factory=lambda: np.zeros(3))

    # Lumen half-width in planar_2d or circular radius in volumetric_3d [um].
    # The run metadata and format-v3 field name make the interpretation explicit.
    radius: float = 1.0

    # how to branch this vessel
    branching_mode: str = "versatile"
    # the role of this vessel in the vascular tree
    role: str = "distribution"

    is_main_trunk: bool = False

    # specified outlet flow; if not specified, it will be handled by terminal_outflow.
    prescribed_outflow: float = 0.0
    # Planar flux [um^2/s] in planar_2d or volume flow [um^3/s] in volumetric_3d.
    flow_rate: float = 0.0
    # Mean strip velocity q/(2R) in 2-D or circular-tube Q/(pi*R^2) in 3-D.
    mean_velocity: float = 0.0
    # Relative local conservation and Murray residuals used for validation.
    flow_conservation_residual: float = 0.0
    murray_residual: float = 0.0
    def length(self) -> float:
        """
        Return the length of the vessel segment, 
        which is used in geometric constraints and volume costs."""
        return _distance(self.x_d, self.x_p)
