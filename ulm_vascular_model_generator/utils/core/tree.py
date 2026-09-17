"""Core planar-or-volumetric DCCO vessel-tree object."""

from __future__ import annotations

import numpy as np

from ..config.yaml_config import DCCOConfig
from ..generation.growth import GrowthMixin
from ..geometry.branching import BranchingMixin
from ..geometry.constraints import ConstraintMixin
from ..hemodynamics.dimensional_models import haemodynamics_for_config
from ..hemodynamics.flow_scaling import FlowScalingMixin
from .models import Vessel
from .tree_index import SegmentIndex
from .tree_queries import TreeQueryMixin


class DCCOTree(
    GrowthMixin,
    ConstraintMixin,
    BranchingMixin,
    FlowScalingMixin,
    TreeQueryMixin,
):
    """Build a 2-D strip-lumen tree or a 3-D circular-tube tree."""

    def __init__(self, cfg: DCCOConfig):
        self.cfg = cfg
        self.haemodynamics = haemodynamics_for_config(cfg)
        self.rng = np.random.default_rng(cfg.seed)
        self.vessels: list[Vessel] = []
        self._fr_runtime = 1.0
        self._segment_index: SegmentIndex | None = None
        self._segment_index_size = -1
