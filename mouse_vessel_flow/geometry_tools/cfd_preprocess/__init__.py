"""Global-to-ROI boundary-condition preprocessing."""

from .config import CFDPreprocessConfig, load_cfd_preprocess_config
from .one_d_flow import GlobalFlowResult, edge_resistance, solve_global_flow

__all__ = [
    "CFDPreprocessConfig",
    "GlobalFlowResult",
    "edge_resistance",
    "load_cfd_preprocess_config",
    "solve_global_flow",
]
