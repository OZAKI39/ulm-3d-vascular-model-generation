"""LOCALITY-STRONG in-domain constrained-support reference; all nodes physical.

This is explicitly not an exterior ghost-material or energy completion method.
The constraint width is one direct interaction horizon, not a 3-delta rule.
"""
import numpy as np


def constrained_support(geometry, family, surrogate):
    mask = geometry.masks(family.X, family.horizon)['direct_horizon_region']
    return mask, surrogate.position(family.X)


def impose(x, mask, prescribed):
    result = np.array(x, copy=True); result[mask] = prescribed[mask]; return result
