"""Frozen surrogate-to-load conversion. The identity is not an acceptance gate."""
from dataclasses import dataclass
from types import SimpleNamespace
import numpy as np
from ..baqp import evaluate
from ..boundary_force_diagnostics import force_terms
from .geometry import influence_support


@dataclass(frozen=True)
class ConvertedLoad:
    nodal_force_N: np.ndarray
    force_density_N_m3: np.ndarray
    support: np.ndarray
    seed: np.ndarray
    omitted_internal_force_N: np.ndarray


def convert_volume_force(geometry, family, surrogate, state=None):
    state = evaluate(family, surrogate.position(family.X), surrogate.material) if state is None else state
    support, seed = influence_support(geometry, family)
    load = np.where(support[:, None], -state.force, 0.)
    load.setflags(write=False)
    return ConvertedLoad(load, load/family.volume[:, None], support, seed,
                         np.where(support[:, None], 0., state.force))


def stress_channel_load(family, state, P):
    """Linear-in-stress diagnostic at frozen geometry; not a new constitutive law."""
    proxy = SimpleNamespace(P=np.broadcast_to(P, state.P.shape), shape_i=state.shape_i, shape_j=state.shape_j)
    return -force_terms(family, proxy)['total']
