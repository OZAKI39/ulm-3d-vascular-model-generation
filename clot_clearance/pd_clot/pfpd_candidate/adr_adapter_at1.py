"""Displacement-controlled tensile fixture and callback for the existing ADR.

The ADR boundary view has displaced grip coordinates. The original BAQP
family reference coordinates never change. There is no fluid force provider.
"""
from dataclasses import replace
import copy
import numpy as np
from ..geometry import make_cloud
from ..regularization.equilibrium import ForceState
from ..regularization import equilibrium_adr
from .kernel import Family
from .energy import Material
from .baqp import evaluate


def fixture(config):
    c = config['geometry']; h = c['spacing_m']
    cloud = make_cloud(dict(particle_spacing_m=h, cells=c['cells'], origin_m=c['origin_m'],
                            fixed_base_thickness_m=h, horizon_ratio=c['horizon_ratio']))
    bottom = cloud.X[:, 2] == cloud.X[:, 2].min()
    top = cloud.X[:, 2] == cloud.X[:, 2].max()
    cloud.fixed = bottom | top
    family = Family.from_cloud(cloud)
    m = config['material']
    material = Material(model=m['model'], mu=m['mu'], bulk=m['bulk'])
    if material.model != 'svk':
        raise ValueError('This benchmark requires SVK')
    return cloud, family, material, top


def pair_weights(family):
    i, j = family.pairs.T
    return 2*family.wbar*family.volume[i]*family.volume[j]


def boundary_view(cloud, top, strain):
    if not np.isfinite(strain) or strain < 0:
        raise ValueError('Nonnegative finite tensile strain required')
    X = cloud.X.copy()
    X[top, 2] += strain * np.ptp(cloud.X[:, 2])
    return replace(cloud, X=X)


class FrozenPhaseForce:
    def __init__(self, family, material, s_pf_bond):
        self.family, self.material = family, material
        s = np.array(s_pf_bond, float, copy=True)
        if s.shape != (len(family.pairs),) or not np.isfinite(s).all() or np.any((s < 0) | (s >= 1)):
            raise ValueError('Valid early phase required')
        self.g = (1-s)**2
        self.g.setflags(write=False)

    def __call__(self, x):
        f = self.family
        state = evaluate(f, x, self.material, g=self.g, h=np.ones(len(f.pairs)))
        i, j = f.pairs.T
        eps = np.maximum(np.linalg.norm(x[j]-x[i], axis=1)/f.length-1, 0.)
        return ForceState(state.force, np.zeros_like(x), eps, state.energy,
                          float(np.linalg.det(state.Ftilde).min()), np.zeros_like(x),
                          np.zeros(len(x), dtype=bool))


def equilibrate_tensile(cloud, family, material, top, x, s_pf_bond, strain, config):
    view = boundary_view(cloud, top, strain)
    initial = np.array(x, copy=True)
    initial[view.fixed] = view.X[view.fixed]
    # Legacy guard interface only: the radius is a loose numerical coordinate
    # envelope, NOT a modeled vessel or fluid. Density is diagnostic metadata.
    adapter_config = dict(mechanical_equilibrium=copy.deepcopy(config['mechanical_equilibrium']),
                          clot=dict(density_kg_m3=1000.), pipe=dict(radius_m=1.),
                          safety=dict(maximum_displacement_m=.01, maximum_radius_ratio=1.))
    return equilibrium_adr.equilibrate(view, initial, np.zeros_like(initial),
        np.ones(len(family.pairs)), adapter_config,
        force_callback=FrozenPhaseForce(family, material, s_pf_bond),
        method='adaptive_dynamic_relaxation')


def tensile_domain(state, tolerance):
    green = .5*(np.swapaxes(state.Ftilde, 1, 2) @ state.Ftilde-np.eye(3))
    minimum = float(np.linalg.eigvalsh(green).min())
    if minimum < tolerance:
        raise ValueError(f'PURE_TENSION_DOMAIN_VIOLATION: minimum Green strain {minimum}')
    return minimum
