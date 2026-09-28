"""An explicit environment adapter; no historical flow identity or bypass."""
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pyvista as pv
from ..field import FrozenFEMField
from ..wall_geometry import WallGeometry
from ..validation_boundary import ValidationBoundaryClassifier
from ..inlet_flux import frozen_boundary_flux
from ..injection_admission import FiniteSizeAdmission
from ..injection_method_c import TruncatedSonoVue
from ..open_boundary_rim import build_rim_topology
from .provenance import digest


def load_environment(spec):
    """All paths/hashes are supplied by the caller's authoritative case manifest."""
    required = ['mesh', 'flow', 'boundaries', 'sha256', 'mu_Pa_s', 'sonovue_root']
    if any(k not in spec for k in required):
        raise ValueError('Incomplete authoritative flow environment')
    paths = dict(mesh=spec['mesh'], flow=spec['flow'], **spec['boundaries'])
    if set(paths) != set(spec['sha256']):
        raise ValueError('Every input requires an expected SHA256')
    for key, path in paths.items():
        if digest(path) != spec['sha256'][key]:
            raise ValueError('Actual flow/mesh/boundary provenance mismatch: ' + key)
    mu = float(spec['mu_Pa_s'])
    if not np.isfinite(mu) or mu <= 0:
        raise ValueError('Positive dynamic viscosity in Pa.s required')
    mesh, flow = pv.read(paths['mesh']), pv.read(paths['flow'])
    field = FrozenFEMField.from_grids(mesh, flow)
    boundaries = {k: pv.read(v) for k, v in spec['boundaries'].items()}
    w = boundaries['WALL']; faces = w.faces.reshape(-1, 4)[:, 1:]
    wall = WallGeometry(np.asarray(w.points)[faces], provenance=dict(
        role='REAL_AUTHORITATIVE_FEM_WALL', sha256=spec['sha256']['WALL'],
        physical_coating_thickness_m=0., wall_rigid=True, wall_stationary=True))
    wall.global_node_ids = np.asarray(w.point_data['GlobalNodeID'])[faces] - 1
    caps = {k: np.asarray(v.point_data['GlobalNodeID'])[v.faces.reshape(-1,4)[:,1:]]-1
            for k, v in boundaries.items() if k != 'WALL'}
    wall.open_boundary_topology = build_rim_topology(wall.global_node_ids, caps, spec['sha256'])
    classifier = ValidationBoundaryClassifier({k: v for k, v in boundaries.items() if k.startswith('OUTLET_')})
    audit, samplers = frozen_boundary_flux(mesh, flow, boundaries)
    env = SimpleNamespace(field=field, mesh=mesh, wall=wall, mu=mu,
                          classifier=classifier, boundaries=boundaries,
                          sampler=samplers['INLET'], flux_audit=audit,
                          distribution=TruncatedSonoVue(Path(spec['sonovue_root'])),
                          flow_sha256=spec['sha256']['flow'], provenance=spec)
    env.checker = FiniteSizeAdmission(None, None, wall=wall, field=field)
    return env
