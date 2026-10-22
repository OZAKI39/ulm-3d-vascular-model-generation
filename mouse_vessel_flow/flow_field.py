"""Load frozen P1 fields in SI; canonical mesh ordering is checked before interpolation."""
from types import SimpleNamespace
import pyvista as pv
from .paths import PROJECT_ROOT, project_path, case_config
from .particles.audit import read_frozen, check_hash
from .particles.field import FrozenFEMField
from .particles.wall_geometry import WallGeometry
from .particles.validation_boundary import ValidationBoundaryClassifier
from .particles.hydrodynamic_resistance import viscosity_from_frozen
from .particles.inlet_flux import frozen_boundary_flux
from .particles.injection_admission import FiniteSizeAdmission

def load_case(name='roi_only', *, with_particles=True):
    """Read a case without solving CFD or changing saved results; coordinates m, U m/s, p Pa."""
    cfg=case_config(name); fem=project_path(cfg['fem_reference'])
    provenance, mesh, reference, boundaries=read_frozen(fem)
    path=project_path(cfg['flow']['path']);check_hash(path,cfg['flow']['sha256'])
    flow=pv.read(path)
    field=FrozenFEMField.from_grids(mesh,flow)
    audit,samplers=frozen_boundary_flux(mesh,flow,boundaries)
    mu,viscosity=viscosity_from_frozen(fem)
    if mu!=cfg['dynamic_viscosity_pa_s']:raise ValueError('Viscosity identity changed')
    env=SimpleNamespace(provenance=provenance,mesh=mesh,flow=flow,field=field,boundaries=boundaries,
        audit=audit,flux_audit=audit,sampler=samplers['INLET'],mu=mu,viscosity=viscosity,
        new_flow_path=path,new_flow_sha256=cfg['flow']['sha256'],repo_root=PROJECT_ROOT,case_config=cfg)
    if with_particles:
        from .particles.injection_method_c import TruncatedSonoVue
        env.wall=WallGeometry.from_frozen(fem)
        env.classifier=ValidationBoundaryClassifier({k:v for k,v in boundaries.items() if k.startswith('OUTLET_')})
        env.distribution=TruncatedSonoVue(project_path('inputs/reference/sonovue'))
        env.checker=FiniteSizeAdmission(None,None,wall=env.wall,field=field)
    return env
