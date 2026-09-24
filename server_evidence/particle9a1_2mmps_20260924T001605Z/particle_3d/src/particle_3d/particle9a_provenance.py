"""Fail-closed current production input and source identity; no historical cache."""
from pathlib import Path
import hashlib,json

FLOW_SHA='129ebb77396550b300b20ce29cde7b03919e6037b7a81de60fe75d6c6efb616d'
OLD_FLOW_SHA='373ff549430cab710d57eb4406f2e7588f7a5bb68ebaa68cf32f8da5662e26f1'
ROLE='MEAN_2P0_MMPS'
REPO=Path(__file__).resolve().parents[3]
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular'
ADMISSION_VERSION='P82A_METHOD_B_FIXED_FLUX_ANCHOR_SIZE_RETRY_512'
MODEL='P9A_SINGLE_NEAREST_PLANAR_WALL_SPHERE_WITH_UNCHANGED_P65_NORMAL'


def sha256(path):
    with open(path,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def require_current_flow(fem=FEM):
    root=Path(fem);frozen=root/'frozen_reference'
    data=json.loads((frozen/'flow/flow_field_manifest.json').read_text())
    if data.get('case_role')!=ROLE or data.get('sha256')!=FLOW_SHA or 'steady_flow_stage_sv1_3q' in data.get('path',''):
        raise ValueError('CURRENT_PRODUCTION_REJECTS_LEGACY_OR_UNBOUND_FLOW')
    if any(frozen.rglob('steady_flow_stage_sv1_3q.vtu')):
        raise ValueError('LEGACY_ACTIVE_FLOW_COPY_PRESENT')
    if sha256(root/data['path'])!=FLOW_SHA:raise ValueError('CURRENT_FLOW_SHA_MISMATCH')
    if abs(data['inlet_mean_velocity_m_s']-.002)>2e-9:raise ValueError('CURRENT_FLOW_MEAN_MISMATCH')
    return dict(case_role=ROLE,flow_sha256=FLOW_SHA,flow_path=data['path'],
        inlet_mean_velocity_m_s=data['inlet_mean_velocity_m_s'])


def source_identity():
    names=['planar_wall_hydrodynamics.py','particle9a_motion.py','particle9a_provenance.py',
           'particle9a_production.py','particle82a_integration.py','particle82a_admission.py',
           'particle81_cache.py','particle81_simulation.py','particle65_motion.py','field.py',
           'nearfield_handoff.py','nearfield_regularization.py','wall_gap.py','wall_geometry.py',
           'resistance_solver.py','inlet_flux.py','injection_admission.py','injection_population.py']
    hashes={name:sha256(Path(__file__).parent/name) for name in names}
    return dict(source_sha256=hashes,p9a_source_sha256=hashlib.sha256(json.dumps(hashes,sort_keys=True).encode()).hexdigest())
