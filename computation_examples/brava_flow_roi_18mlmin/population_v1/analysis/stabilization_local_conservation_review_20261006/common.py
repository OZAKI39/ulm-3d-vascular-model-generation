"""Independent read-only field/source review; no production solver mutation."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
PREVIOUS = HERE.parent/'field_reliability_event_review'
CAUSE = HERE.parent/'speed_event_cause_review'
Q = Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q')
SOURCE = Q/'external/svMultiPhysics-reuse'
NAMES = ['M0', 'M1_fixedplc_medium', 'M1_fixedplc_fine']
PIDS = [1189, 1329, 642, 822]
sys.path.insert(0, str(HERE/'runtime/rt_dependencies'))

def read(path):
    return json.loads(Path(path).read_text())

def put(name, value):
    (HERE/name).write_text(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def paths(name):
    if name == 'M0':
        return BASE/'mesh/SV_MESH/mesh-complete.mesh.vtu', BASE/'cases/balanced_pressure_final/frozen_flow/steady_flow.vtu'
    return (BASE/'mesh_particle_validation_v2/02_mesh_candidates'/name/'mesh/SV_MESH/mesh-complete.mesh.vtu',
            BASE/'mesh_particle_validation_v2/03_cfd_runs'/name/'frozen_flow/steady_flow.vtu')
