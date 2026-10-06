"""Read-only event field diagnostics; no CFD solve or particle integration."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
OLD = HERE.parent / 'speed_event_cause_review'
NAMES = ['M0', 'M1_fixedplc_medium', 'M1_fixedplc_fine']
PIDS = [1189, 1329, 642, 822]
SOURCE = BASE / 'population_v1/source/particle_3d/src'
if not SOURCE.exists():
    SOURCE = BASE.parent / '_shared/particle_workflow/particle_3d/src'
sys.path.insert(0, str(SOURCE))
sys.path.insert(0, str(OLD))

def read(p):
    return json.loads(Path(p).read_text())

def put(name, value):
    (HERE / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def paths(name):
    if name == 'M0':
        return BASE/'mesh', BASE/'cases/balanced_pressure_final'
    return (BASE/'mesh_particle_validation_v2/02_mesh_candidates'/name/'mesh',
            BASE/'mesh_particle_validation_v2/03_cfd_runs'/name)

def meshpath(name):
    return paths(name)[0]/'SV_MESH/mesh-complete.mesh.vtu'

def flowpath(name):
    return paths(name)[1]/'frozen_flow/steady_flow.vtu'

def unit(x):
    x=np.asarray(x, float)
    return x/np.linalg.norm(x)

def norm(x):
    return np.linalg.norm(x, axis=-1)

def finite_json(value):
    if isinstance(value, np.ndarray): return finite_json(value.tolist())
    if isinstance(value, dict): return {str(k):finite_json(v) for k,v in value.items()}
    if isinstance(value, (tuple,list)): return [finite_json(v) for v in value]
    if isinstance(value, (float,np.floating)): return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer): return int(value)
    if isinstance(value, np.bool_): return bool(value)
    return value
