"""Validation-only common-geometry study; no solver or particle jobs."""
from pathlib import Path
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parents[2]
PREVIOUS = HERE.parent / 'common_curved_geometry_review_20261006'
NATIVE = HERE.parent / 'native_p2p1_review_20261006'
IDENTITY = HERE.parent / 'stabilization_local_conservation_review_20261006'
RETAINED = BASE / 'gpu_solver_fix/source'
SOURCE = RETAINED / 'Code/Source/solver'

def read(p):
    return json.loads(Path(p).read_text())

def put(name, value):
    (HERE/name).write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def native_env():
    import os
    env = dict(os.environ)
    old = read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json')
    env['LD_LIBRARY_PATH'] = ':'.join([str(Path(old['PETSc']['path']).parent),
        '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib',env.get('LD_LIBRARY_PATH','')])
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    return env
