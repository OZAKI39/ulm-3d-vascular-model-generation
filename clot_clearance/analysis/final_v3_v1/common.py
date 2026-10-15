from pathlib import Path
import json,hashlib
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path(C['server']);O=R/'verification'/C['task']/C['run'];F=R/'results'/C['task']/C['run'];P=O/'C1_final';OLD=R/'verification/y_clot_point_cloud_and_early_flow_coupling_v1/cases/C1_GAUSSIAN_A200';PRIOR=R/'verification/lps_native_state_identity_and_fracture_bugfix_v1/20261011T203015Z'
S=PRIOR/'source_owner_fix';B=PRIOR/'build_owner_fix';PROD=R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so'
EXPECTED='00b246ebe4224e648c27b74afb64fc22b9cea4ba6ffb23f7c40d96c37ca22ae9'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False))
def loaded(L):
 import sys,lammps,mpi4py
 paths=sorted({s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps.so' in s or 'libmpi.so' in s})
 libs={p:sha(p) for p in paths};ll=[p for p in paths if 'liblammps' in p]
 assert len(ll)==1 and libs[ll[0]]==EXPECTED and Path(ll[0]).parent==B,(ll,libs)
 return dict(python=sys.executable,python_version=sys.version,lammps_python=lammps.__file__,mpi4py=mpi4py.__file__,lammps_version=L.version(),loaded_libraries=libs)
