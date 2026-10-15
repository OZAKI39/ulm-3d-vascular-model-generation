from pathlib import Path
import json,hashlib
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path(C['server']);O=R/'verification'/C['task']/C['run'];F=R/'results'/C['task']/C['run'];MAP=R/C['mapping'];CFD=R/C['cfd'];PRIOR=R/C['prior_bugfix'];V3=R/C['prior_v3'];S=PRIOR/'source_owner_fix';B=PRIOR/'build_owner_fix';PROD=R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so';EXPECTED=C['library_sha256'];GDIR=R/'verification/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1';DX=50e-6;DELTA=150.5e-6;V=DX**3;K=100000/3;G=940.;N=3488;NC=1927

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for a in iter(lambda:f.read(4194304),b''):h.update(a)
 return h.hexdigest()
def dump(p,obj):Path(p).write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False))
def loaded(L):
 import sys,lammps,mpi4py
 paths=sorted({s.split()[-1] for s in Path('/proc/self/maps').read_text().splitlines() if 'liblammps.so' in s or 'libmpi.so' in s});libs={p:sha(p) for p in paths};ll=[p for p in paths if 'liblammps' in p];assert len(ll)==1 and libs[ll[0]]==EXPECTED and Path(ll[0]).parent==B
 return dict(python=sys.executable,python_version=sys.version,lammps_python=lammps.__file__,mpi4py=mpi4py.__file__,lammps_version=L.version(),actual_libraries=libs)
