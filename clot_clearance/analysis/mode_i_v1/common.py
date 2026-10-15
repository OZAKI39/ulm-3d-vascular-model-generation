from pathlib import Path
import json,hashlib
import numpy as np
C=json.loads(Path(__file__).with_name('config.json').read_text());R=Path(C['server']);O=R/'verification'/C['task']/C['run'];F=R/'results'/C['task']/C['run'];W=R/C['wide'];PRIOR=R/C['prior_bugfix'];S=PRIOR/'source_owner_fix';B=PRIOR/'build_owner_fix';PROD=R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so';EXPECTED=C['library_sha256'];DX=50e-6;DELTA=150.5e-6;K=100000/3;G=940.;S00=.15;V=DX**3;PERIOD=.0004

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as fp:
  for v in iter(lambda:fp.read(4194304),b''):h.update(v)
 return h.hexdigest()
def dump(p,x):
 def cv(z):
  if isinstance(z,np.ndarray):return z.tolist()
  if isinstance(z,np.generic):return z.item()
  if isinstance(z,Path):return str(z)
  raise TypeError(type(z).__name__)
 Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,default=cv,allow_nan=False))

SURFACE=R/C['prior_surface']
