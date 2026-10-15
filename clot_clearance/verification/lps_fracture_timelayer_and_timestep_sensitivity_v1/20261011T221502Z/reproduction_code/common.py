from pathlib import Path
import json,hashlib
C=json.loads(Path(__file__).with_name('config.json').read_text())
R=Path(C['server']);O=R/'verification'/C['task']/C['run'];F=R/'results'/C['task']/C['run']
PRIOR=R/C['prior_bugfix'];V3=R/C['prior_v3'];S=PRIOR/'source_owner_fix';B=PRIOR/'build_owner_fix'
EXPECTED=C['library_sha256'];PROD=R/'tools/lammps-native-y-wall-shared-v1-build/liblammps.so'
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False))
