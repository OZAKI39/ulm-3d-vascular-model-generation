from pathlib import Path
import hashlib,json,subprocess,numpy as np
from scripts.straight_tube_exp1.f10s.common import context,policy
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/straight_tube_exp1/f10ab'
CONFIG=ROOT/'configs/straight_tube_exp1/f10ab/protocol.json'
OLD=ROOT/'results/straight_tube_exp1/f10n'
LABELS=['EXPLORATORY','HOMOGENEOUS_WATER_CONTROL_PROXY','PHYSICS_UNCALIBRATED','ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED']
R=1.78e-6

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def read(p):return json.loads(Path(p).read_text())
def clean(x):
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,complex):return {'real':x.real,'imag':x.imag}
 if isinstance(x,float) and not np.isfinite(x):return None
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [clean(v) for v in x]
 return x
def write(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(clean(x),ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def identity():
 marker=ROOT/'SOURCE_F10AB_COMMIT.txt'
 commit=marker.read_text().strip() if marker.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 return dict(source_commit=commit,source_sha256={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'scripts/straight_tube_exp1/f10ab').glob('*.py')})
def verify():
 d=read(BASE/'PREREGISTRATION.json')
 for group in ['model_sha256','input_sha256']:
  for name,digest in d[group].items():
   if sha(ROOT/name)!=digest:raise ValueError('Frozen input changed: '+name)
 return d
