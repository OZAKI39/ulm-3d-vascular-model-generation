from pathlib import Path
import os,sys,json,hashlib,time,subprocess
import numpy as np
from collections import deque
ROOT=Path(__file__).resolve().parents[1]
# Read-only F10X computation imports; all generated files redirected to F10Y.
os.environ.setdefault('F10X_OUTPUT',str(ROOT/'runtime/f10x_import'))
sys.path.append(str(ROOT.parent/'collective_motion_v4/src'))
def clean(x):
 if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
 if isinstance(x,(list,tuple,deque)):return [clean(v) for v in x]
 if isinstance(x,np.ndarray):return clean(x.tolist())
 if isinstance(x,np.generic):return clean(x.item())
 if isinstance(x,float) and not np.isfinite(x):return None
 return x

def dump(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(clean(x),indent=2,allow_nan=False))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def stats(a):
 a=np.asarray(a,float);a=a[np.isfinite(a)]
 return dict(n=len(a),p50=float(np.median(a)),p95=float(np.percentile(a,95)),max=float(a.max()),mean=float(a.mean())) if len(a) else dict(n=0,p50=None,p95=None,max=None,mean=None)
def identity():
 p=ROOT/'EXECUTION.json'
 return json.loads(p.read_text()) if p.exists() else {'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'code_sha':{str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'src').glob('*.py')}}
