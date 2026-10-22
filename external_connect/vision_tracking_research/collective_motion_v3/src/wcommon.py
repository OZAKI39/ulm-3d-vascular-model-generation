"""F10W paths/identity. Existing modules are read-only dependencies."""
import os,sys,json,hashlib,subprocess,datetime
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('F10W_OUTPUT_ROOT',ROOT))
VIDEO=Path(os.environ.get('F10W_VIDEO','/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4'))
THIRD=Path(os.environ.get('F10W_THIRD_PARTY',ROOT/'third_party'))
MODELS=Path(os.environ.get('F10W_MODELS',ROOT/'models'))
P=json.loads((ROOT/'configs/protocol.json').read_text())
sys.path.insert(0,str(ROOT.parent/'collective_motion_v2/src'))
import legacy as old
from common import NAMES,ResourceMonitor
from methods import preprocess
from state_v2 import new_state,to_v1,VersionedGate
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(tuple,list,np.ndarray)):return [clean(v) for v in x]
    if isinstance(x,np.bool_):return bool(x)
    if isinstance(x,np.integer):return int(x)
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    if isinstance(x,Path):return str(x)
    return x
def dump(path,x):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(clean(x),indent=2,ensure_ascii=False,allow_nan=False)+'\n');tmp.replace(path)
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def identity():
    try:
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,stderr=subprocess.DEVNULL,text=True).strip()
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    except subprocess.CalledProcessError:
        v=json.loads((ROOT/'EXECUTION.json').read_text());head=v['git_head'];branch=v['git_branch']
    return dict(git_head=head,git_branch=branch,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_sha256={str(p.relative_to(ROOT)):sha(p) for folder in ['src','scripts','configs'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts})
def stats(x):
    v=np.asarray(x,float);v=v[np.isfinite(v)]
    return dict(n=len(v),p50=float(np.median(v)),p95=float(np.percentile(v,95)),max=float(v.max()),mean=float(v.mean())) if len(v) else dict(n=0,p50=None,p95=None,max=None,mean=None)
