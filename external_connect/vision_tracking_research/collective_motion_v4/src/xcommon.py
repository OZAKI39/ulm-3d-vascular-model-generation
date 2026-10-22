"""F10X paths, immutable execution identity, and read-only F10U math reuse."""
import os,sys,json,time,hashlib,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=Path(os.environ.get('F10X_OUTPUT',ROOT))
SERVER=Path(os.environ.get('F10X_SERVER','/workspace/f10x_piv_20261017'))
THIRD=SERVER/'third_party' if SERVER.exists() else ROOT/'third_party'
VIDEO=Path(os.environ.get('F10X_VIDEO','/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4'))
VIDEO_SHA='3fc13aaa611a4eeb94b27804524318b3a8411cc38967804c6620913cfb9c484c'
sys.path.insert(0,str(ROOT.parent/'collective_motion_v1/src'))
from common import clean,dump,sha,stats,NAMES,ResourceMonitor
from affine import fit_affine,parameters,empty
for d in ['results','figures','videos','logs','evidence','data']:(OUT/d).mkdir(parents=True,exist_ok=True)
def identity():
 p=ROOT/'EXECUTION.json'
 meta=json.loads(p.read_text()) if p.exists() else {'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()}
 return dict(**meta,source_sha256={str(p.relative_to(ROOT)):sha(p) for d in ['src','scripts','configs'] for p in (ROOT/d).rglob('*') if p.is_file() and '__pycache__' not in str(p)},utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
def frames(indices=None):
 import av
 with av.open(str(VIDEO)) as c:
  for i,f in enumerate(c.decode(video=0)):
   if indices is None or i in indices:yield i,float(f.pts*f.time_base),f.to_ndarray(format='gray')
   if indices is not None and i>=max(indices):break
