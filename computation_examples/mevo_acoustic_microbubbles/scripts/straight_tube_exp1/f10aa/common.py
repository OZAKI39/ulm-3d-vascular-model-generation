from pathlib import Path
import subprocess, numpy as np
from scripts.straight_tube_exp1.f10s.common import policy, context
from scripts.straight_tube_exp1.f10q.common import read, write, sha
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/straight_tube_exp1/f10aa'
CONFIG=ROOT/'configs/straight_tube_exp1/f10aa/protocol.json'
SCALING=ROOT/'results/straight_tube_exp1/f10z_scaling'
OLD=ROOT/'results/straight_tube_exp1/f10n'
LABELS=['EXPLORATORY','PHYSICS_UNCALIBRATED','HOMOGENEOUS_WATER_CONTROL_PROXY','ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED']
R=1.78e-6
def identity():
 p=ROOT/'SOURCE_F10AA_COMMIT.txt'
 commit=p.read_text().strip() if p.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 paths=list((ROOT/'scripts/straight_tube_exp1/f10aa').glob('*.py'))+list((ROOT/'src/mevo_mb/fast_sim').glob('accelerated*.py'))+[CONFIG]
 return dict(commit=commit,sha256={str(p.relative_to(ROOT)):sha(p) for p in paths})
def verify():
 d=read(BASE/'PREREGISTRATION.json')
 for k in ['input_sha256','frozen_source_sha256']:
  for p,h in d[k].items(): assert sha(ROOT/p)==h,p
 return d
def stats(values):
 a=np.asarray(values)
 return dict(mean=float(a.mean()),p50=float(np.median(a)),p95=float(np.quantile(a,.95)),max=float(a.max()),samples=len(a))
