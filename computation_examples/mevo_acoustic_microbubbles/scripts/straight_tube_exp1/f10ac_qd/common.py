"""Isolated F10AC/AD records; historical mechanics remain read-only."""
from scripts.straight_tube_exp1.f10ab.common import ROOT,sha,read,write,clean,LABELS,R,policy
from pathlib import Path
import subprocess
BASE=ROOT/'results/straight_tube_exp1/f10ac_qd'
AD=ROOT/'results/straight_tube_exp1/f10ad_multifrequency'
AB=ROOT/'results/straight_tube_exp1/f10ab'
OLD=ROOT/'results/straight_tube_exp1/f10n'
CONFIG=ROOT/'configs/straight_tube_exp1/f10ac_qd/protocol.json'
def identity():
 m=ROOT/'SOURCE_F10AC_COMMIT.txt'
 return dict(source_commit=m.read_text().strip() if m.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),training_ready=False)
def verify():
 d=read(BASE/'PREREGISTRATION.json')
 for name,digest in d['protected_sha256'].items():
  if sha(ROOT/name)!=digest:raise ValueError('Frozen baseline changed: '+name)
 return d
