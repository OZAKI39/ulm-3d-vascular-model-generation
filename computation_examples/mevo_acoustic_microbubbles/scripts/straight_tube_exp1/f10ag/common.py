"""Isolated F10AG orchestration: original physical functions stay unchanged."""
from pathlib import Path
import datetime,subprocess
from scripts.straight_tube_exp1.f10ab.common import sha,read,write,clean,R,policy
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/straight_tube_exp1/f10ag'
AF=ROOT/'results/straight_tube_exp1/f10af'
CONFIG=ROOT/'configs/straight_tube_exp1/f10ag/protocol.json'
LABELS=['EXPLORATORY','HOMOGENEOUS_WATER_CONTROL_PROXY','PHYSICS_UNCALIBRATED','ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED','STREAMING_NOT_MODELED','ROTATIONAL_ACOUSTIC_FORCE_DIAGNOSTIC','COLLECTIVE_ROTATION_UNCONFIRMED']
def identity():
 p=ROOT/'SOURCE_F10AG_COMMIT.txt'
 return dict(source_commit=p.read_text().strip() if p.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),labels=LABELS,training_ready=False)
def verify():
 for p,h in read(BASE/'PREREGISTRATION.json')['protected_sha256'].items():
  if sha(ROOT/p)!=h:raise RuntimeError('Protected input changed: '+p)
def bind():
 # These module globals only select inputs/output/identity; no physical code replacement.
 from scripts.straight_tube_exp1.f10af import run,slow,audit,paired
 from scripts.straight_tube_exp1.f10aa import trajectory
 for m in [run,slow,audit,paired]:
  m.BASE=BASE;m.CONFIG=CONFIG;m.identity=identity;m.verify=verify
 trajectory.BASE=BASE
 return run,slow,audit,paired
