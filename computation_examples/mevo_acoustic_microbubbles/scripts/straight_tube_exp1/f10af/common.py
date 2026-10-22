"""F10AF immutable identities; old physics and histories remain separate."""
from pathlib import Path
import datetime, subprocess
from scripts.straight_tube_exp1.f10ab.common import sha, read, write, clean, R, policy
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/straight_tube_exp1/f10af'
AE=ROOT/'results/straight_tube_exp1/f10ae'
OLD=ROOT/'results/straight_tube_exp1/f10n'
AD=ROOT/'results/straight_tube_exp1/f10ad_multifrequency'
CONFIG=ROOT/'configs/straight_tube_exp1/f10af/protocol.json'
LABELS=['EXPLORATORY','HOMOGENEOUS_WATER_CONTROL_PROXY','PHYSICS_UNCALIBRATED','ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED','REAL_CAMERA_OBSERVABILITY_UNVERIFIED']
def identity():
    p=ROOT/'SOURCE_F10AF_COMMIT.txt'
    return dict(source_commit=p.read_text().strip() if p.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),labels=LABELS,training_ready=False)
def verify():
    for p,h in read(BASE/'PREREGISTRATION.json')['protected_sha256'].items():
        if sha(ROOT/p)!=h:raise RuntimeError('Protected baseline changed: '+p)
