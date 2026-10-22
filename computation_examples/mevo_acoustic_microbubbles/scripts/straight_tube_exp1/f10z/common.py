from pathlib import Path
import hashlib,json,subprocess,os
import numpy as np
from scripts.straight_tube_exp1.f10q.common import read,write,sha
from scripts.straight_tube_exp1.f10s.common import policy
ROOT=Path(__file__).resolve().parents[3]
BASE=ROOT/'results/straight_tube_exp1/f10z'
OLD=ROOT/'results/straight_tube_exp1/f10n'
CONFIG=ROOT/'configs/straight_tube_exp1/f10z/protocol.json'
LABELS=['EXPLORATORY','HOMOGENEOUS_WATER_CONTROL_PROXY','PHYSICS_UNCALIBRATED','ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED']
R=1.78e-6
PROTECTED=['src/mevo_mb/control/scheduled_pzt.py','src/mevo_mb/acoustics/water_disc_proxy.py','src/mevo_mb/fast_sim/water_sequential.py','src/mevo_mb/fast_sim/sequential_hard_sphere.py','src/mevo_mb/fast_sim/water_b2.py','src/mevo_mb/fast_sim/f10g_contact.py','src/mevo_mb/fast_sim/nearfield.py','src/mevo_mb/fast_sim/water_control_engine.py','src/mevo_mb/coupling/pdms_cavity.py','scripts/straight_tube_exp1/f10s/run.py','configs/straight_tube_exp1/f10q/protocol.json','configs/straight_tube_exp1/f10q/algorithm_revision_1.json']
INPUTS=['cache/coarse.h5','cache/fine.h5','inputs/response.h5','initial_states.npz','CONFIG_SNAPSHOT.json']
def identity():
 p=ROOT/'SOURCE_F10Z_COMMIT.txt'
 commit=p.read_text().strip() if p.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 paths=[ROOT/q for q in PROTECTED]+list((ROOT/'scripts/straight_tube_exp1/f10z').glob('*.py'))+[CONFIG]
 return dict(execution_commit=commit,source_sha256={str(q.relative_to(ROOT)):sha(q) for q in paths if q.exists()})
def verify_inputs():
 reg=read(BASE/'PREREGISTRATION.json')
 for rel,digest in reg['input_sha256'].items():assert sha(ROOT/rel)==digest,rel
 for rel,digest in reg['protected_source_sha256'].items():assert sha(ROOT/rel)==digest,rel
 return reg
