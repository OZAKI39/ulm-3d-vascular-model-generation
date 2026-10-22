from pathlib import Path
import subprocess
from scripts.straight_tube_exp1.f10z.common import ROOT, OLD, R, LABELS, INPUTS, PROTECTED
from scripts.straight_tube_exp1.f10s.common import policy
from scripts.straight_tube_exp1.f10q.common import read, write, sha

BASE = ROOT / 'results/straight_tube_exp1/f10z_scaling'
CONFIG = ROOT / 'configs/straight_tube_exp1/f10z_scaling/protocol.json'
HISTORICAL = ROOT / 'results/straight_tube_exp1/f10z'

def identity():
    marker = ROOT / 'SOURCE_F10Z_SCALING_COMMIT.txt'
    commit = marker.read_text().strip() if marker.exists() else subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    paths = [ROOT / q for q in PROTECTED] + list(Path(__file__).parent.glob('*.py')) + [CONFIG]
    return dict(execution_commit=commit, source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths})

def verify():
    reg = read(BASE / 'PREREGISTRATION.json')
    for collection in ['input_sha256', 'protected_source_sha256']:
        for rel, digest in reg[collection].items():
            assert sha(ROOT / rel) == digest, rel
    return reg
