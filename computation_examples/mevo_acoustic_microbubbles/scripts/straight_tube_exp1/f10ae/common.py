"""F10AE isolated paths, immutable source identity and strict input hashes."""
import subprocess, datetime
from scripts.straight_tube_exp1.f10ab.common import ROOT, sha, read, write, clean, LABELS, R, policy
BASE = ROOT/'results/straight_tube_exp1/f10ae'
AC = ROOT/'results/straight_tube_exp1/f10ac_qd'
AD = ROOT/'results/straight_tube_exp1/f10ad_multifrequency'
AB = ROOT/'results/straight_tube_exp1/f10ab'
OLD = ROOT/'results/straight_tube_exp1/f10n'
CONFIG = ROOT/'configs/straight_tube_exp1/f10ae/protocol.json'
def identity():
    p=ROOT/'SOURCE_F10AE_COMMIT.txt'
    commit=p.read_text().strip() if p.exists() else subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    return dict(source_commit=commit,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),labels=LABELS,training_ready=False)
def verify():
    for p,h in read(BASE/'PREREGISTRATION.json')['protected_sha256'].items():
        if sha(ROOT/p)!=h: raise RuntimeError('Protected input changed: '+p)
