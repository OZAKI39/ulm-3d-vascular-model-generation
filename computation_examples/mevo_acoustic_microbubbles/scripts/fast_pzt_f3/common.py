from pathlib import Path
import json
import hashlib
import numpy as np
from scripts.fast_pzt_f2.common import write as write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f3"
F2 = ROOT / "results/fast_pzt/f2"
PROTOCOL = json.loads((ROOT / "configs/fast_pzt/f3/protocol.json").read_text())


def js(path):
    return json.loads(Path(path).read_text())


def state():
    return np.load(F2 / "mature_v1/MATURE_STATE.npz")["positions_m"]


def identity():
    files = sorted((ROOT / "src/mevo_mb/fast_sim").glob("*.py")) + sorted(
        (ROOT / "scripts/fast_pzt_f3").glob("*.py")
    )
    mark = ROOT / "SOURCE_F3_COMMIT.txt"
    return dict(
        execution_commit=mark.read_text().strip() if mark.exists() else None,
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
        f2_archive=PROTOCOL["f2_archive"],
        mature_state_sha256=hashlib.sha256((F2 / "mature_v1/MATURE_STATE.npz").read_bytes()).hexdigest(),
    )
