from pathlib import Path
import subprocess
import numpy as np
from scripts.straight_tube_exp1.common import read, write as json_write, sha  # noqa: F401

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/straight_tube_exp1/f10p"
OLD = ROOT / "results/straight_tube_exp1/f10n"
CONFIG = ROOT / "configs/straight_tube_exp1/f10p/protocol.json"


def identity():
    marker = ROOT / "SOURCE_F10P_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    files = list((ROOT / "scripts/straight_tube_exp1/f10p").glob("*.py")) + list(CONFIG.parent.glob("*.json"))
    files += [
        ROOT / p
        for p in [
            "src/mevo_mb/fast_sim/hard_sphere.py",
            "src/mevo_mb/fast_sim/water_hard_sphere.py",
            "src/mevo_mb/fast_sim/water_b2.py",
            "src/mevo_mb/fast_sim/f10h_contact.py",
            "src/mevo_mb/fast_sim/f10g_contact.py",
            "src/mevo_mb/fast_sim/water_control_engine.py",
            "src/mevo_mb/fast_sim/pdms_cavity_engine.py",
            "src/mevo_mb/fast_sim/vessel_kernels.py",
            "src/mevo_mb/acoustics/water_disc_proxy.py",
        ]
    ]
    return dict(execution_commit=commit, source_sha256={str(p.relative_to(ROOT)): sha(p) for p in files})


def event(name, duration):
    p = OLD / "campaign_v1/runs" / name / "INPUT.json"
    e = read(p)["events"][0].copy()
    e["duration_s"] = duration
    return e


def write(path, value):
    def encode(v):
        if isinstance(v, (complex, np.complexfloating)):
            return dict(real=float(v.real), imag=float(v.imag))
        if isinstance(v, np.ndarray):
            return encode(v.tolist())
        if isinstance(v, dict):
            return {str(k): encode(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [encode(x) for x in v]
        return v

    json_write(path, encode(value))
