from pathlib import Path
import subprocess
import yaml
from scripts.straight_tube_exp1.common import sha as sha, read as read, write as write

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/straight_tube_exp1/f10m_r3"
CONFIG = ROOT / "configs/straight_tube_exp1/f10m_r3"
BASELINE = "a494cdddd71f0f7bbd09f4ef37e445a8b89d8548"


def configs():
    return {p.stem: yaml.safe_load(p.read_text()) for p in sorted(CONFIG.glob("*.yaml"))}


def identity():
    marker = ROOT / "SOURCE_F10M_R3_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    files = [*CONFIG.glob("*.yaml"), *(ROOT / "scripts/straight_tube_exp1/f10m_r3").glob("*.py")]
    files += [
        ROOT / p
        for p in [
            "src/mevo_mb/coupling/pdms_cavity.py",
            "src/mevo_mb/fast_sim/pdms_cavity_engine.py",
            "src/mevo_mb/acoustics/kwave_backend.py",
            "src/mevo_mb/acoustics/spatial_medium.py",
            "src/mevo_mb/fast_sim/f10c_engine.py",
        ]
        if (ROOT / p).exists()
    ]
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        files_sha256={str(p.relative_to(ROOT)): sha(p) for p in files},
    )
