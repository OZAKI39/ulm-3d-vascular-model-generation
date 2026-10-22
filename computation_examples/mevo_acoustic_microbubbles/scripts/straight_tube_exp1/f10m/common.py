from pathlib import Path
import subprocess
import yaml
from scripts.straight_tube_exp1.common import sha as sha, read as read, write as write

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/straight_tube_exp1/f10m"
CONFIG = ROOT / "configs/straight_tube_exp1/f10m"
BASELINE = "674482e4f76ecbcff4efe98432fb0eaeb9df5153"


def configs():
    return {p.stem: yaml.safe_load(p.read_text()) for p in sorted(CONFIG.glob("*.yaml"))}


def identity():
    marker = ROOT / "SOURCE_F10M_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    files = [
        *CONFIG.glob("*.yaml"),
        *(ROOT / "scripts/straight_tube_exp1/f10m").glob("*.py"),
        ROOT / "src/mevo_mb/acoustics/kwave_backend.py",
        ROOT / "src/mevo_mb/acoustics/spatial_medium.py",
        ROOT / "src/mevo_mb/coupling/pdms_block.py",
        ROOT / "src/mevo_mb/fast_sim/pdms_engine.py",
    ]
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        files_sha256={str(p.relative_to(ROOT)): sha(p) for p in files},
    )
