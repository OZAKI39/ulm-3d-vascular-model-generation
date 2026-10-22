from pathlib import Path
import subprocess
from scripts.fast_pzt_f9.common import (
    clean as clean,
    write as write,
    js as js,
    sha as sha,
    objects as objects,
)

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10d"
F10C = ROOT / "results/fast_pzt/f10c"
F5 = ROOT / "results/fast_pzt/f5"
BASELINE = "2fc90fd4ad57d16addfe5bed2878c7c3085efd61"


def identity():
    marker = ROOT / "SOURCE_F10D_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10d/protocol.json"),
        inputs=js(F5 / "INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["src/mevo_mb/fast_sim", "scripts/fast_pzt_f10c", "scripts/fast_pzt_f10d"]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )


def case(key):
    return js(ROOT / "results/fast_pzt/f10a/extended_10_v1" / key / "CONFIG.json")
