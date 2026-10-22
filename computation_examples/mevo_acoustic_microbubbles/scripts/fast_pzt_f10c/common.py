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
BASE = ROOT / "results/fast_pzt/f10c"
F10A = ROOT / "results/fast_pzt/f10a"


def identity():
    marker = ROOT / "SOURCE_F10C_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit="81fe1206d623f04ea7e3d3a78021251890185c9e",
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10c/protocol.json"),
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["src/mevo_mb/fast_sim", "scripts/fast_pzt_f10c"]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )


def case(key, phase="extended_10_v1"):
    return js(F10A / phase / key / "CONFIG.json")
