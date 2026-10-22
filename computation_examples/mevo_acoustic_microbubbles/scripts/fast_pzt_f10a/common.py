from pathlib import Path
import subprocess
from scripts.fast_pzt_f9.common import clean as clean, write as write, js, sha, objects as objects

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10a"


def identity():
    marker = ROOT / "SOURCE_F10A_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit="4bf369f4f65efb5909e9aedb69e460eb692584da",
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10a/protocol.json"),
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["src/mevo_mb/fast_sim", "scripts/fast_pzt_f10a"]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )
