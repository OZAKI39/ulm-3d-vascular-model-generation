from pathlib import Path
import subprocess
from scripts.fast_pzt_f9.common import clean, write, js, sha, objects  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10f"
HISTORY = ROOT / "results/fast_pzt/f10e"
BASELINE = "741ad7d8683bfbb2fc9a270d9e92a218dead5e54"


def identity():
    marker = ROOT / "SOURCE_F10F_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10f/protocol.json"),
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in [
                "src/mevo_mb/fast_sim",
                "scripts/fast_pzt_f10c",
                "scripts/fast_pzt_f10d",
                "scripts/fast_pzt_f10e",
                "scripts/fast_pzt_f10f",
            ]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )
