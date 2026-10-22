"""F10K isolated paths and scientific execution identity (CPU-safe)."""

from pathlib import Path
import subprocess
from scripts.fast_pzt_f9.common import clean, js, sha, write  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10k"
INPUT = ROOT / "results/fast_pzt/f5/inputs"
BASELINE = "91b09993167625fc103ebd6a026773753a12d3f4"


def identity():
    marker = ROOT / "SOURCE_F10K_COMMIT.txt"
    return dict(
        execution_commit=marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        baseline_commit=BASELINE,
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10k/protocol.json"),
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["scripts/fast_pzt_f10k", "src/mevo_mb/fast_sim", "configs/fast_pzt/f10k"]
            for p in sorted((ROOT / folder).glob("*"))
            if p.is_file()
        },
    )
