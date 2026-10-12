"""Isolated F10I analysis identity; no imports of GPU-dependent modules."""

from pathlib import Path
import subprocess
from scripts.fast_pzt_f9.common import js, sha, write, clean  # noqa: F401
from scripts.fast_pzt_f10h.analyze import gz, read_steps  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10i"
H = ROOT / "results/fast_pzt/f10h"
BASELINE = "9419837e0ae1ac0a2c6629c63f9a69027b2fef8f"
TRIALS = {
    "DT10": H / "b2_short_v1/O_1000_1102_B2_GEOM_K2_DT10",
    "DT5": H / "b2_dt5_v1/O_1000_1102_B2_GEOM_K2_DT5",
}


def identity():
    marker = ROOT / "SOURCE_F10I_COMMIT.txt"
    return dict(
        execution_commit=marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        baseline_commit=BASELINE,
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10i/protocol.json"),
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["scripts/fast_pzt_f10i", "src/mevo_mb/fast_sim"]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )
