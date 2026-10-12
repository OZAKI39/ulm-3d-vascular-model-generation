"""Isolated F10H paths and reproducible input identity; CPU-safe imports."""

from pathlib import Path
import gzip
import json
import subprocess
from scripts.fast_pzt_f9.common import clean, js, sha, write  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10h"
BASELINE = "8df3f063af6e0099ac630cd350bdd8a6a22ed4bd"
CASES = ["O_1000_1102", "prepare_1104"]


def state(case):
    p = ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0/CONTACT_TRACE.json.gz"
    with gzip.open(p, "rt") as s:
        data = json.load(s)
    return data, data["traces"][-1], p


def identity():
    p = ROOT / "SOURCE_F10H_COMMIT.txt"
    commit = (
        p.read_text().strip()
        if p.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10h/protocol.json"),
        traces={c: sha(state(c)[2]) for c in CASES},
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["src/mevo_mb/fast_sim", "scripts/fast_pzt_f10h"]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )
