from pathlib import Path
import subprocess
import numpy as np
from scripts.straight_tube_exp1.f10q.common import write, read, sha  # noqa: F401

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results/straight_tube_exp1/f10s"
OLD = ROOT / "results/straight_tube_exp1/f10n"
Q = ROOT / "results/straight_tube_exp1/f10q"
R = ROOT / "results/straight_tube_exp1/f10r"
CONFIG = ROOT / "configs/straight_tube_exp1/f10s/protocol.json"


def policy():
    p = read(ROOT / "configs/straight_tube_exp1/f10q/protocol.json")
    p.update(read(ROOT / "configs/straight_tube_exp1/f10q/algorithm_revision_1.json"))
    return p


def identity():
    marker = ROOT / "SOURCE_F10S_COMMIT.txt"
    commit = marker.read_text().strip() if marker.exists() else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    paths = list((ROOT / "scripts/straight_tube_exp1/f10s").glob("*.py")) + [
        ROOT / "src/mevo_mb/control/scheduled_pzt.py", CONFIG,
        ROOT / "src/mevo_mb/fast_sim/water_sequential.py",
        ROOT / "src/mevo_mb/fast_sim/sequential_hard_sphere.py"]
    return dict(execution_commit=commit, source_sha256={str(p.relative_to(ROOT)): sha(p) for p in paths if p.exists()})


def context():
    from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis, WaterDiscs
    from mevo_mb.fast_sim.models import ResponseSection
    from mevo_mb.fast_sim.nearfield import F2Physics
    from mevo_mb.coupling.pdms_cavity import apparatus
    physics = F2Physics(radius_m=policy()["radius_m"], frequency_Hz=2.5e6)
    geometry = read(OLD / "CONFIG_SNAPSHOT.json")["geometry"]
    cavity, _ = apparatus(geometry)
    return dict(physics=physics, response=ResponseSection(OLD / "inputs/response.h5", physics),
                basis=CachedWaterBasis(OLD / "cache/coarse.h5"), cavity=cavity,
                discs=WaterDiscs(geometry["pzt"]), x0=np.load(OLD / "initial_states.npz")["aggregate_m"])


def original(axis):
    return read(Q / f"b_navigation_v1/runs/NAV_{axis}_DT1/INPUT.json")["control"]
