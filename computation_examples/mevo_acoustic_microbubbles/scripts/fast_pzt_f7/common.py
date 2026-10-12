from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f7"
F5 = ROOT / "results/fast_pzt/f5"
INPUT = F5 / "inputs"
BASELINE = "a8324747da00ae00d4c7e194aa81a3677289374d"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [clean(v) for v in value]
    if isinstance(value, np.ndarray):
        return clean(value.tolist())
    if isinstance(value, np.generic):
        return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False) + "\n")


def js(path):
    return json.loads(Path(path).read_text())


def identity():
    marker = ROOT / "SOURCE_F7_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ("src/mevo_mb/fast_sim", "scripts/fast_pzt_f6", "scripts/fast_pzt_f7")
            for p in sorted((ROOT / folder).glob("*.py"))
        },
        inputs=js(F5 / "INPUT_MANIFEST.json"),
    )


def objects():
    from mevo_mb.fast_sim.nearfield import F2Physics
    from mevo_mb.fast_sim.models import Basis, ResponseSection
    from mevo_mb.fast_sim.f4_engine import GPUVessel
    from mevo_mb.fast_sim.flow_gpu import GPUFlow

    p = F2Physics()
    return (
        p,
        Basis(INPUT / "f4/basis.h5", p.frequency_Hz),
        ResponseSection(INPUT / "f4/response.h5", p),
        GPUVessel(INPUT / "f4/prepared_geometry.npz", INPUT / "f4/labelled_lumen_m.vtp"),
        GPUFlow(INPUT / "result_121.vtu"),
    )


MODES = {"M0": "f5_local", "M1": "diagonal_stokes", "M2": "rpy_free_space"}
