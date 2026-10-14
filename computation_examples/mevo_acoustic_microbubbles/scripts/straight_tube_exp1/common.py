from pathlib import Path
import hashlib
import json
import subprocess
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/straight_tube_exp1/f10l"
CONFIG = ROOT / "configs/straight_tube_exp1"
BASELINE = "850fd2df3225351f4fe9564f39cbe72a38662b30"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(v):
    if isinstance(v, dict):
        return {str(k): clean(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [clean(x) for x in v]
    if isinstance(v, np.ndarray):
        return clean(v.tolist())
    if isinstance(v, np.generic):
        return clean(v.item())
    if isinstance(v, float) and not np.isfinite(v):
        return None
    return v


def write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(clean(data), indent=2, allow_nan=False) + "\n")


def read(path):
    return json.loads(Path(path).read_text())


def configs():
    return {p.stem: yaml.safe_load(p.read_text()) for p in sorted(CONFIG.glob("*.yaml"))}


def identity():
    marker = ROOT / "SOURCE_F10L_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        config_sha256={p.name: sha(p) for p in sorted(CONFIG.glob("*.yaml"))},
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["src/mevo_mb", "scripts/straight_tube_exp1"]
            for p in sorted((ROOT / folder).rglob("*.py"))
        },
    )


def validate_saved_acoustic_inputs():
    """Refuse a stale three-source field after any registered configuration change."""
    import h5py
    from mevo_mb.acoustics.kwave_backend import acoustic_identity

    recorded = read(BASE / "INPUT_IDENTITY.json")
    current = {p.name: sha(p) for p in sorted(CONFIG.glob("*.yaml"))}
    if current != recorded["config_sha256"]:
        raise ValueError("Configuration changed after preflight; register a new result directory")
    resolved = read(BASE / "RESOLVED_ACOUSTIC_CONFIG.json")
    with h5py.File(BASE / "acoustic/basis.h5") as h:
        metadata = json.loads(h.attrs["metadata_json"])
    if metadata["acoustic_identity"] != acoustic_identity(resolved):
        raise ValueError("Stored basis does not match the resolved source, medium and grid")
    if sha(BASE / "acoustic/basis.h5") != read(BASE / "acoustic/ACOUSTIC_GATE.json")["basis_sha256"]:
        raise ValueError("Acoustic basis checksum mismatch")
    return True
