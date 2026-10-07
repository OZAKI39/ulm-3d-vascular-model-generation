from pathlib import Path
import json
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f2"
F1 = ROOT / "results/fast_pzt/f1"
PROTOCOL = json.loads((ROOT / "configs/fast_pzt/f2/protocol.json").read_text())


def write(path, data):
    def convert(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, np.generic):
            return x.item()
        if isinstance(x, Path):
            return str(x)
        raise TypeError(type(x).__name__)

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, default=convert, allow_nan=False) + "\n")


def source_identity():
    paths = sorted((ROOT / "src/mevo_mb/fast_sim").glob("*.py")) + sorted(
        (ROOT / "scripts/fast_pzt_f2").glob("*.py")
    )
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def quality_pass(q, p):
    limits = PROTOCOL["quality"]
    return (
        q["error_bits"] == 0
        and q["max_accepted_penetration_m"] <= limits["max_accepted_penetration_over_radius"] * p.radius_m
        and q["max_contact_pass_correction_m"]
        <= limits["max_contact_pass_correction_over_radius"] * p.radius_m
    )
