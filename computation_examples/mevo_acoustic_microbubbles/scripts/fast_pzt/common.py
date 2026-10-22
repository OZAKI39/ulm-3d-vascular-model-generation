from pathlib import Path
import json
import hashlib
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f1"


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
        (ROOT / "scripts/fast_pzt").glob("*.py")
    )
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def inputs():
    return BASE / "inputs/response.h5", BASE / "inputs/basis.h5"
