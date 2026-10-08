from pathlib import Path
import hashlib
import json
from scripts.fast_pzt_f2.common import write as write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f5"
INPUT = BASE / "inputs"
F4_ARCHIVE = "e247fa0c61f668a9b0b2e26c909dd687eedde389"


def js(path):
    return json.loads(Path(path).read_text())


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity():
    marker = ROOT / "SOURCE_F5_COMMIT.txt"
    files = sorted((ROOT / "src/mevo_mb/fast_sim").glob("*.py")) + sorted(
        (ROOT / "scripts/fast_pzt_f5").glob("*.py")
    )
    return dict(
        execution_commit=marker.read_text().strip() if marker.exists() else None,
        base_f4_archive=F4_ARCHIVE,
        source_sha256={str(p.relative_to(ROOT)): sha(p) for p in files},
    )
