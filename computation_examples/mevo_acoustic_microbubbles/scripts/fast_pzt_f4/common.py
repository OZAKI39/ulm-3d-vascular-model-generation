from pathlib import Path
import hashlib
import json
from scripts.fast_pzt_f2.common import write as write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f4"


def js(path):
    return json.loads(Path(path).read_text())


def identity():
    mark = ROOT / "SOURCE_F4_COMMIT.txt"
    files = sorted((ROOT / "src/mevo_mb/fast_sim").glob("*.py")) + sorted(
        (ROOT / "scripts/fast_pzt_f4").glob("*.py")
    )
    return dict(
        execution_commit=mark.read_text().strip() if mark.exists() else None,
        base_f3_archive="3eeb94eb885f5df4f0e214ba9b7904dbdecec744",
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
    )
