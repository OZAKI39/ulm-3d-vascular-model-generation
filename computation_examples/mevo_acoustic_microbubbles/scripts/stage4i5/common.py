"""Versioned Stage I.5 evidence, separate from immutable scientific gold."""
from datetime import datetime, timezone
import json
import importlib.metadata
import os
import platform
import subprocess
import zipfile

from scripts.stage4.common import ROOT, sha, write

BASE = ROOT / "results/stage4/performance_refactor"
GOLD = BASE / "gold_baseline"


def new_run(kind, config):
    out = BASE / kind / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True)
    sources = {str(p.relative_to(ROOT)): sha(p)
               for folder in ("src", "scripts/stage4", "scripts/stage4i5", "configs/stage4")
               for p in (ROOT / folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    write(out / "IDENTITY.json", dict(created_UTC=datetime.now(timezone.utc).isoformat(),
          config=config, source_sha256=sources,
          git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
          working_tree=subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True),
          environment=dict(python=platform.python_version(),platform=platform.platform(),
              packages={name:importlib.metadata.version(name) for name in ["numpy","scipy","vtk","h5py"]},
              threads={name:os.environ.get(name) for name in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]}),
          gold_manifest_sha256=sha(GOLD / "GOLD_BASELINE_MANIFEST.json")))
    with zipfile.ZipFile(out / "source_snapshot.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for path in sources:
            archive.write(ROOT / path, path)
    return out


def finish(out, summary):
    result = dict(summary, run=str(out.relative_to(ROOT)))
    write(out / "SUMMARY.json", result)
    write(out.parent / "STATE.json", result)
    print(json.dumps(result, allow_nan=False), flush=True)
    return result
