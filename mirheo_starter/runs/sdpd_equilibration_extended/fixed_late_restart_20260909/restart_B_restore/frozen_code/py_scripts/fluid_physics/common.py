"""Reuse stage-one safe paths, exclusive JSON, checksums. Never import a solver."""
from datetime import datetime, timezone
from pathlib import Path
import hashlib
import json
import os
import subprocess
import yaml
from py_scripts.vessel_geometry.io import PROJECT_ROOT, require_file, sha256_file, read_json, write_json, safe_output

LEGACY = PROJECT_ROOT.parent / "ulm_3D_vascular"


def now():
    return datetime.now(timezone.utc).isoformat()


def fingerprint(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, allow_nan=False, separators=(",", ":")).encode()).hexdigest()


def output(path):
    p = Path(path)
    return safe_output(p if p.is_absolute() else PROJECT_ROOT / p, LEGACY)


def config_read(path):
    path = require_file(Path(path))
    c = yaml.safe_load(path.read_text())
    if c.get("schema_version") != 1:
        raise ValueError("unsupported schema")
    b = c["budget"]
    if not (0 < b["task_limit_s"] <= 600 and 0 < b["campaign_limit_s"] <= 3600 and b["gpu_concurrency"] == 1):
        raise ValueError("budget exceeds user authorization")
    if not c["campaign_id"] or "/" in c["campaign_id"] or c["campaign_id"] in (".", ".."):
        raise ValueError("invalid explicit campaign id")
    c["_config_path"] = str(path)
    c["_config_sha256"] = sha256_file(path)
    for key in ("data_root", "runs_root", "review_root"):
        output(c[key])
    return c


def atomic_state(path, value):
    """Only mutable scheduler state is replaced; immutable artifacts use write_json."""
    path = output(path)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    with tmp.open("x") as f:
        json.dump(value, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def resource_snapshot():
    gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.total,memory.used,memory.free", "--format=csv,noheader,nounits"],
                         capture_output=True, text=True, timeout=5)
    if gpu.returncode:
        raise RuntimeError(gpu.stderr)
    rows = []
    for line in gpu.stdout.strip().splitlines():
        name, total, used, free = [s.strip() for s in line.split(",")]
        rows.append(dict(name=name, total_MiB=int(total), used_MiB=int(used), free_MiB=int(free)))
    mem = {s.split(":")[0]: int(s.split()[1])*1024 for s in Path("/proc/meminfo").read_text().splitlines() if s.startswith(("MemTotal:", "MemAvailable:", "SwapTotal:", "SwapFree:"))}
    return {"recorded_at": now(), "gpus": rows, "host_memory_bytes": mem,
            "gpu_memory_semantics": "device-wide sampled values, includes other applications; not exact process peak"}


def environment_identity():
    lib = require_file(PROJECT_ROOT / ".venv/lib/python3.12/site-packages/libmirheo.cpython-312-x86_64-linux-gnu.so")
    vendor = PROJECT_ROOT / "vendor/Mirheo"
    commit = subprocess.run(["git", "-C", str(vendor), "rev-parse", "HEAD"], capture_output=True, text=True, check=True).stdout.strip()
    paths = [lib, PROJECT_ROOT/"scripts/activate_mirheo.sh", PROJECT_ROOT/"scripts/run_official_case.sh", PROJECT_ROOT/"metadata/CMakeCache.txt"]
    paths += [p for p in (PROJECT_ROOT/"metadata/compat_patches").rglob("*") if p.is_file()]
    paths += [vendor/p for p in (
        "tests/flow/double_poiseuille.py", "tests/stress/pressure.py", "tests/doc_scripts/basic.py", "tests/doc_scripts/walls.py", "LICENSE",
        "src/mirheo/core/interactions/pairwise/kernels/dpd.h", "src/mirheo/core/interactions/pairwise/kernels/stress_wrapper.h",
        "src/mirheo/core/integrators/forcing_terms/periodic_poiseuille.h", "src/mirheo/plugins/stats.cu", "src/mirheo/plugins/virial_pressure.cu",
        "src/mirheo/bindings/plugins.cpp", "src/mirheo/bindings/particle_vectors.cpp")]
    cache = (PROJECT_ROOT/"metadata/CMakeCache.txt").read_text()
    precision = "double" if "MIR_DOUBLE_PRECISION:BOOL=ON" in cache else "single"
    return {"mirheo_commit": commit, "library": str(lib), "library_sha256": sha256_file(lib), "precision": precision,
            "python": str(PROJECT_ROOT/".venv/bin/python"), "mpi": "/usr/bin/mpirun.openmpi",
            "sha256": {str(p): sha256_file(require_file(p)) for p in paths}, "patches_preserved": True}
