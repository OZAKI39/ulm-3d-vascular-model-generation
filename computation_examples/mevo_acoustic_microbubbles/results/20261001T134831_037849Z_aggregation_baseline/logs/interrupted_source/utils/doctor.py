"""Read-only environment audit; never installs or modifies system components."""
from __future__ import annotations

import importlib.metadata as md
import importlib.util
import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import psutil


def command(args: list[str]) -> dict:
    if not shutil.which(args[0]):
        return {"available": False}
    result = subprocess.run(args, capture_output=True, text=True, timeout=20, check=False)
    return {"returncode": result.returncode, "stdout": result.stdout, "stderr": result.stderr}


def doctor() -> dict:
    packages = {}
    for name in ["k-wave-python", "cupy-cuda13x", "cupy-cuda12x", "numpy", "scipy", "h5py", "PyYAML", "pyvista", "vtk", "matplotlib", "pandas", "pytest", "ruff"]:
        try:
            packages[name] = md.version(name)
        except md.PackageNotFoundError:
            packages[name] = None
    result = {"UTC": datetime.now(timezone.utc).isoformat(), "platform": platform.platform(),
              "os_release": Path("/etc/os-release").read_text(), "python": sys.version,
              "executable": sys.executable, "logical_cpus": os.cpu_count(),
              "affinity_cpus": len(os.sched_getaffinity(0)), "host_RAM_bytes": psutil.virtual_memory().total,
              "packages": packages, "tools": {s: shutil.which(s) for s in ["uv", "conda", "gcc", "g++", "nvcc", "srun", "sbatch", "docker", "apptainer"]}}
    for name, args in {"lscpu": ["lscpu"], "nvidia_smi": ["nvidia-smi"], "cuda_toolkit": ["nvcc", "--version"], "gcc": ["gcc", "--version"], "g++": ["g++", "--version"]}.items():
        result[name] = command(args)
    result["cgroup"] = {p: (Path("/sys/fs/cgroup") / p).read_text().strip() for p in ["memory.max", "memory.current", "cpu.max"] if (Path("/sys/fs/cgroup") / p).exists()}
    if importlib.util.find_spec("cupy"):
        try:
            import cupy as cp
            props = cp.cuda.runtime.getDeviceProperties(0)
            result["cupy"] = {"version": cp.__version__, "runtime_version": cp.cuda.runtime.runtimeGetVersion(),
                               "driver_version": cp.cuda.runtime.driverGetVersion(), "device": props["name"].decode(),
                               "total_memory_bytes": props["totalGlobalMem"], "FP64_test": float(cp.sum(cp.arange(10, dtype=cp.float64)).get())}
        except Exception as error:
            result["cupy"] = {"error": repr(error)}
    return result
