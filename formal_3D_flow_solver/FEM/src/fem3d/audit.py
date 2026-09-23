"""Small independent provenance utilities; no reference-project imports."""

import hashlib
import json
import os
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def command(args, timeout=30):
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=timeout)
        return {"command": args, "returncode": p.returncode,
                "stdout": p.stdout, "stderr": p.stderr}
    except (OSError, subprocess.TimeoutExpired) as exc:
        return {"command": args, "returncode": None, "error": str(exc)}


def environment():
    commands = [["uname", "-a"], ["cat", "/etc/os-release"], ["lscpu"],
                ["free", "-h"], ["df", "-h"], ["nvidia-smi"],
                ["nvcc", "--version"], ["mpirun", "--version"],
                ["python", "--version"], [sys.executable, "--version"]]
    packages = {}
    for name in ("dolfinx", "petsc4py", "mpi4py", "basix", "ufl", "ffcx", "gmsh",
                 "numpy", "scipy", "pyvista", "vtk", "pytest", "yaml", "matplotlib"):
        code = ("import importlib; m=importlib.import_module(" + repr(name) + "); "
                "print(getattr(m, '__version__', 'version attribute unavailable'))")
        packages[name] = command([sys.executable, "-B", "-c", code])
    petsc = command([sys.executable, "-B", "-c", "from petsc4py import PETSc; "
                     "print('PETSc',PETSc.Sys.getVersion()); "
                     "print('CUDA',PETSc.Sys.hasExternalPackage('cuda')); "
                     "print('scalar',PETSc.ScalarType); "
                     "print('MUMPS',PETSc.Sys.hasExternalPackage('mumps'))"])
    mpi = command([sys.executable, "-B", "-c", "from mpi4py import MPI; "
                   "print(MPI.Get_library_version())"])
    return {"timestamp": timestamp(), "hostname": platform.node(),
            "python_executable": sys.executable, "python_version": sys.version,
            "commands": [command(c) for c in commands], "packages": packages,
            "petsc": petsc, "mpi_library": mpi}


def reference_snapshot(roots):
    files = {}
    for root in roots:
        root = Path(root).resolve()
        for base, dirs, names in os.walk(root, followlinks=False):
            for name in sorted(dirs + names):
                p = Path(base) / name
                key = str(p)
                s = p.lstat()
                if p.is_symlink():
                    files[key] = {"symlink": os.readlink(p), "mode": s.st_mode}
                elif p.is_file():
                    files[key] = {"size": s.st_size, "mtime_ns": s.st_mtime_ns,
                                  "mode": s.st_mode, "sha256": sha256(p)}
                elif p.is_dir():
                    files[key] = {"directory": True, "mode": s.st_mode}
    return {"timestamp": timestamp(), "roots": list(map(str, roots)), "files": files,
            "policy": "All regular files including .git, caches, inputs and results; symlinks recorded without following; atime excluded."}
