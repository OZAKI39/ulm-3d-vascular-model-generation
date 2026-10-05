"""Same-input backend benchmark; isolated runtime libraries, no global symlinks."""

from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess

import numpy as np

from mevo_mb.acoustics.kwave_backend import KWaveCW
from mevo_mb.utils.config import load_config


def main():
    root = Path.cwd()
    out = (
        root
        / "results/stage2/performance"
        / ("kwave_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ"))
    )
    out.mkdir(parents=True)
    lib = root / "third_party/build/kwave_runtime/usr/lib/x86_64-linux-gnu"
    os.environ["LD_LIBRARY_PATH"] = str(lib) + ":" + os.environ.get("LD_LIBRARY_PATH", "")
    import kwave

    binary = Path(kwave.__file__).parent / "bin/linux/kspaceFirstOrder-CUDA"
    evidence = {}
    for name, args in (("ldd", ["ldd", str(binary)]), ("binary_version", [str(binary), "--version"])):
        z = subprocess.run(args, capture_output=True, text=True, timeout=30)
        (out / (name + ".txt")).write_text(z.stdout + z.stderr)
        evidence[name] = z.returncode
    results = {}
    fields = {}
    for size in (32, 64):
        c = load_config("configs/smoke.yaml")
        c["grid"]["shape"] = [size] * 3
        for backend in ("python", "cpp"):
            c["acoustics"].update(backend=backend, device="gpu", allow_python_gpu_fallback=False)
            path = out / f"{size}_{backend}"
            path.mkdir()
            try:
                model = KWaveCW(c, path)
                p, g, meta = model.run(np.array([1000.0 + 0j]), "run")
                fields[size, backend] = (p, g)
                results[f"{size}_{backend}"] = meta
            except Exception as e:
                results[f"{size}_{backend}"] = {"error": repr(e)}
        if (size, "cpp") in fields and (size, "python") in fields:
            p, g = fields[size, "python"]
            q, h = fields[size, "cpp"]
            results[f"{size}_comparison"] = {
                "pressure_relative_L2": float(np.linalg.norm(p - q) / np.linalg.norm(p)),
                "gradient_relative_L2": float(np.linalg.norm(g - h) / np.linalg.norm(g)),
                "speedup_python_over_cpp": results[f"{size}_python"]["wall_clock_s"]
                / results[f"{size}_cpp"]["wall_clock_s"],
            }
    summary = dict(
        status="BENCHMARK_COMPLETE",
        results=results,
        dependency_checks=evidence,
        package_version=importlib.metadata.version("k-wave-python"),
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        isolated_library_path=str(lib),
        hardware="results/stage2/system/SERVER_DOCTOR.json",
        random_seed=None,
        notes=[
            "No baseline acoustic cache recomputed; new backend equivalence fixtures only",
            "Full time record in C++ can have much higher IO/RAM cost; benchmark size limited to 64³",
            "Backend timings include initial compilation/loading and file IO; not a universal speedup",
            "Device peak is sampled whole-device VRAM, includes other allocations",
        ],
    )
    (out / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out.parent / "KWAVE_STATE.json").write_text(
        json.dumps(dict(summary, run=str(out.relative_to(root))), indent=2) + "\n"
    )
    print({k: v for k, v in results.items() if "comparison" in k}, flush=True)


if __name__ == "__main__":
    main()
