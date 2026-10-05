"""Measured crossover and explicit resource fallback; never change force physics."""

import json
from pathlib import Path
import warnings

import numpy as np
from .forces import ClassicalFarFieldBjerknes


class AutoPairBackend:
    def __init__(self, benchmark_path=None, *, crossover=None, **kwargs):
        self.decisions = []
        self.gpu = None
        self.cpu = ClassicalFarFieldBjerknes("numpy", **kwargs)
        self.kwargs = kwargs
        if crossover is None:
            if benchmark_path is None:
                raise ValueError("Auto selection needs benchmark evidence or explicit crossover")
            data = json.loads(Path(benchmark_path).read_text())
            rows = data["rows"]
            if data["status"] != "PASS":
                raise ValueError("Pair benchmark did not pass")
            wins = []
            for n in sorted({r["N"] for r in rows}):
                by = {r["backend"]: r for r in rows if r["N"] == n}
                if (
                    "numpy" in by
                    and "cupy" in by
                    and by["cupy"]["synthetic_step_s"] < by["numpy"]["synthetic_step_s"]
                ):
                    wins.append(n)
            crossover = min(wins) if wins else np.inf
        if crossover <= 0:
            raise ValueError("Positive crossover required")
        self.crossover = crossover
        self.gpu_disabled = False
        self.decisions.append(
            dict(
                crossover=None if not np.isfinite(crossover) else int(crossover),
                benchmark=str(benchmark_path),
            )
        )

    def _gpu(self):
        import cupy as cp

        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError("GPU unavailable")
        return ClassicalFarFieldBjerknes("cupy", **self.kwargs)

    def forces(self, *args):
        if len(args[0]) < self.crossover or self.gpu_disabled:
            return self.cpu.forces(*args)
        try:
            if self.gpu is None:
                self.gpu = self._gpu()
        except (ImportError, RuntimeError) as e:
            return self._fallback(args, e, "GPU_UNAVAILABLE")
        try:
            result = self.gpu.forces(*args)
            return result.get() if hasattr(result, "get") else result
        except MemoryError as e:
            return self._fallback(args, e, "CUDA_OOM_FALLBACK")
        except Exception as e:
            # CuPy OutOfMemoryError is not consistently a built-in MemoryError.
            if type(e).__name__ == "OutOfMemoryError":
                return self._fallback(args, e, "CUDA_OOM_FALLBACK")
            raise

    def _fallback(self, args, error, reason):
        self.gpu_disabled = True
        self.decisions.append(dict(backend="numpy", reason=reason, error=repr(error)))
        warnings.warn(
            reason + ": explicitly falling back to the same NumPy force model", RuntimeWarning, stacklevel=2
        )
        values = [a.get() if hasattr(a, "get") else a for a in args]
        return self.cpu.forces(*values)
