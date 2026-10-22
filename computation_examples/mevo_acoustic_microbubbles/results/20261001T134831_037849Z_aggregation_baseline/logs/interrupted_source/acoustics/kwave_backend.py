"""Verified k-wave-python >=0.6 unified API and linear continuous-wave extraction."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import subprocess
import threading
import time
from pathlib import Path

import numpy as np

from mevo_mb.acoustics.basis import AcousticBasisDatabase
from mevo_mb.control.pzt import PZTArray
from mevo_mb.coupling.geometry import grid_axes
from mevo_mb.io.artifacts import git_commit, write_json
from mevo_mb.utils.config import fingerprint


def acoustic_identity(config: dict) -> str:
    """Controls excluded: geometry, propagation and extraction define reusable H."""
    elements = [{k: v for k, v in e.items() if k not in ("amplitude_Pa", "phase_rad", "enabled")} for e in config["pzt"]]
    return fingerprint({"grid": config["grid"], "medium": config["medium"], "acoustics": config["acoustics"], "elements": elements})


def extract_complex(data: np.ndarray, times_s: np.ndarray, frequency_Hz: float) -> np.ndarray:
    """Exact coherent DFT on integer cycles, absolute phase Re[P exp(+i wt)]."""
    kernel = np.exp(-2j * np.pi * frequency_Hz * times_s)
    out = np.empty(data.shape[0], np.complex64)
    for start in range(0, len(out), 4096):
        out[start:start+4096] = (2.0 / len(times_s)) * (data[start:start+4096] @ kernel)
    return out


class GPUTelemetry:
    """Sample whole-device VRAM; report attribution limits explicitly."""
    def __init__(self):
        self.stop_event = threading.Event()
        self.samples = []
        self.errors = []

    def _sample(self):
        while not self.stop_event.is_set():
            try:
                r = subprocess.run(["nvidia-smi", "--query-gpu=name,memory.used,utilization.gpu", "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=5, check=True)
                name, memory, utilization = r.stdout.strip().splitlines()[0].rsplit(",", 2)
                self.samples.append({"monotonic_s": time.monotonic(), "gpu": name, "used_MiB": float(memory), "utilization_percent": float(utilization)})
            except (OSError, subprocess.SubprocessError, ValueError) as error:
                self.errors.append(repr(error))
            self.stop_event.wait(0.5)

    def __enter__(self):
        self.thread = threading.Thread(target=self._sample, daemon=True)
        self.thread.start()
        return self

    def __exit__(self, *args):
        self.stop_event.set()
        self.thread.join(timeout=6)


class KWaveCW:
    """One prepared source geometry reusable across basis and direct excitations."""
    def __init__(self, config: dict, output: Path):
        from kwave.kgrid import kWaveGrid
        from kwave.utils.kwave_array import kWaveArray

        self.config, self.output = config, output
        a, grid = config["acoustics"], config["grid"]
        if a.get("BonA") is not None:
            raise ValueError("Linear basis requires BonA unset")
        self.axes = grid_axes(grid["shape"], grid["spacing_m"])
        self.kgrid = kWaveGrid(grid["shape"], grid["spacing_m"])
        f = a["frequency_Hz"]
        c = config["medium"]["sound_speed_m_s"]
        self.samples_per_period = int(np.ceil(c / (a["CFL"] * min(grid["spacing_m"]) * f)))
        self.dt_s = 1.0 / (self.samples_per_period * f)
        diagonal = np.linalg.norm(np.array(grid["shape"]) * grid["spacing_m"])
        settle_periods = int(np.ceil(a["transit_multiples"] * diagonal / c * f))
        self.Nt = int((a["ramp_periods"] + settle_periods + a["recording_periods"]) * self.samples_per_period)
        self.kgrid.setTime(self.Nt, self.dt_s)
        self.record_samples = int(a["recording_periods"] * self.samples_per_period)
        self.record_start = self.Nt - self.record_samples
        self.times = np.arange(self.Nt) * self.dt_s
        self.roi_indices = tuple(np.where(np.abs(axis) <= half + 1e-14)[0] for axis, half in zip(self.axes, a["roi_half_width_m"], strict=True))
        self.roi_axes = tuple(axis[idx] for axis, idx in zip(self.axes, self.roi_indices, strict=True))
        self.sensor_mask = np.zeros(grid["shape"], bool)
        self.sensor_mask[np.ix_(*self.roi_indices)] = True
        self.roi_shape = tuple(len(x) for x in self.roi_axes)
        if min(self.roi_shape) < 3:
            raise ValueError("ROI needs at least three samples per axis")
        array = PZTArray.from_config(config)
        self.source_weights = []
        self.source_mask = np.zeros(grid["shape"], bool)
        kw_array = kWaveArray(bli_tolerance=a["bli_tolerance"], upsampling_rate=a["upsampling_rate"], single_precision=True)
        for e in array.elements:
            if e.geometry != "disc":
                raise NotImplementedError("Initial source adapter supports 3D disc elements")
            normal = np.asarray(e.normal, float)
            normal /= np.linalg.norm(normal)
            kw_array.add_disc_element(e.position_m, e.size_m[0], (np.asarray(e.position_m) + normal * e.size_m[0]).tolist())
        cache = output / "source_weights.npz"
        for i in range(len(array.elements)):
            w = np.asarray(kw_array.get_element_grid_weights(self.kgrid, i), np.float32)
            self.source_weights.append(w)
            self.source_mask |= w != 0
        self.weights = np.stack([w[self.source_mask] for w in self.source_weights])
        self.source_weights = None
        pml = a["pml_size"]
        occupied = np.argwhere(self.source_mask | self.sensor_mask)
        if a["pml_inside"] and (np.any(occupied.min(axis=0) < pml) or np.any(occupied.max(axis=0) >= np.asarray(grid["shape"]) - pml)):
            raise ValueError("Source BLI support or ROI intersects PML; enlarge domain or move source")
        np.savez_compressed(cache, weights=self.weights, source_indices_C=np.flatnonzero(self.source_mask))
        self.weights_sha256 = hashlib.sha256(cache.read_bytes()).hexdigest()

    def run(self, control_Pa: np.ndarray, label: str) -> tuple[np.ndarray, np.ndarray, dict]:
        from kwave.kmedium import kWaveMedium
        from kwave.ksensor import kSensor
        from kwave.ksource import kSource
        from kwave.kspaceFirstOrder import kspaceFirstOrder

        c = self.config
        a, med = c["acoustics"], c["medium"]
        folder = self.output / label
        folder.mkdir(exist_ok=False)
        ramp_end = a["ramp_periods"] / a["frequency_Hz"]
        ramp = 0.5 * (1-np.cos(np.pi * np.minimum(self.times / ramp_end, 1.0)))
        harmonics = np.real(np.asarray(control_Pa)[:, None] * np.exp(2j*np.pi*a["frequency_Hz"]*self.times)) * ramp
        source = kSource()
        source.p_mask = self.source_mask
        source.p = np.asarray(self.weights.T @ harmonics.astype(np.float32), np.float32)
        source.p_mode = "additive"
        sensor = kSensor(mask=self.sensor_mask)
        sensor.record = ["p"]
        sensor.record_start_index = self.record_start + 1  # public k-Wave convention: 1-based
        medium = kWaveMedium(sound_speed=med["sound_speed_m_s"], density=med["density_kg_m3"])
        backend, device = a["backend"], a["device"]
        failure = None
        start = time.perf_counter()
        with GPUTelemetry() as telemetry:
            try:
                result = kspaceFirstOrder(self.kgrid, medium, source, sensor, backend=backend, device=device,
                                          pml_size=a["pml_size"], pml_alpha=a["pml_alpha"], pml_inside=a["pml_inside"],
                                          smooth_p0=False, data_path=str(folder), quiet=True, num_threads=1,
                                          dtype="float32" if backend == "python" else None)
            except (OSError, RuntimeError, subprocess.SubprocessError) as error:
                failure = repr(error)
                write_json(folder / "CPP_FAILURE.json", {"exception": failure, "backend": backend, "device": device})
                if not (backend == "cpp" and device == "gpu" and a["allow_python_gpu_fallback"]):
                    raise
                backend = "python"
                result = kspaceFirstOrder(self.kgrid, medium, source, sensor, backend=backend, device="gpu", dtype="float32",
                                          pml_size=a["pml_size"], pml_alpha=a["pml_alpha"], pml_inside=a["pml_inside"], smooth_p0=False, quiet=True)
        elapsed = time.perf_counter() - start
        raw = result["p"]
        if hasattr(raw, "get"):
            raw = raw.get()
        raw = np.asarray(raw)
        expected = (int(self.sensor_mask.sum()), self.record_samples)
        if raw.shape != expected or not np.isfinite(raw).all():
            raise RuntimeError(f"Unexpected sensor output {raw.shape}; expected finite {expected}")
        times = self.times[self.record_start:]
        p = extract_complex(raw, times, a["frequency_Hz"]).reshape(self.roi_shape)
        half = len(times)//2
        p_first = extract_complex(raw[:, :half], times[:half], a["frequency_Hz"])
        p_last = extract_complex(raw[:, half:], times[half:], a["frequency_Hz"])
        steady_delta = float(np.linalg.norm(p_last-p_first)/max(np.linalg.norm(p_last), np.finfo(float).tiny))
        gradient = np.stack(np.gradient(p, *self.roi_axes, edge_order=2), axis=-1).astype(np.complex64)
        np.savez_compressed(folder / "field.npz", P_complex=p, gradP_complex=gradient)
        meta = {"k_wave_python_version": importlib.metadata.version("k-wave-python"), "backend": backend, "device": device,
                "fallback_failure": failure, "GPU": telemetry.samples[0]["gpu"] if telemetry.samples else None,
                "sampled_device_peak_VRAM_MiB": max((s["used_MiB"] for s in telemetry.samples), default=None),
                "VRAM_scope": "whole-device sampled peak, includes other processes; not a guaranteed allocation upper bound",
                "wall_clock_s": elapsed, "frequency_Hz": a["frequency_Hz"], "source_control_real_Pa": np.real(control_Pa).tolist(),
                "source_control_imag_Pa": np.imag(control_Pa).tolist(), "grid_dimensions": c["grid"]["shape"],
                "grid_spacing_m": c["grid"]["spacing_m"], "CFL_requested": a["CFL"],
                "CFL_actual": med["sound_speed_m_s"]*self.dt_s/min(c["grid"]["spacing_m"]),
                "dt_s": self.dt_s, "Nt": self.Nt, "samples_per_period": self.samples_per_period,
                "sensor_count": int(self.sensor_mask.sum()), "recorded_samples": self.record_samples,
                "record_start_index_1based": self.record_start+1, "PML": {k: a[k] for k in ["pml_size", "pml_inside", "pml_alpha"]},
                "medium": med, "BonA": None, "git_commit": git_commit(), "acoustic_identity": acoustic_identity(c),
                "source_weight_sha256": self.weights_sha256, "source_adapter": "kWaveArray disc BLI weights; explicit C-order mask rows",
                "harmonic_convention": "Re[P exp(+i omega t)]; source cos(omega t + phase)",
                "local_ROI_PNP_Pa": float(np.abs(p).max()), "recording_half_window_relative_L2": steady_delta,
                "output_bytes": sum(f.stat().st_size for f in folder.iterdir() if f.is_file())}
        write_json(folder / "metadata.json", meta)
        write_json(folder / "GPU_TELEMETRY.json", {"samples": telemetry.samples, "errors": telemetry.errors})
        print(json.dumps({"stage": label, "elapsed_s": elapsed, "ROI_peak_Pa": meta["local_ROI_PNP_Pa"], "steady_delta": steady_delta}), flush=True)
        return p, gradient, meta


def precompute(config: dict, run: Path) -> Path:
    adapter = KWaveCW(config, run / "acoustic")
    n = len(config["pzt"])
    amplitude = config["acoustics"]["basis_source_amplitude_Pa"]
    pressures, gradients, metadata = [], [], []
    for i in range(n):
        u = np.zeros(n, complex)
        u[i] = amplitude
        p, g, meta = adapter.run(u, f"basis_{i:02d}")
        pressures.append(p/amplitude)
        gradients.append(g/amplitude)
        metadata.append(meta)
    info = {"frequency_Hz": config["acoustics"]["frequency_Hz"], "BonA": None,
            "acoustic_identity": acoustic_identity(config), "basis_source_amplitude_Pa": amplitude,
            "source_ids": [e["id"] for e in config["pzt"]], "simulations": metadata}
    path = run / "acoustic/basis.h5"
    AcousticBasisDatabase(adapter.roi_axes, pressures, gradients, info).save_basis(path)
    write_json(run / "metadata.json", {"status": "BASIS_COMPLETE", **info})
    return path
