"""File-based wrapper of unmodified APECSS v1.8. No radial ODE implemented here."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.interpolate import PchipInterpolator

from .response import BubbleResponseModel

APECSS_COMMIT = "a0ee071ca247b5489e997d90e83b5f53c40e4ec7"


@dataclass(frozen=True)
class RadialCase:
    radius_m: float
    frequency_Hz: float
    pressure_Pa: float
    elasticity_N_m: float
    shell_viscosity_kg_s: float
    sigma_initial_N_m: float
    density_kg_m3: float = 1000.0
    sound_speed_m_s: float = 1500.0
    viscosity_Pa_s: float = 0.002
    sigma_water_N_m: float = 0.072
    polytropic_exponent: float = 1.07
    ambient_pressure_Pa: float = 101300.0
    cycles: int = 48
    samples_per_cycle: int = 512
    tolerance: float = 1e-10
    provenance: str = "Unspecified; caller must provide parameter source before scientific use"

    def __post_init__(self):
        for name, v in asdict(self).items():
            if name == "provenance":
                continue
            if not np.isfinite(v) or v < 0:
                raise ValueError(f"Invalid radial parameter {name}")
        for name in (
            "radius_m",
            "frequency_Hz",
            "density_kg_m3",
            "sound_speed_m_s",
            "elasticity_N_m",
            "polytropic_exponent",
            "ambient_pressure_Pa",
            "tolerance",
        ):
            if getattr(self, name) == 0:
                raise ValueError(f"Positive {name} required")
        if not 0 <= self.sigma_initial_N_m <= self.sigma_water_N_m:
            raise ValueError("Initial surface tension outside clipped shell law")
        if self.cycles < 16 or self.samples_per_cycle < 64:
            raise ValueError("At least 16 cycles and 64 samples per cycle required")


def complex_harmonics(time_s, values, frequency_Hz, harmonic_count=8):
    """Peak phasors for Re[X_h exp(i h omega t)], on complete uniform cycles.

    The endpoint is excluded. Input must have an integer number of periods.
    Units unchanged; DC is excluded, so mean-volume offset cannot leak.
    """
    t, y = np.asarray(time_s), np.asarray(values)
    dt = np.diff(t)
    if len(t) < 4 or y.shape != t.shape or not np.isfinite(y).all():
        raise ValueError("Finite one-dimensional samples required")
    if not np.allclose(dt, dt[0], rtol=1e-7, atol=0) or dt[0] <= 0:
        raise ValueError("Uniform increasing times required")
    cycles = len(t) * dt[0] * frequency_Hz
    if not np.isclose(cycles, round(cycles), atol=1e-6, rtol=0):
        raise ValueError("Harmonic window must span integer cycles, endpoint excluded")
    h = np.arange(1, harmonic_count + 1)
    return 2 * np.mean((y - y.mean())[:, None] * np.exp(-2j * np.pi * frequency_Hz * t[:, None] * h), axis=0)


def analyze_cycles(data, case: RadialCase, retained_cycles=8, harmonic_count=8):
    t, radius, speed = data[:, 1], data[:, 3], data[:, 4]
    tend = case.cycles / case.frequency_Hz
    if not np.isfinite(data).all() or np.min(radius) <= 0 or np.any(np.diff(t) <= 0):
        raise ValueError("Invalid APECSS trajectory")
    if abs(t[-1] - tend) > tend * 1e-10:
        raise ValueError("Incomplete APECSS trajectory")
    n = case.samples_per_cycle
    q = (case.cycles - retained_cycles + np.arange(retained_cycles * n) / n) / case.frequency_Hz
    r = PchipInterpolator(t, radius)(q)
    u = PchipInterpolator(t, speed)(q)
    volume = 4 * np.pi * r**3 / 3
    dot_volume = 4 * np.pi * r * r * u
    # Upstream pressure is -pa sin(wt): P_hat=i*pa, not a real cosine.
    raw = complex_harmonics(q, volume, case.frequency_Hz, harmonic_count)
    h = np.arange(1, harmonic_count + 1)
    normalized = raw / (1j**h)
    periodic_error = float(np.max(np.abs(r[-n:] - r[-2 * n : -n])) / case.radius_m)
    derivative_energy = float(np.mean(dot_volume**2))
    harmonic_energy = float(0.5 * np.sum(np.abs(2 * np.pi * case.frequency_Hz * h * raw) ** 2))
    unresolved = max(0.0, 1 - harmonic_energy / max(derivative_energy, 1e-300))
    metrics = dict(
        V_hat_real_m3=float(normalized[0].real),
        V_hat_imag_m3=float(normalized[0].imag),
        phase_rad=float(np.angle(normalized[0])),
        mean_volume_m3=float(volume.mean()),
        mean_dot_volume_squared_m6_s2=derivative_energy,
        mean_Rdot_squared_m2_s2=float(np.mean(u * u)),
        max_R_over_R0=float(r.max() / case.radius_m),
        min_R_over_R0=float(r.min() / case.radius_m),
        maximum_wall_mach=float(np.max(np.abs(u)) / case.sound_speed_m_s),
        periodic_error_R0=periodic_error,
        unresolved_dot_volume_energy_fraction=unresolved,
        periodic=periodic_error < 5e-4,
        harmonic_energy_resolved=unresolved < 0.01,
        weakly_compressible=float(np.max(np.abs(u)) / case.sound_speed_m_s) < 0.05,
    )
    metrics["valid"] = metrics["periodic"] and metrics["weakly_compressible"]
    return (
        metrics,
        normalized,
        dict(time_s=q, radius_m=r, radial_velocity_m_s=u, volume_m3=volume, dot_volume_m3_s=dot_volume),
    )


class ApecssBubbleResponseModel(BubbleResponseModel):
    """Offline only: production trajectories use ResponseTableBubbleModel."""

    def __init__(self, executable: str | Path):
        self.executable = Path(executable).resolve()
        if not self.executable.is_file():
            raise FileNotFoundError(f"Build APECSS first: {self.executable}")

    def response(self, R0_m, frequency_Hz, local_complex_pressure_Pa, medium_properties):
        raise RuntimeError(
            "APECSS is offline only; use solve() then ResponseTableBubbleModel in trajectories"
        )

    def solve(self, case: RadialCase, directory: str | Path):
        out = Path(directory).resolve()
        out.mkdir(parents=True, exist_ok=False)
        config = f"""BUBBLE
InitialRadius {case.radius_m:.17g}
RPModel RPAR
PressureAmbient {case.ambient_pressure_Pa:.17g}
END
GAS
EoS IG
ReferencePressure {case.ambient_pressure_Pa:.17g}
PolytropicExponent {case.polytropic_exponent:.17g}
END
LIQUID
ReferenceDensity {case.density_kg_m3:.17g}
ReferenceSoundSpeed {case.sound_speed_m_s:.17g}
Viscosity {case.viscosity_Pa_s:.17g}
END
INTERFACE
SurfaceTensionCoeff {case.sigma_water_N_m:.17g}
LipidCoatingModel Marmottant
SigmaInit {case.sigma_initial_N_m:.17g}
Elasticity {case.elasticity_N_m:.17g}
DilatationalViscosity {case.shell_viscosity_kg_s:.17g}
END
ODESOLVER
tolerance {case.tolerance:.17g}
minTimeStep 1e-15
maxTimeStep {1 / (case.frequency_Hz * case.samples_per_cycle):.17g}
END
RESULTS
Bubble
OutputFreqRP 1
OutputDigits 16
END
"""
        (out / "config.apecss").write_text(config)
        cmd = [
            str(self.executable),
            "-options",
            "config.apecss",
            "-freq",
            str(case.frequency_Hz),
            "-amp",
            str(case.pressure_Pa),
            "-tend",
            str(case.cycles / case.frequency_Hz),
        ]
        start = time.perf_counter()
        metadata = dict(
            case=asdict(case),
            command=cmd,
            apecss_commit=APECSS_COMMIT,
            executable_sha256=hashlib.sha256(self.executable.read_bytes()).hexdigest(),
            physics_status="partially validated",
            model="APECSS RPAR + reversible Marmottant",
            pressure_convention="p_inf=p0-pa*sin(wt); exported V_h divided by i**h",
        )
        try:
            process = subprocess.run(cmd, cwd=out, capture_output=True, text=True, timeout=120, check=False)
            (out / "stdout.log").write_text(process.stdout)
            (out / "stderr.log").write_text(process.stderr)
            metadata["returncode"] = process.returncode
            if process.returncode != 0:
                raise RuntimeError(f"APECSS returned {process.returncode}: see {out}")
            paths = list(out.glob("RPAR_*.txt"))
            if len(paths) != 1:
                raise RuntimeError("Expected exactly one APECSS RPAR trajectory")
            raw = np.loadtxt(paths[0])
            metrics, harmonics, cycles = analyze_cycles(raw, case)
            metadata["metrics"] = metrics
            np.savez_compressed(out / "cycles.npz", harmonics_cosine_m3=harmonics, **cycles)
            metadata["status"] = "COMPLETE"
            return metadata, harmonics, cycles
        except Exception as exc:
            metadata.update(status="FAILED", error=repr(exc))
            raise
        finally:
            metadata["elapsed_s"] = time.perf_counter() - start
            (out / "RESULT.json").write_text(json.dumps(metadata, indent=2) + "\n")
