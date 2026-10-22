"""Stage 2F: primary-force identities and four-state published reference."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.bubbles.forces import primary_bjerknes
from mevo_mb.bubbles.translation import primary_cycle_average
from mevo_mb.validation.doinikov_dayton import DoinikovDaytonParameters, DoinikovDaytonSingleBubbleReference


def calculate(job):
    key, p, directory, fine = job
    out = Path(directory) / key
    out.mkdir()
    start = time.perf_counter()
    meta = {
        "parameters": asdict(p),
        "status": "FAILED",
        "source_doi": "10.1121/1.2215228",
        "equations": [6, 31, 32, 33, 36, 40],
        "physics_status": "partially validated",
        "pulse": "20 cycles, sin(omega*t-k*x), phase assumed; no postpulse drift included",
        "assumptions": "gamma=1.4 air; gas density=1.2 kg/m3; no fluid acoustic velocity, no history drag",
    }
    try:
        model = DoinikovDaytonSingleBubbleReference(p)
        z = model.solve(
            cycles=20,
            samples_per_cycle=512 if fine else 256,
            max_steps_per_cycle=256 if fine else 128,
            rtol=1e-10 if fine else 1e-9,
        )
        t, r = z["time_s"], z["outer_radius_m"]
        force = -4 * np.pi / 3 * r**3 * model.gradient(z["displacement_m"], t)
        n = (len(t) - 1) // 20
        cycle_force = np.array(
            [
                np.trapezoid(force[i * n : (i + 1) * n + 1], t[i * n : (i + 1) * n + 1]) * p.frequency_Hz
                for i in range(20)
            ]
        )
        fast_velocity = cycle_force / (6 * np.pi * p.viscosity_liquid_Pa_s * p.outer_radius_m)
        fast_x = np.r_[0.0, np.cumsum(fast_velocity) / p.frequency_Hz]
        late = t >= 10 / p.frequency_Hz
        late_force = float(np.trapezoid(force[late], t[late]) / (t[-1] - t[late][0]))
        analytic = model.linear_force()
        metrics = dict(
            displacement_m=float(z["displacement_m"][-1]),
            quasisteady_displacement_m=float(fast_x[-1]),
            approximation_error_m=float(fast_x[-1] - z["displacement_m"][-1]),
            late_mean_force_N=late_force,
            linear_force_equation63_N=analytic,
            relative_linear_force_difference=abs(late_force - analytic) / max(abs(analytic), 1e-30),
            max_Re=float(
                np.max(2 * p.density_liquid_kg_m3 * r * np.abs(z["velocity_m_s"]) / p.viscosity_liquid_Pa_s)
            ),
            max_radial_fraction=float(np.max(np.abs(r - p.outer_radius_m)) / p.outer_radius_m),
        )
        np.savez_compressed(
            out / "trajectory.npz",
            **z,
            force_N=force,
            cycle_force_N=cycle_force,
            cycle_time_s=np.arange(21) / p.frequency_Hz,
            quasisteady_position_m=fast_x,
            quasisteady_cycle_velocity_m_s=fast_velocity,
        )
        meta.update(status="COMPLETE", metrics=metrics)
    except Exception as e:
        meta["error"] = repr(e)
    meta["runtime_s"] = time.perf_counter() - start
    (out / "RESULT.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(key, meta["status"], flush=True)
    return key, meta


def main():
    root = Path.cwd()
    out = root / "results/stage2/translation" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True)
    jobs = []
    for r in (1.5, 2.0, 3.0):
        for pressure in (5.0, 40.0, 150.0):
            p = DoinikovDaytonParameters(outer_radius_m=r * 1e-6, pressure_Pa=pressure * 1e3)
            jobs.append((f"r{r:g}_p{pressure:g}", p, out, False))
    for thickness in (15.0, 200.0):
        for f in (1.0, 2.0, 3.0, 4.0, 6.0, 8.0):
            p = DoinikovDaytonParameters(shell_thickness_m=thickness * 1e-9, frequency_Hz=f * 1e6)
            jobs.append((f"s{thickness:g}_f{f:g}", p, out, False))
    jobs.append(("r2_p150_fine", DoinikovDaytonParameters(), out, True))
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=6, mp_context=multiprocessing.get_context("spawn")) as pool:
        results = dict(pool.map(calculate, jobs))
    coarse = np.load(out / "r2_p150/trajectory.npz")
    fine = np.load(out / "r2_p150_fine/trajectory.npz")
    time_error = float(np.max(np.abs(coarse["displacement_m"] - fine["displacement_m"][::2])))
    # All 60 coarse-table waveforms, using their own physical input phasor i*pa.
    radial_root = root / "results/stage2/single_bubble/20261001T201429472239Z"
    force_checks = []
    for path in sorted(radial_root.glob("node_*/RESULT.json")):
        meta = json.loads(path.read_text())
        z = np.load(path.parent / "cycles.npz")
        t, v = z["time_s"], z["volume_m3"]
        f = meta["case"]["frequency_Hz"]
        omega = 2 * np.pi * f
        # Append periodic endpoint for full trapezoidal cycle integration.
        tt = np.r_[t, t[-1] + (t[1] - t[0])]
        vv = np.r_[v, v[0]]
        G = np.array([8e7 + 2e7j, 1e7 - 4e7j, -3e7 + 5e7j])
        raw_G = 1j * G
        gradient = np.real(np.exp(1j * omega * tt[:, None]) * raw_G)
        averaged = primary_cycle_average(tt, vv, gradient)
        harmonic = primary_bjerknes(z["harmonics_cosine_m3"][0], G)
        error = float(np.linalg.norm(averaged - harmonic) / max(np.linalg.norm(harmonic), 1e-30))
        G2 = 0.25 * G
        gradient2 = gradient + np.real(np.exp(2j * omega * tt[:, None]) * (-G2))
        full2 = primary_cycle_average(tt, vv, gradient2)
        harmonic2 = harmonic + primary_bjerknes(z["harmonics_cosine_m3"][1], G2)
        force_checks.append(
            dict(
                case=path.parent.name,
                radial_periodic=meta["metrics"]["periodic"],
                monochromatic_relative_error=error,
                two_harmonic_relative_error=float(
                    np.linalg.norm(full2 - harmonic2) / max(np.linalg.norm(full2), 1e-30)
                ),
                fundamental_only_missing_harmonic_fraction=float(
                    np.linalg.norm(full2 - harmonic) / max(np.linalg.norm(full2), 1e-30)
                ),
            )
        )
    completed = all(r["status"] == "COMPLETE" for r in results.values())
    low_errors = [
        results[f"r{r:g}_p5"]["metrics"]["relative_linear_force_difference"] for r in (1.5, 2.0, 3.0)
    ]
    passed = (
        completed
        and time_error < 1e-10
        and max(low_errors) < 0.03
        and max(r["monochromatic_relative_error"] for r in force_checks) < 1e-8
    )
    summary = dict(
        status="PASS_REFERENCE_NUMERICS" if passed else "VALIDATION_FAILED",
        physics_status="partially validated",
        run=str(out.relative_to(root)),
        elapsed_s=time.perf_counter() - start,
        source_doi="10.1121/1.2215228",
        random_seed=None,
        git_commit=subprocess.check_output(
            ["git", "-c", f"safe.directory={root}", "rev-parse", "HEAD"], text=True
        ).strip(),
        reference_source_sha256=hashlib.sha256(
            (root / "src/mevo_mb/validation/doinikov_dayton.py").read_bytes()
        ).hexdigest(),
        hardware="results/stage2/system/SERVER_DOCTOR.json",
        workers=6,
        timestep_max_displacement_difference_m=time_error,
        small_signal_force_relative_errors=low_errors,
        primary_force_checks=force_checks,
        cases=results,
        limitations=[
            "No author reference trajectory supplied; Eq63 is an independent analytic cross-check",
            "Solid-shell reference is not SonoVue; pulse phase and air properties explicitly assumed",
            "Fast-model errors are measured, not fitted away; history force absent in both 2006 and fast model",
            "This benchmark does not establish quantitative validity of quasisteady transport in blood",
        ],
    )
    (out / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out.parent / "STATE.json").write_text(json.dumps(summary, indent=2) + "\n")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    z = coarse
    ax[0, 0].plot(z["time_s"] * 1e6, z["displacement_m"] * 1e6, label="Published coupled Eqs31–32")
    ax[0, 0].plot(
        z["cycle_time_s"] * 1e6, z["quasisteady_position_m"] * 1e6, label="Same forcing, cycle mean / Stokes"
    )
    ax[0, 1].plot(z["time_s"] * 1e6, z["velocity_m_s"], label="Coupled velocity")
    ax[0, 1].step(z["cycle_time_s"][1:] * 1e6, z["quasisteady_cycle_velocity_m_s"], label="Fast velocity")
    for r in (1.5, 2.0, 3.0):
        pressure = [5.0, 40.0, 150.0]
        ax[1, 0].plot(
            pressure,
            [results[f"r{r:g}_p{p:g}"]["metrics"]["displacement_m"] * 1e6 for p in pressure],
            "o-",
            label=f"R0={r:g} um",
        )
    for s in (15.0, 200.0):
        f = [1.0, 2.0, 3.0, 4.0, 6.0, 8.0]
        ax[1, 1].plot(
            f,
            [results[f"s{s:g}_f{v:g}"]["metrics"]["displacement_m"] * 1e6 for v in f],
            "o-",
            label=f"Shell={s:g} nm",
        )
    for a, xl, yl in zip(
        ax.flat,
        ["Time (us)", "Time (us)", "Pressure (kPa)", "Frequency (MHz)"],
        ["Displacement (um)", "Velocity (m/s)", "Displacement (um)", "Displacement (um)"],
    ):
        a.set(xlabel=xl, ylabel=yl)
        a.legend(fontsize=8)
        a.grid(alpha=0.25)
    fig.suptitle("Figure B — Doinikov–Dayton solid-shell reference; 20-cycle pulse, no history force")
    fig.savefig(out / "Figure_B_translation.png", dpi=180)
    plt.close(fig)
    print(
        json.dumps({k: v for k, v in summary.items() if k not in ("cases", "primary_force_checks")}, indent=2)
    )


if __name__ == "__main__":
    main()
