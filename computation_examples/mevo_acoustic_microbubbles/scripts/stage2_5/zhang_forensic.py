"""Predeclared Zhang matrix; independent SI residual vs Stage2 acceleration matrix."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
import csv
import json
from pathlib import Path
import sys
import time

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from references.zhang2011_minimal import Parameters, solve, acceleration, residual
from mevo_mb.validation.zhang2011 import Zhang2011Parameters, Zhang2011EncapsulatedPairReference
from common import resumable_run, atomic_json, finish, watermark


def production_parameters(p):
    return Zhang2011Parameters(
        radius1_m=p.R10,
        radius2_m=p.R20,
        distance_m=p.D0,
        frequency_Hz=p.f,
        pressure_Pa=p.Pa,
        density_kg_m3=p.rho,
        sound_speed_m_s=p.c,
        ambient_pressure_Pa=p.P0,
        surface_tension_N_m=p.sigma,
        gas_exponent=p.gamma,
        viscosity_Pa_s=p.mu,
        elasticity_N_m=p.chi,
        shell_viscosity_Pa_s=p.mus,
        shell_thickness_m=p.epsilon,
    )


def coefficient(t, y, p, window):
    use = (t * p.f >= window[0] - 1e-7) & (t * p.f <= window[1] + 1e-7)
    dotv = 4 * np.pi * y[:, [0, 4]] ** 2 * y[:, [1, 5]]
    return float(
        p.rho / (4 * np.pi) * np.trapezoid(np.prod(dotv[use], axis=1), t[use]) / (t[use][-1] - t[use][0])
    )


def worker(job):
    name, p, expected, cycles, out = job
    out = Path(out) / name
    out.mkdir(exist_ok=True)
    result_path = out / "RESULT.json"
    if result_path.exists():
        return json.loads(result_path.read_text())
    start = time.perf_counter()
    a = solve(p, cycles=cycles)
    b = Zhang2011EncapsulatedPairReference(production_parameters(p)).solve(
        cycles=cycles, samples_per_cycle=128, rtol=2e-10, max_steps_per_cycle=96
    )
    by = np.column_stack(
        [
            b["radius_m"][:, 0],
            b["radial_speed_m_s"][:, 0],
            b["position_m"][:, 0],
            b["velocity_m_s"][:, 0],
            b["radius_m"][:, 1],
            b["radial_speed_m_s"][:, 1],
            b["position_m"][:, 1],
            b["velocity_m_s"][:, 1],
        ]
    )
    rows = []
    for impl, y in [("independent_residual", a["state"]), ("stage2_matrix", by)]:
        value = coefficient(a["time_s"], y, p, (cycles - 8, cycles))
        rows.append(
            dict(
                case=name,
                implementation=impl,
                expected_sign=expected,
                calculated_sign=int(np.sign(value)),
                calculated_N_m2=value,
                paper_N_m2=6.6e-18 if name == "C_fig3_failed" else None,
                relative_difference=abs(value - 6.6e-18) / 6.6e-18 if name == "C_fig3_failed" else None,
                averaging_window_cycles=[cycles - 8, cycles],
                rtol=2e-10,
                result="UNSPECIFIED_IN_PAPER"
                if expected is None
                else ("PASS_SIGN_ONLY" if int(np.sign(value)) == expected else "FAIL_SIGN"),
            )
        )
    window_coefficients = {}
    for w in [(0, 1), (0, 8), (8, 16), (24, 32), (56, 64), (cycles - 8, cycles)]:
        if w[1] <= cycles:
            window_coefficients[str(w)] = coefficient(a["time_s"], a["state"], p, w)
    np.savez_compressed(out / "trajectories.npz", time_s=a["time_s"], independent=a["state"], stage2=by)
    comparison = dict(
        max_relative_radius_difference=float(
            np.max(np.abs(a["state"][:, [0, 4]] - by[:, [0, 4]]) / np.array([p.R10, p.R20]))
        ),
        max_position_difference_m=float(np.max(np.abs(a["state"][:, [2, 6]] - by[:, [2, 6]]))),
    )
    meta = dict(
        case=name,
        parameters=asdict(p),
        rows=rows,
        comparison=comparison,
        windows=window_coefficients,
        runtime_s=time.perf_counter() - start,
        nfev=a["nfev"],
        solver="DOP853",
        atol_SI=a["atol"],
    )
    atomic_json(result_path, meta)
    print(name, rows[0]["calculated_N_m2"], comparison, flush=True)
    return meta


def debug_artifact(out, p):
    z = np.load(out / "C_fig3_failed/trajectories.npz")
    t, y = z["time_s"], z["independent"]
    details = {}
    acc = []
    for ti, yi in zip(t, y):
        ai = acceleration(ti, yi, p)
        acc.append(ai)
        for k, v in residual(ti, yi, ai, p, True).items():
            details.setdefault(k, []).append(v)
    radius = y[:, [0, 4]]
    speed = y[:, [1, 5]]
    V = 4 * np.pi * radius**3 / 3
    dv = 4 * np.pi * radius**2 * speed
    product = np.prod(dv, axis=1)
    distance = y[:, 6] - y[:, 2]
    integrand = p.rho / (4 * np.pi) * product / distance**2
    coeff_inst = p.rho / (4 * np.pi) * product
    cumulative = cumulative_trapezoid(coeff_inst, t, initial=0) / np.where(t > 0, t, 1)
    cycles = int(round(t[-1] * p.f))
    n = 128
    cf = np.array(
        [
            np.trapezoid(coeff_inst[i * n : (i + 1) * n + 1], t[i * n : (i + 1) * n + 1]) * p.f
            for i in range(cycles)
        ]
    )
    with h5py.File(out / "ZHANG_FIG3_FAILED_CASE.h5", "w") as h:
        h.attrs.update(
            parameters_json=json.dumps(asdict(p)),
            units="SI",
            sign="positive attraction",
            plotted_quantity="f_B / 1e-18 N m^2; dimensional coefficient, not a dimensionless force",
            physics_status="UNRESOLVED",
            source="Zhang2011 Eqs3-11",
        )
        data = dict(
            t=t,
            R=radius,
            Rdot=speed,
            V=V,
            Vdot=dv,
            center_positions=y[:, [2, 6]],
            separation=distance,
            accelerations=acc,
            Vdot_product=product,
            force_integrand_N=integrand,
            coefficient_integrand_N_m2=coeff_inst,
            cycle_average_coefficient_N_m2=cf,
            accumulated_average_N_m2=cumulative,
            paper_scaled_cycle_coefficient=cf / 1e-18,
            **details,
        )
        for k, v in data.items():
            h[k] = v
    fig, ax = plt.subplots(4, 2, figsize=(12, 12), constrained_layout=True)
    sel = t * p.f >= cycles - 3
    tx = (t[sel] - t[sel][0]) * 1e6
    ax[0, 0].plot(tx, radius[sel, 0] / p.R10)
    ax[0, 0].set_ylabel("R1 / R10")
    ax[0, 1].plot(tx, radius[sel, 1] / p.R20)
    ax[0, 1].set_ylabel("R2 / R20")
    ax[1, 0].plot(tx, np.asarray(details["drive_Pa"])[sel] / 1e3)
    ax[1, 0].set_ylabel("Drive (kPa)")
    ax[1, 1].plot(radius[sel, 0] / p.R10 - 1, radius[sel, 1] / p.R20 - 1)
    ax[1, 1].set(xlabel="R1/R10 - 1", ylabel="R2/R20 - 1", title="Radial phase relation")
    ax[2, 0].plot(tx, product[sel])
    ax[2, 0].set_ylabel("Vdot1 Vdot2 (m6/s2)")
    ax[2, 1].plot(tx, integrand[sel])
    ax[2, 1].set_ylabel("Force integrand (N)")
    ax[3, 0].plot(np.arange(cycles) + 0.5, cf / 1e-18, label="Per cycle")
    ax[3, 0].plot(t * p.f, cumulative / 1e-18, label="From initial rest")
    ax[3, 0].axhline(6.6, color="k", ls=":", label="Fig3 approximate digitization")
    ax[3, 0].legend(fontsize=8)
    ax[3, 0].set(xlabel="Elapsed cycles", ylabel="fB / 1e-18 N m2")
    ax[3, 1].plot(t * 1e6, distance * 1e6)
    ax[3, 1].set(xlabel="Time (us)", ylabel="Separation (um)")
    for i, j in [(0, 0), (0, 1), (1, 0), (2, 0), (2, 1)]:
        ax[i, j].set_xlabel("Time in final three cycles (us)")
    fig.suptitle("Zhang Figure 3 forensic: R1=3 um, R2=6 um, D0=120 um")
    watermark(fig, "PARTIALLY VALIDATED — FIGURE 3 SIGN UNRESOLVED")
    fig.savefig(out / "ZHANG_FIG3_DEBUG.png", dpi=170)
    plt.close(fig)


def main():
    start = time.perf_counter()
    base = Parameters()
    cases = [
        ("A_attraction", replace(base, R10=4e-6, D0=50e-6), 1, 64),
        ("B_repulsion", replace(base, D0=50e-6), -1, 64),
        ("C_fig3_failed", base, 1, 256),
        ("D_nearby_distance", replace(base, D0=70e-6), 1, 64),
        ("E_nearby_frequency", replace(base, f=2.20e6), None, 64),
        ("F_resonance", replace(base, R10=2.9e-6, R20=1.5e-6, D0=40e-6), 1, 64),
        ("weak_drive_diagnostic", replace(base, Pa=1000), None, 64),
    ]
    out = resumable_run(
        "zhang_forensic",
        {"cases": [(n, asdict(p), s, c) for n, p, s, c in cases], "workers": 6},
        ["10.1088/1674-1056/20/11/114302 Eqs3-11, Figs3-5"],
    )
    with ProcessPoolExecutor(max_workers=6) as pool:
        results = list(pool.map(worker, [(n, p, s, c, str(out)) for n, p, s, c in cases]))
    rows = [r for z in results for r in z["rows"]]
    with (out / "ZHANG_REPRODUCTION_MATRIX.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    debug_artifact(out, base)
    finish(
        out,
        dict(
            status="UNRESOLVED",
            physics_status="partially validated",
            cases=results,
            runtime_s=time.perf_counter() - start,
            independent_implementation_agrees=all(
                z["comparison"]["max_relative_radius_difference"] < 1e-6 for z in results
            ),
            fig3_sign_pass=next(z for z in results if z["case"] == "C_fig3_failed")["rows"][0]["result"]
            == "PASS_SIGN_ONLY",
        ),
    )


if __name__ == "__main__":
    main()
