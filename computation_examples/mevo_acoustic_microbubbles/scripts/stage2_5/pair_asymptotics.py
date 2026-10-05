"""Matched uncoated CW limits, no substitution of lipid/solid-shell constitutive laws."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import csv
import json
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.validation.zhang2011 import Zhang2011Parameters, Zhang2011EncapsulatedPairReference
from mevo_mb.validation.doinikov_bouakaz_pair import (
    DoinikovBouakazParameters,
    DoinikovBouakazPairReferenceModel,
)
from common import resumable_run, atomic_json, finish, watermark


class ContinuousDoinikov(DoinikovBouakazPairReferenceModel):
    """Replace only the published Gaussian envelope with declared uniform CW.

    Same pressure algebra, mass system, damping and force; benchmark boundary
    condition, not a claim that the publication used continuous excitation.
    """

    def rhs(self, t, y):
        base = super().rhs(t, y)
        p = self.p
        duration = p.pulse_cycles / p.frequency_Hz
        env = np.exp(-0.5 * ((t - duration / 2) / (duration / p.gaussian_duration_sigmas)) ** 2)
        dp = p.pressure_Pa * (1 - env) * np.sin(2 * np.pi * p.frequency_Hz * t)
        correction = np.r_[4 * np.pi * y[:2] ** 2 * dp, 0.0, 0.0]
        base[4:] += np.linalg.solve(self.mass_matrix(y[:4]), correction)
        return base


def isolated_linear_volumes(p):
    r = np.array([p.radius1_m, p.radius2_m])
    w = 2 * np.pi * p.frequency_Hz
    pg = p.ambient_pressure_Pa + 2 * p.surface_tension_N_m / r
    stiffness = (3 * p.gas_exponent * pg - 2 * p.surface_tension_N_m / r) / (p.density_kg_m3 * r * r)
    damping = 4 * p.viscosity_Pa_s / (p.density_kg_m3 * r * r) + 3 * p.gas_exponent * pg / (
        p.density_kg_m3 * p.sound_speed_m_s * r
    )
    return -4 * np.pi * r * p.pressure_Pa / p.density_kg_m3 / (stiffness - w * w + 1j * w * damping)


def work(job):
    kind, ratio, f, pressure, c, out = job
    key = f"{kind}_q{ratio}_f{f:g}_p{pressure:g}_c{c:g}"
    start = time.perf_counter()
    path = Path(out) / key
    path.mkdir(exist_ok=True)
    if (path / "RESULT.json").exists():
        saved = json.loads((path / "RESULT.json").read_text())
        return saved["result"]
    dp = DoinikovBouakazParameters(
        radius1_m=1.5e-6,
        radius2_m=2e-6,
        distance_m=ratio * 2e-6,
        frequency_Hz=f,
        pressure_Pa=pressure,
        sound_speed_m_s=c,
        pulse_cycles=64,
    )
    if kind == "doinikov":
        model = ContinuousDoinikov(dp, multipole_order=6)
    else:
        zp = Zhang2011Parameters(
            radius1_m=dp.radius1_m,
            radius2_m=dp.radius2_m,
            distance_m=dp.distance_m,
            frequency_Hz=f,
            pressure_Pa=pressure,
            density_kg_m3=dp.density_kg_m3,
            sound_speed_m_s=c,
            ambient_pressure_Pa=dp.ambient_pressure_Pa,
            surface_tension_N_m=dp.surface_tension_N_m,
            gas_exponent=dp.gas_exponent,
            viscosity_Pa_s=dp.viscosity_Pa_s,
            elasticity_N_m=0,
            shell_viscosity_Pa_s=0,
            shell_thickness_m=0,
        )
        model = Zhang2011EncapsulatedPairReference(zp)
    z = model.solve(cycles=64, samples_per_cycle=64, rtol=2e-9, max_steps_per_cycle=48)
    t = z["time_s"]
    use = t * f >= 56 - 1e-7
    tt = t[use]
    r = z["radius_m"][use]
    u = z["radial_speed_m_s"][use]
    v = z["velocity_m_s"][use]
    x = z["position_m"][use]
    d = z["distance_m"][use]
    dv = 4 * np.pi * r * r * u
    coeff = float(dp.density_kg_m3 / (4 * np.pi) * np.trapezoid(dv[:, 0] * dv[:, 1], tt) / (tt[-1] - tt[0]))
    vol = isolated_linear_volumes(dp)
    classic = float(
        dp.density_kg_m3 / (8 * np.pi) * (2 * np.pi * f) ** 2 * np.real(vol[0] * vol[1].conjugate())
    )
    # Recover hydrodynamic interaction from isolated added-mass momentum and
    # actual external drag. This is independent of evaluating Eq11 on volumes.
    acc = np.array([model.rhs(ti, np.r_[ri, xi, ui, vi])[6:] for ti, ri, xi, ui, vi in zip(tt, r, x, u, v)])
    if kind == "zhang":
        incident = np.column_stack(
            [
                -(r[:, 1] ** 2) * u[:, 1] / d**2 + r[:, 1] ** 3 * v[:, 1] / d**3,
                r[:, 0] ** 2 * u[:, 0] / d**2 + r[:, 0] ** 3 * v[:, 0] / d**3,
            ]
        )
    else:
        incident = np.zeros_like(v)
    drag = -12 * np.pi * dp.viscosity_Pa_s * r * (v - incident)
    force = (
        2 * np.pi / 3 * dp.density_kg_m3 * r**3 * acc + 2 * np.pi * dp.density_kg_m3 * r * r * u * v - drag
    )
    hyd = float(np.trapezoid((force[:, 0] - force[:, 1]) / 2, tt) / (tt[-1] - tt[0]))
    classic_force = classic / dp.distance_m**2
    row = dict(
        model=kind,
        ratio=ratio,
        R_over_d=1 / ratio,
        frequency_Hz=f,
        pressure_Pa=pressure,
        sound_speed_m_s=c,
        radius1_m=dp.radius1_m,
        radius2_m=dp.radius2_m,
        separation_m=dp.distance_m,
        shell_model="uncoated",
        classical_coefficient_N_m2=classic,
        coefficient_N_m2=coeff,
        classical_force_N=classic_force,
        hydrodynamic_force_N=hyd,
        coefficient_relative_error=abs(coeff - classic) / max(abs(classic), 1e-40),
        hydrodynamic_relative_error=abs(hyd - classic_force) / max(abs(classic_force), 1e-40),
        calculated_sign=int(np.sign(coeff)),
        classical_sign=int(np.sign(classic)),
        kd=2 * np.pi * f * dp.distance_m / c,
        max_radial_fraction=float(np.max(np.abs(z["radius_m"] / [dp.radius1_m, dp.radius2_m] - 1))),
        relative_distance_change=float(np.max(np.abs(z["distance_m"] / dp.distance_m - 1))),
        runtime_s=time.perf_counter() - start,
        rtol=2e-9,
        averaging_cycles=[56, 64],
    )
    np.savez_compressed(path / "trajectory.npz", **{k: v for k, v in z.items() if isinstance(v, np.ndarray)})
    atomic_json(path / "RESULT.json", dict(parameters=asdict(dp), result=row))
    print(key, row["coefficient_relative_error"], row["hydrodynamic_relative_error"], flush=True)
    return row


def main():
    start = time.perf_counter()
    ratios = [10, 20, 40, 80, 160]
    freqs = [0.15e6, 0.3e6, 0.6e6, 2.15e6, 2.95e6]
    config = dict(
        ratios=ratios,
        frequencies_Hz=freqs,
        pressures_Pa=[100, 500],
        sound_speeds=[1482.0, 1e9],
        shell_model="uncoated",
        cycles=64,
        workers=6,
        force_definition="Eq11 plus independent added-mass momentum residual",
        limitations="c=1e9 is formal incompressible limit, not physical water. finite-c large kd is excluded from regime map.",
    )
    out = resumable_run(
        "pair_asymptotics", config, ["Zhang2011 Eqs3-11", "DoinikovBouakaz2015 Eqs8-13,39-43,67-68"]
    )
    jobs = [
        (k, q, f, p, c, str(out))
        for k in ["zhang", "doinikov"]
        for c in [1482.0, 1e9]
        for p in [100, 500]
        for f in freqs
        for q in ratios
    ]
    with ProcessPoolExecutor(max_workers=6) as pool:
        rows = list(pool.map(work, jobs))
    with (out / "COMPARISON.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    scaling = []
    for k in ["zhang", "doinikov"]:
        for f in freqs:
            rr = [
                r
                for r in rows
                if r["model"] == k
                and r["frequency_Hz"] == f
                and r["pressure_Pa"] == 100
                and r["sound_speed_m_s"] == 1e9
            ]
            power = float(
                np.polyfit(
                    np.log([r["R_over_d"] for r in rr]),
                    np.log([r["coefficient_relative_error"] for r in rr]),
                    1,
                )[0]
            )
            scaling.append(dict(model=k, frequency_Hz=f, observed_error_power=power))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True)
    for k in ["zhang", "doinikov"]:
        rr = [
            r
            for r in rows
            if r["model"] == k
            and r["frequency_Hz"] == 0.3e6
            and r["pressure_Pa"] == 100
            and r["sound_speed_m_s"] == 1e9
        ]
        axes[0].loglog(
            [r["ratio"] for r in rr],
            [abs(r["hydrodynamic_force_N"] / rr[0]["classical_force_N"]) for r in rr],
            "o-",
            label=k,
        )
        axes[1].loglog(
            [r["R_over_d"] for r in rr],
            [r["coefficient_relative_error"] for r in rr],
            "o-",
            label=k + " Eq11",
        )
        axes[1].loglog(
            [r["R_over_d"] for r in rr],
            [r["hydrodynamic_relative_error"] for r in rr],
            "s--",
            label=k + " momentum",
        )
    axes[0].loglog(ratios, (np.array(ratios) / 10.0) ** -2, ":", label="Classical d^-2")
    axes[0].set(
        xlabel="d/max(R0)", ylabel="|F| / |F_classical(d=10R)|", title="Formal incompressible, 0.3MHz,100Pa"
    )
    axes[1].set(xlabel="max(R0)/d", ylabel="Relative error vs classical")
    for ax in axes:
        ax.legend(fontsize=7)
        ax.grid(alpha=0.2)
    watermark(fig)
    fig.savefig(out / "FORCE_AND_ERROR_SCALING.png", dpi=180)
    plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    for ax, k in zip(axes, ["classical", "zhang", "doinikov"]):
        rr = [
            r
            for r in rows
            if r["model"] == ("zhang" if k == "classical" else k)
            and r["pressure_Pa"] == 100
            and r["sound_speed_m_s"] == 1482.0
        ]
        data = np.array(
            [
                [
                    (
                        next(r for r in rr if r["frequency_Hz"] == f and r["ratio"] == q)[
                            "classical_sign" if k == "classical" else "calculated_sign"
                        ]
                    )
                    for q in ratios
                ]
                for f in freqs
            ]
        )
        ax.imshow(data, origin="lower", cmap="coolwarm", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(5), ratios)
        ax.set_yticks(range(5), [f / 1e6 for f in freqs])
        ax.set(xlabel="d/Rmax", ylabel="Frequency (MHz)", title=k)
    watermark(fig, "PARTIALLY VALIDATED — +1 ATTRACTION; LARGE kd IS FORMAL ONLY")
    fig.savefig(out / "SIGN_MAP.png", dpi=180)
    plt.close(fig)
    finish(
        out,
        dict(
            status="COMPARISON_COMPLETE",
            physics_status="partially validated",
            rows=rows,
            scaling=scaling,
            runtime_s=time.perf_counter() - start,
        ),
    )


if __name__ == "__main__":
    main()
