"""Three SonoVue sensitivity cases, matched radius/force; not a production dependency."""

from dataclasses import asdict
import gzip
import hashlib
from pathlib import Path
import time

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import PchipInterpolator

from mevo_mb.bubbles.adaptive_table import ResponseNodeCache
from mevo_mb.bubbles.apecss import RadialCase
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution
from validation.history_force_sensitivity.model import solve
from scripts.stage3.common import new_run, finish, atomic_json, watermark

EXE = "third_party/build/apecss/ultrasound_apecss"
PARAM = "configs/stage2/sonovue_parameters.json"


def main():
    parameters = SonoVueParameterDistribution(PARAM).parameter_set(1.78e-6)
    case_parameters = dict(
        radius_m=1.78e-6,
        frequency_Hz=2.5e6,
        elasticity_N_m=parameters.shell_elasticity_N_m.value,
        shell_viscosity_kg_s=parameters.shell_viscosity_kg_s.value,
        sigma_initial_N_m=0.036,
    )
    config = dict(
        cases=["H1", "H2", "H3"],
        pressures_Pa=[5000, 20000, 150000],
        short_cycles=20,
        long_s=0.0001,
        samples_per_cycle=[64, 128, 256],
        convergence_endpoint_relative_limit=0.02,
        forcing="local standing-wave gradient: dp/dx=-k P sin(wt); fixed-position prescribed field",
        translation_viscosity_Pa_s=0.001,
        radial_effective_viscosity_Pa_s=0.002,
        production_averaging_control="Aprod=<F>steady_cycle/(6pi mu R0); Aresolved=F(t)/(6pi mu R(t)); B adds conservative added mass; C adds history to B",
    )
    out = new_run(
        "history",
        config,
        [
            "10.1063/1.3227903",
            "10.1103/p475-xggm; full equation source arXiv:2506.16413v2",
            "10.1121/1.3242346",
        ],
    )
    cache = ResponseNodeCache("results/stage2_5/response_table/node_cache", EXE)
    rows = []
    curves = {}
    convergence = []
    for label, pressure in zip(config["cases"], config["pressures_Pa"]):
        cp = dict(
            case_parameters,
            pressure_Pa=pressure,
            provenance="Tu2009 Fig4 parameter manifold, sigma_init=0.036 assumption",
        )
        if label == "H3":
            cp.update(
                elasticity_N_m=0.30,
                shell_viscosity_kg_s=3.2e-9,
                provenance="Tu2009 Table2 SonoVue case; reversible Marmottant sigma_init=0.036 is an assumption, not fitted here",
            )
        case = RadialCase(**cp)
        node = cache.get(case)
        if node["status"] != "COMPLETE" or not node["metrics"]["valid"]:
            raise RuntimeError("Invalid APECSS response " + label)
        folder = cache.root / node["key"] / Path(node["payload"]).parent
        raw_path = next(folder.glob("RPAR_*.txt.gz"))
        with gzip.open(raw_path, "rt") as stream:
            raw = np.loadtxt(stream)
        radial = PchipInterpolator(raw[:, 1], raw[:, 3])
        # Keep the full startup, then repeat only the final verified periodic cycle.
        period = 1 / case.frequency_Hz
        end = case.cycles * period
        previous = None
        for spc in config["samples_per_cycle"]:
            start = time.perf_counter()
            n = round(config["long_s"] / period) * spc
            t = np.arange(n + 1) * period / spc
            query = np.where(t <= end, t, (case.cycles - 1 + np.mod(t / period, 1)) * period)
            r = radial(query)
            force = (
                4 * np.pi * r**3 / 3 * (2 * np.pi / 1500 / period) * pressure * np.sin(2 * np.pi * t / period)
            )
            force[0] = 0.0
            v_res = force / (6 * np.pi * 0.001 * r)
            ares = dict(
                velocity_m_s=v_res,
                displacement_m=cumulative_trapezoid(v_res, t, initial=0),
                quasisteady_drag_N=-force,
                added_mass_N=np.zeros(len(t)),
                history_force_N=np.zeros(len(t)),
            )
            # Production assumes periodic radial response from t=0. Use its steady-cycle force;
            # Aresolved/B/C alone isolate inertia and history with matched startup forcing.
            cycle_force = np.array(
                [np.trapezoid(force[j : j + spc + 1], dx=period / spc) / period for j in range(0, n, spc)]
            )
            production_v = np.full(n, cycle_force[-1] / (6 * np.pi * 0.001 * case.radius_m))
            production_x = np.r_[0, np.cumsum(production_v) * period / spc]
            production_v = np.r_[production_v, production_v[-1]]
            models = {
                "Aprod": dict(
                    velocity_m_s=production_v,
                    displacement_m=production_x,
                    force_N=np.full(len(t), cycle_force[-1]),
                    quasisteady_drag_N=-6 * np.pi * 0.001 * case.radius_m * production_v,
                    added_mass_N=np.zeros(len(t)),
                    history_force_N=np.zeros(len(t)),
                ),
                "Aresolved": ares,
                "B": solve(t, r, force, history=False),
                "C": solve(t, r, force, history=True),
            }
            endpoints = {}
            for window, idx in [("short", 20 * spc), ("long", n)]:
                endpoint = {k: float(v["displacement_m"][idx]) for k, v in models.items()}
                endpoints[window] = endpoint
                if previous is not None:
                    convergence.append(
                        dict(
                            case=label,
                            window=window,
                            samples_per_cycle=spc,
                            changes={
                                k: abs(v - previous[window][k]) / max(abs(v), 1e-15)
                                for k, v in endpoint.items()
                            },
                        )
                    )
                if spc == config["samples_per_cycle"][-1]:
                    comparison = {}
                    for k, model in models.items():
                        sl = slice(0, idx + 1)
                        v = model["velocity_m_s"][sl]
                        re = 2 * 1000 * (case.radius_m if k == "Aprod" else r[sl]) * abs(v) / 0.001
                        target = 0.5 * endpoint["C"]
                        hits = np.flatnonzero(np.sign(target) * model["displacement_m"][sl] >= abs(target))
                        comparison[k] = dict(
                            displacement_m=endpoint[k],
                            mean_drift_velocity_m_s=endpoint[k] / t[idx],
                            peak_speed_m_s=float(max(abs(v))),
                            Re_max=float(max(re)),
                            Re_mean=float(np.trapezoid(re, t[sl]) / t[idx]),
                            selected_displacement_m=target,
                            time_to_selected_displacement_s=float(t[hits[0]]) if len(hits) else None,
                            force_budget={
                                name: dict(
                                    mean_N=float(np.trapezoid(values[sl], t[sl]) / t[idx]),
                                    rms_N=float(np.sqrt(np.trapezoid(values[sl] ** 2, t[sl]) / t[idx])),
                                )
                                for name, values in dict(
                                    primary_acoustic_N=model.get("force_N", force),
                                    quasisteady_drag_N=model["quasisteady_drag_N"],
                                    added_mass_N=model["added_mass_N"],
                                    history_force_N=model["history_force_N"],
                                ).items()
                            },
                        )
                    err = abs(endpoint["Aprod"] - endpoint["C"]) / max(abs(endpoint["C"]), 1e-15)
                    category = (
                        "small effect in tested condition"
                        if err < 0.05
                        else "modest approximation error"
                        if err < 0.10
                        else "non-negligible approximation"
                        if err < 0.20
                        else "quasi-steady approximation requires explicit limitation"
                    )
                    rows.append(
                        dict(
                            case=label,
                            window=window,
                            duration_s=float(t[idx]),
                            radial_case=asdict(case),
                            cache_key=node["key"],
                            raw_sha256=hashlib.sha256(raw_path.read_bytes()).hexdigest(),
                            max_R_over_R0=float(max(r[: idx + 1]) / case.radius_m),
                            min_R_over_R0=float(min(r[: idx + 1]) / case.radius_m),
                            models=comparison,
                            total_production_vs_C_relative_error=err,
                            mean_drift_relative_error=err,
                            averaging_Aprod_vs_Aresolved_relative_error=abs(
                                endpoint["Aprod"] - endpoint["Aresolved"]
                            )
                            / max(abs(endpoint["Aresolved"]), 1e-15),
                            added_mass_Aresolved_vs_B_relative_error=abs(
                                endpoint["Aresolved"] - endpoint["B"]
                            )
                            / max(abs(endpoint["B"]), 1e-15),
                            history_B_vs_C_relative_error=abs(endpoint["B"] - endpoint["C"])
                            / max(abs(endpoint["C"]), 1e-15),
                            classification=category,
                        )
                    )
            previous = endpoints
            if spc == config["samples_per_cycle"][-1]:
                with h5py.File(out / (label + ".h5"), "x") as h:
                    h["time_s"] = t
                    h["radius_m"] = r
                    h["instantaneous_force_N"] = force
                    for k, model in models.items():
                        for key, val in model.items():
                            h[f"{k}/{key}"] = val
                curves[label] = (t, r, force, models)
            atomic_json(
                out / "PROGRESS.json",
                dict(case=label, samples_per_cycle=spc, elapsed_s=time.perf_counter() - start),
            )
    final_conv = [v for v in convergence if v["samples_per_cycle"] == config["samples_per_cycle"][-1]]
    passed = all(
        max(v["changes"].values()) < config["convergence_endpoint_relative_limit"] for v in final_conv
    )
    summary = dict(
        status="PASS" if passed else "FAIL_CONVERGENCE",
        rows=rows,
        convergence=convergence,
        comparison_scope="Prescribed fixed-position standing field and repeated radial cycles; no finite-Re corrections, gas inertia, dissolution or pulse-off tail. Aresolved/B/C share exact R(t), F(t), t=0 rest, step and medium. Aprod uses the steady APECSS cycle-average force and R0, as production does. It is algebraic, with no independent velocity initial condition; its averaging/startup difference is reported separately, never called a history effect.",
        experimental_reproduction=False,
        production_history_dependency=False,
    )
    atomic_json(out / "SUMMARY.json", summary)
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.8), layout="constrained")
    for ax, (label, (t, r, force, models)) in zip(axs, curves.items()):
        for k in ("Aprod", "B", "C"):
            ax.plot(t * 1e6, models[k]["displacement_m"] * 1e6, label=k)
        ax.set(title=label, xlabel="Time (us)", ylabel="Displacement (um)")
        ax.legend()
    watermark(fig, "REFERENCE SENSITIVITY; NOT EXPERIMENTAL CALIBRATION")
    fig.savefig(out / "HF1_displacement.png", dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for offset, window in [(-0.2, "short"), (0.2, "long")]:
        ax.bar(
            np.arange(3) + offset,
            [100 * v["total_production_vs_C_relative_error"] for v in rows if v["window"] == window],
            0.4,
            label=window,
        )
    ax.set(xticks=np.arange(3), xticklabels=config["cases"], ylabel="Production vs C displacement error (%)")
    ax.legend()
    watermark(fig)
    fig.savefig(out / "HF2_error.png", dpi=220)
    plt.close(fig)
    t, r, f, models = curves["H2"]
    sel = t <= 4e-6
    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    for name, values in dict(
        primary=f,
        drag=models["C"]["quasisteady_drag_N"],
        added_mass=models["C"]["added_mass_N"],
        history=models["C"]["history_force_N"],
    ).items():
        ax.plot(t[sel] * 1e6, values[sel] * 1e9, label=name)
    ax.set(xlabel="Time (us)", ylabel="Force (nN)", title="H2 matched instantaneous force budget")
    ax.legend()
    watermark(fig)
    fig.savefig(out / "HF3_force_budget.png", dpi=220)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for window in ["short", "long"]:
        rr = [v for v in rows if v["window"] == window]
        ax.scatter(
            [max(v["max_R_over_R0"] - 1, 1 - v["min_R_over_R0"]) for v in rr],
            [100 * v["total_production_vs_C_relative_error"] for v in rr],
            label=window,
        )
    ax.set(xlabel="Maximum absolute radial excursion / R0", ylabel="Production vs C displacement error (%)")
    ax.legend()
    watermark(fig)
    fig.savefig(out / "HF4_excursion.png", dpi=220)
    plt.close(fig)
    finish(out, summary)


if __name__ == "__main__":
    main()
