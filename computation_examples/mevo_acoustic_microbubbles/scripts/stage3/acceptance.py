"""Focused force, clustering and FlowField acceptance with immutable run evidence."""

import argparse
import copy
from dataclasses import replace
import itertools

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist, squareform

from mevo_mb.bubbles.adaptive_table import ResponseNodeCache, harmonics
from mevo_mb.bubbles.apecss import RadialCase
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes, primary_bjerknes
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.coupling.experiment import flow_from_config
from mevo_mb.flow.fields import FlowField, ZeroFlowField
from mevo_mb.io.artifacts import new_run as particle_run
from mevo_mb.utils.config import load_config
from scripts.stage3.common import new_run, finish, atomic_json, watermark


class UniformField:
    frequency_Hz = 2.5e6

    def __init__(self, pressure):
        self.pressure = pressure

    def interpolate_complex_pressure(self, x):
        return np.full(len(x), self.pressure, dtype=complex)

    def interpolate_complex_gradient(self, x):
        return np.zeros((len(x), 3), dtype=complex)


def secondary():
    out = new_run(
        "secondary",
        dict(
            weak_pressures_Pa=[100, 200, 400, 800],
            distance_slope_tolerance=1e-10,
            pressure_slope_tolerance=0.05,
        ),
        ["classical cycle-averaged monopole identity", "Tu2009 10.1121/1.3242346"],
    )
    model = ClassicalFarFieldBjerknes()
    r = np.full(2, 1.78e-6)
    x = np.array([[0.0, 0, 0], [40e-6, 0, 0]])
    v = np.array([1 + 2j, 1 + 2j]) * 1e-19
    f = model.forces(x, r, v, 2.5e6, 1000)
    opposite = model.forces(x, r, v * np.array([1, -1]), 2.5e6, 1000)
    distances = np.array([20, 30, 40, 60, 80]) * 1e-6
    forces = [
        float(model.forces(np.array([[0.0, 0, 0], [d, 0, 0]]), r, v, 2.5e6, 1000)[0, 0]) for d in distances
    ]
    ds = float(np.polyfit(np.log(distances), np.log(forces), 1)[0])
    cache = ResponseNodeCache(
        "results/stage2_5/response_table/node_cache", "third_party/build/apecss/ultrasound_apecss"
    )
    param = SonoVueParameterDistribution("configs/stage2/sonovue_parameters.json").parameter_set(r[0])
    case = RadialCase(
        r[0],
        2.5e6,
        100,
        param.shell_elasticity_N_m.value,
        param.shell_viscosity_kg_s.value,
        0.036,
        provenance="Tu2009 Fig4, fixed Stage2.5 manifold",
    )
    pressure_rows = []
    for pressure in [100, 200, 400, 800]:
        node = cache.get(replace(case, pressure_Pa=pressure))
        if node["status"] != "COMPLETE" or not node["metrics"]["valid"]:
            raise RuntimeError("Weak-field radial failure")
        h = harmonics(node)
        force = float(model.forces(x, r, np.stack([h, h]), case.frequency_Hz, 1000)[0, 0])
        pressure_rows.append(
            dict(
                pressure_Pa=pressure,
                force_N=force,
                cache_key=node["key"],
                max_R_over_R0=node["metrics"]["max_R_over_R0"],
            )
        )
    ps = float(
        np.polyfit(
            np.log([q["pressure_Pa"] for q in pressure_rows]),
            np.log([q["force_N"] for q in pressure_rows]),
            1,
        )[0]
    )
    zero = model.forces(x, r, np.array([v[0], 0]), 2.5e6, 1000)
    # Independent full-cycle primary quadrature, including arbitrary DC volume.
    phase = 2 * np.pi * np.arange(8192) / 8192
    grad = np.array([1 + 3j, -2 + 0.2j, 0.4 - 0.3j]) * 1e7
    vh = (1 - 2j) * 1e-18
    direct = -np.mean(
        (1e-17 + np.real(vh * np.exp(1j * phase)))[:, None] * np.real(np.exp(1j * phase[:, None]) * grad),
        axis=0,
    )
    algebra = primary_bjerknes(np.array([vh]), grad[None])[0]
    pe = float(np.linalg.norm(algebra - direct) / np.linalg.norm(direct))
    tests = dict(
        SB1=bool(f[0, 0] > 0 and f[1, 0] < 0),
        SB2=bool(opposite[0, 0] < 0),
        SB3=abs(ds + 2) < 1e-10,
        SB4=abs(ps - 2) < 0.05,
        SB5=bool(np.allclose(f[0], -f[1], rtol=1e-14, atol=0)),
        SB6=bool(np.max(abs(zero)) == 0),
    )
    finish(
        out,
        dict(
            status="PASS" if all(tests.values()) and pe < 1e-12 else "FAIL",
            tests=tests,
            distance_slope=ds,
            pressure_slope=ps,
            pressure_rows=pressure_rows,
            distances_m=distances.tolist(),
            distance_force_N=forces,
            primary_full_cycle_relative_error=pe,
            primary_status="PASS" if pe < 1e-12 else "FAIL",
            ClassicalSecondaryBjerknes="ACCEPTED_FOR_SIMULATOR" if all(tests.values()) else "FAIL",
            scope="Mathematical and basic physical implementation acceptance; not near-contact validation",
        ),
    )


def clustering():
    c = load_config("configs/stage3/cylinder.yaml")
    c["particles"].update(
        count=8,
        duration_s=0.003,
        dt_s=1e-5,
        primary_enabled=False,
        secondary_enabled=True,
        pair_backend="numpy",
    )
    x = np.array(list(itertools.product([-10e-6, 10e-6], repeat=3)))
    pop = BubblePopulation(x, np.full(8, 1.78e-6), np.arange(8))
    config = dict(
        pressures_Pa=[5000, 10000, 15000],
        early_window_s=0.0001,
        duration_s=0.003,
        dt_s=[1e-5, 5e-6],
        initial_cube_edge_m=20e-6,
        far_field_fraction_min=0.80,
        convergence_relative_limit=0.02,
        seed=c["seed"],
    )
    out = new_run(
        "clustering",
        config,
        ["Lazarus2017 10.1121/1.5010170: trend benchmark, not parameter-identical reproduction"],
    )
    rows = []
    trajectories = {}
    for pressure in config["pressures_Pa"]:
        for enabled, dt in [(False, 1e-5), (True, 1e-5), (True, 5e-6)]:
            cc = copy.deepcopy(c)
            cc["particles"].update(secondary_enabled=enabled, dt_s=dt)
            run = particle_run(out, f"P{pressure}_SB{int(enabled)}_dt{dt}", cc)
            tracker = ParticleTracker(cc, UniformField(pressure), ZeroFlowField())
            summary = tracker.run(pop, run)
            with h5py.File(run / "particles/trajectories.h5") as h:
                t = h["time_s"][:]
                positions = h["positions_m"][:]
            nn = np.array(
                [np.min(squareform(pdist(p)) + np.eye(len(p)) * 1e10, axis=1).mean() for p in positions]
            )
            pair = np.array([pdist(p).mean() for p in positions])
            early = t <= config["early_window_s"] + 1e-14
            checks = [v for v in tracker.monitor.history if v["time_s"] <= config["early_window_s"] + 1e-14]
            far = sum(v["valid_pair_count"] for v in checks) / sum(v["pair_count"] for v in checks)
            threshold = 0.95 * nn[0]
            hits = np.flatnonzero(nn <= threshold)
            row = dict(
                pressure_Pa=pressure,
                secondary_enabled=enabled,
                dt_s=dt,
                run=str(run),
                early_NN_slope_m_s=float(np.polyfit(t[early], nn[early], 1)[0]),
                early_pair_contraction_m=float(pair[0] - pair[np.flatnonzero(early)[-1]]),
                attractive_pair_fraction=1.0 if enabled else 0.0,
                attractive_fraction_definition="Identical monodisperse responses in uniform common phase; all 28 pair correlations positive when enabled",
                time_to_5percent_NN_contraction_s=float(t[hits[0]]) if len(hits) else None,
                early_classical_policy_fraction=far,
                final_NN_m=float(nn[-1]),
                NEAR_CONTACT_PHYSICS_APPROXIMATED=summary["NEAR_CONTACT_PHYSICS_APPROXIMATED"],
                final_pair_regime=summary["physics_validity"]["pair_regime"],
                timings=summary["timings"],
            )
            trajectories[(pressure, enabled, dt)] = (t, positions, nn)
            rows.append(row)
            atomic_json(out / "ROWS.json", rows)
    convergence = []
    for pressure in config["pressures_Pa"]:
        tc, xc, nc = trajectories[(pressure, True, 1e-5)]
        tf, xf, nf = trajectories[(pressure, True, 5e-6)]
        # Compare early-time endpoint; contact-topology differences do not decide acceptance.
        ic = np.argmin(abs(tc - config["early_window_s"]))
        iff = np.argmin(abs(tf - config["early_window_s"]))
        error = float(np.linalg.norm(xc[ic] - xf[iff]) / max(np.linalg.norm(xf[iff] - xf[0]), 1e-15))
        convergence.append(
            dict(
                pressure_Pa=pressure,
                early_displacement_relative_error=error,
                final_NN_relative_difference=float(abs(nc[-1] - nf[-1]) / nf[-1]),
            )
        )
    on = [v for v in rows if v["secondary_enabled"] and v["dt_s"] == 1e-5]
    off = [v for v in rows if not v["secondary_enabled"]]
    tests = dict(
        ON_OFF=all(
            a["early_NN_slope_m_s"] < -1e-8 and abs(b["early_NN_slope_m_s"]) < 1e-12 for a, b in zip(on, off)
        ),
        pressure_trend=all(b["early_NN_slope_m_s"] < a["early_NN_slope_m_s"] for a, b in zip(on, on[1:])),
        timestep=all(q["early_displacement_relative_error"] < 0.02 for q in convergence),
        early_validity=all(v["early_classical_policy_fraction"] >= 0.80 for v in on),
    )
    fig, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for p in config["pressures_Pa"]:
        t, x, nn = trajectories[(p, True, 1e-5)]
        ax.plot(t * 1e3, nn * 1e6, label=f"{p / 1000:g} kPa, ON")
    t, x, nn = trajectories[(5000, False, 1e-5)]
    ax.plot(t * 1e3, nn * 1e6, "k--", label="Secondary OFF")
    ax.set(xlabel="Time (ms)", ylabel="Mean nearest neighbour (um)")
    ax.legend()
    watermark(fig, "EARLY TREND ACCEPTANCE; LATE NEAR-FIELD APPROXIMATION")
    fig.savefig(out / "clustering_pressure_trend.png", dpi=240)
    plt.close(fig)
    finish(
        out,
        dict(
            status="PASS" if all(tests.values()) else "FAIL",
            tests=tests,
            rows=rows,
            convergence=convergence,
            benchmark="trend benchmark, not parameter-identical reproduction",
            late_interpretation="Geometric topology and engineering trend; near-contact is not validated pair physics",
        ),
    )


class AnalyticFlow(FlowField):
    def velocity(self, x, time_s):
        x = np.asarray(x)
        return np.column_stack([0.001 + 2 * x[:, 1] + 0.1 * time_s, -0.3 * x[:, 2], 0.2 * x[:, 1]])


def flow():
    import pyvista as pv

    out = new_run(
        "flow",
        dict(duration_s=0.001, dt_s=2e-5, absolute_trajectory_tolerance_m=1e-10),
        ["FlowField API; linear tetrahedral space/time interpolation identity"],
    )
    c = load_config("configs/stage3/cylinder.yaml")
    c["particles"].update(
        count=8, duration_s=0.001, dt_s=2e-5, primary_enabled=False, secondary_enabled=False
    )
    pop = BubblePopulation.initialize(c)
    points = np.array(list(itertools.product([-0.002, 0.002], repeat=3)))
    mesh = pv.PolyData(points).delaunay_3d()
    files = []
    analytic = AnalyticFlow()
    for idx, t in enumerate([0.0, 0.001]):
        mesh.point_data["Velocity"] = analytic.velocity(mesh.points, t)
        path = out / f"flow_{idx}.vtu"
        mesh.save(path)
        files.append(str(path))
    providers = {
        "analytic": analytic,
        "vtk": flow_from_config(dict(c, flow=dict(backend="vtk", files=files, times_s=[0.0, 0.001]))),
    }
    trajectories = {}
    timings = {}
    runs = {}
    for name, provider in providers.items():
        run = particle_run(out, name, c)
        summary = ParticleTracker(c, UniformField(0), provider).run(pop, run)
        with h5py.File(run / "particles/trajectories.h5") as h:
            trajectories[name] = h["positions_m"][:]
        timings[name] = summary["timings"]
        runs[name] = str(run)
    error = float(np.max(np.linalg.norm(trajectories["analytic"] - trajectories["vtk"], axis=-1)))
    finish(
        out,
        dict(
            status="PASS" if error < 1e-10 else "FAIL",
            maximum_trajectory_error_m=error,
            runs=runs,
            timings=timings,
            scope="Synthetic linear spatial and temporal velocity on tetrahedral VTU; does not validate a patient CFD solution",
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("task", choices=["secondary", "clustering", "flow"])
    a = p.parse_args()
    globals()[a.task]()
