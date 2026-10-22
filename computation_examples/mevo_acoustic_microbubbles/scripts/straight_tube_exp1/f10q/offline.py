"""B1 exact saved common times; B2 bounded cached-field proposal design."""

import time
import itertools
import numpy as np
from scipy.optimize import minimize
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from scripts.straight_tube_exp1.f10o.analyze import load, metrics, topology_record
from .common import ROOT, BASE, OLD, CONFIG, read, write, sha, identity, event


def comparison():
    start = time.perf_counter()
    paths = {
        "M0": OLD / "campaign_v1/runs/OPTIMIZED_REPLAY",
        "B2": ROOT / "results/straight_tube_exp1/f10o/validation_v1/runs/STRONG_DT1",
        "HARD": ROOT / "results/straight_tube_exp1/f10p/validation_v1/runs/STRONG_DT1",
    }
    raw = {k: load(v) for k, v in paths.items()}
    common = []
    for t in raw["M0"][0]:
        if t <= 0.000174 + 1e-13 and all(np.min(abs(a[0] - t)) < 1e-13 for a in raw.values()):
            common.append(t)
    ts = np.array(common)
    xs = {k: np.array([x[np.argmin(abs(t - v))] for v in ts]) for k, (t, x) in raw.items()}
    assert len(ts) == 2 and np.isclose(ts[-1], 0.0001)
    m = {k: metrics(ts, x) for k, x in xs.items()}
    contrasts = {}
    for a, b in itertools.combinations(paths, 2):
        contrasts[a + "_VS_" + b] = dict(
            COM_m=float(np.linalg.norm(m[a]["centroid_m"] - m[b]["centroid_m"], axis=1).max()),
            fixed_id_m=float(np.linalg.norm(xs[a] - xs[b], axis=2).max()),
            Rg_m=float(abs(m[a]["rg_m"] - m[b]["rg_m"]).max()),
            Q_relative=float(
                (
                    np.linalg.norm(m[a]["tensor_m2"] - m[b]["tensor_m2"], axis=(1, 2))
                    / np.linalg.norm(m[a]["tensor_m2"], axis=(1, 2))
                ).max()
            ),
            xy_difference=float(abs(m[a]["directional_xy"] - m[b]["directional_xy"]).max()),
        )
    # Controls and initial states are checked at full precision, before any model contrast.
    original = event("OPTIMIZED_REPLAY", 0.002)
    for k, path in paths.items():
        inp = read(path / "INPUT.json")
        e = inp["events"][0] if k == "M0" else inp["control"]
        for key in ["amplitude_Pa", "phase_rad", "enabled"]:
            assert np.array_equal(e[key], original[key]), (k, key)
        np.testing.assert_array_equal(xs[k][0], xs["M0"][0])
    within = {}
    refs = {
        "M0": [(OLD / "campaign_v1/runs/OPTIMIZED_HALF_DT"), (OLD / "campaign_v1/runs/OPTIMIZED_FINE_CACHE")],
        "B2": [
            ROOT / "results/straight_tube_exp1/f10o/validation_v1/runs" / n
            for n in ["STRONG_DT05", "STRONG_FINE"]
        ],
        "HARD": [
            ROOT / "results/straight_tube_exp1/f10p/diagnostic_v1/runs" / n
            for n in ["STRONG_PREFIX_DT05", "STRONG_PREFIX_FINE"]
        ],
    }
    for k, others in refs.items():
        within[k] = {}
        for label, path in zip(["dt", "cache"], others):
            t, x = load(path)
            yy = np.array([x[np.argmin(abs(t - v))] for v in ts])
            mm = metrics(ts, yy)
            within[k][label] = dict(
                COM_m=float(np.linalg.norm(m[k]["centroid_m"] - mm["centroid_m"], axis=1).max()),
                ID_m=float(np.linalg.norm(xs[k] - yy, axis=2).max()),
                Rg_m=float(abs(m[k]["rg_m"] - mm["rg_m"]).max()),
            )
    quality = {}
    for k, path in paths.items():
        topology, labels, _, _ = topology_record(ts, xs[k])
        quality[k] = dict(
            topology_at_saved_times=topology,
            historical_result=read(path / "RESULT.json"),
            interpretation="Historical model-specific quality; no hard-sphere complementarity defined for M0 or B2",
        )
    write(
        BASE / "analysis/THREE_MODEL_COMPARISON.json",
        dict(
            identity=identity(),
            input_sha256={
                str((v / f).relative_to(ROOT)): sha(v / f)
                for v in [*paths.values(), *[z for val in refs.values() for z in val]]
                for f in ["trajectory.h5", "INPUT.json"]
            },
            common_saved_time_s=ts,
            no_interpolation=True,
            contrasts=contrasts,
            within_model_numerical=within,
            models={
                k: dict(
                    COM_m=m[k]["centroid_m"],
                    Rg_m=m[k]["rg_m"],
                    Q_m2=m[k]["tensor_m2"],
                    xy=m[k]["directional_xy"],
                    quality=quality[k],
                )
                for k in paths
            },
            scope="Only identical saved times 0 and 0.1 ms, inside common legal 0–0.174 ms window. Different contact/mobility approximations. No 2ms or 30ms inference.",
            wall_s=time.perf_counter() - start,
        ),
    )
    np.savez_compressed(
        BASE / "analysis/THREE_MODEL_INPUT.npz", time_s=ts, **{k + "_positions_m": v for k, v in xs.items()}
    )
    print("common comparison", contrasts)


def candidates():
    start = time.perf_counter()
    p = read(CONFIG)
    x = np.load(OLD / "initial_states.npz")["aggregate_m"]
    basis = CachedWaterBasis(OLD / "cache/coarse.h5")
    pb, gb = basis.query_basis(x)
    phys = F2Physics(radius_m=p["radius_m"], frequency_Hz=p["frequency_Hz"])
    response = ResponseSection(OLD / "inputs/response.h5", phys)
    zeta = 6 * np.pi * 0.001 * p["radius_m"]
    weak = event("SHAPE_REPLAY", 0.0005)
    u0 = basis.vector(weak["amplitude_Pa"], weak["phase_rad"], weak["enabled"])
    rows = []
    for axis, target in p["XY_targets_m"].items():
        for cap in p["candidate_design"]["pressure_caps_Pa"]:
            seed = u0 * min(0.8 * cap / abs(u0 @ pb).max(), 1)
            q0 = np.r_[seed.real, seed.imag] / 1500
            evals = []
            best = [np.inf, None]

            def evaluate(q):
                u = 1500 * (q[:4] + 1j * q[4:])
                pressure = u @ pb
                grad = np.einsum("i,inj->nj", u, gb)
                ap = abs(pressure)
                try:
                    vol = response.query(pressure)
                    primary = -0.5 * np.real(vol[:, 0, None] * grad.conj())
                    disp = 0.002 * primary.mean(0) / zeta
                    value = float(np.sum((disp - np.array(target)) ** 2) / (2e-7) ** 2)
                    valid = bool(ap.max() <= cap * (1 + 1e-9) and abs(u).max() <= 1500 * (1 + 1e-9))
                except ValueError:
                    disp = np.full(3, np.nan)
                    value = 1e6
                    valid = False
                rec = dict(
                    q=q.copy(),
                    predicted_displacement_m=disp,
                    objective=value,
                    valid=valid,
                    pressure_min_max_Pa=[float(ap.min()), float(ap.max())],
                )
                evals.append(rec)
                if valid and value < best[0]:
                    best[:] = [value, rec]
                return value

            def cons(q):
                u = 1500 * (q[:4] + 1j * q[4:])
                return np.r_[1 - abs(u) / 1500, 1 - abs(u @ pb) / cap]

            result = minimize(
                evaluate,
                q0,
                method="SLSQP",
                bounds=[(-1, 1)] * 8,
                constraints=[dict(type="ineq", fun=cons)],
                options=dict(maxiter=60, ftol=1e-9, eps=1e-6),
            )
            assert best[1] is not None
            q = np.array(best[1]["q"])
            u = 1500 * (q[:4] + 1j * q[4:])
            control = dict(weak, amplitude_Pa=abs(u), phase_rad=np.angle(u), enabled=[True] * 4)
            rows.append(
                dict(
                    name=f"{axis}_P{cap}",
                    axis=axis,
                    target_m=target,
                    pressure_cap_Pa=cap,
                    control=control,
                    proposal=best[1],
                    optimizer_status=str(result.message),
                    objective_evaluations=len(evals),
                    evaluations=evals,
                )
            )
    write(
        BASE / "CANDIDATES.json",
        dict(
            identity=identity(),
            count=len(rows),
            candidates=rows,
            cache_sha256=sha(OLD / "cache/coarse.h5"),
            initial_sha256=sha(OLD / "initial_states.npz"),
            response_sha256=sha(OLD / "inputs/response.h5"),
            wall_s=time.perf_counter() - start,
            scope="Frozen-state primary COM predictor ONLY; all candidates require actual independent 32-particle trajectories. Max initial pressure is a proposal cap; actual APECSS range checked every physical step. No clipping.",
        ),
    )
    print(
        "candidates",
        [
            (v["name"], v["proposal"]["predicted_displacement_m"] * 1e6, v["proposal"]["pressure_min_max_Pa"])
            for v in rows
        ],
    )


if __name__ == "__main__":
    comparison()
    candidates()
