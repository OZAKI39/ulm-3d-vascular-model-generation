"""Bounded GPU force screening; no integrated or fabricated trajectories."""

import argparse
import time
import numpy as np
from scipy.optimize import minimize
from mevo_mb.fast_sim.f8_predictor import MatchedPredictor
from mevo_mb.fast_sim.nearfield import pair_reference
from scripts.fast_pzt_f7.fixed import local_numpy, probe
from .common import ROOT, BASE, js, write, identity, objects


def stats(values):
    v = np.asarray(values)
    return dict(
        min=float(v.min()),
        mean=float(v.mean()),
        median=float(np.median(v)),
        max=float(v.max()),
        quantiles=np.quantile(v, [0, 0.05, 0.25, 0.5, 0.75, 0.95, 1]),
    )


def event(u):
    return dict(
        time_s=0.0,
        amplitude_Pa=abs(u),
        phase_rad=np.angle(u),
        enabled=[True] * 4,
        secondary_enabled=True,
        flow_enabled=True,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter()
    proto = js(ROOT / "configs/fast_pzt/f10a/protocol.json")
    write(out / "REGISTERED_PROTOCOL.json", proto)
    write(out / "IDENTITY.json", identity())
    radial = js(BASE / "radial_v1/RESULTS.json")
    obj = objects()
    p, b, r, g, flow = obj
    cfg = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    x0 = np.array(cfg["environment_Y"]["positions_m"])
    mature = np.array(cfg["environment_A"]["positions_m"])
    original = cfg["environment_Y"]["events"]["unchanged"]
    original_u = np.array(original["amplitude_Pa"]) * np.exp(1j * np.array(original["phase_rad"]))
    upper = b.origin + b.spacing * (np.array(b.shape) - 1)
    radius = p.radius_m

    def support(x):
        geo = g.query(x)
        cf = flow.sample(x)
        gap = geo["distance_m"] - radius
        field = np.all((x >= b.origin + radius) & (x <= upper - radius), axis=1)
        return dict(
            valid=bool(np.all(geo["inside"]) and np.all(gap >= 0) and np.all(cf["valid"]) and field.all()),
            gap=gap,
            inside=geo["inside"],
            field=field,
            flow=cf,
        )

    # Finite geometry-only candidates along the same mother segment; no force-based cherry-picking.
    fractions = np.linspace(0, 1, 41)
    shifts = np.c_[np.zeros(41), np.zeros(41), -x0.mean(0)[2] * fractions]
    geometry_rows = []
    for frac, shift in zip(fractions, shifts):
        s = support(x0 + shift)
        geometry_rows.append(
            dict(fraction=frac, translation_m=shift, valid=s["valid"], minimum_gap_m=s["gap"].min())
        )
    legal = [q for q in geometry_rows if q["valid"]]
    assert legal
    interior = max(legal, key=lambda q: q["minimum_gap_m"])
    middle = min(legal, key=lambda q: abs(q["fraction"] - interior["fraction"] / 2))
    states = {
        "near": x0,
        "middle": x0 + middle["translation_m"],
        "interior": x0 + interior["translation_m"],
        "mature": mature,
    }
    write(
        out / "GEOMETRY_SELECTION.json",
        dict(
            candidates=geometry_rows,
            chosen={"middle": middle, "interior": interior},
            field_lower_m=b.origin,
            field_upper_m=upper,
        ),
    )
    predictor = MatchedPredictor(b, r, p, mode="diagonal_stokes")
    rpy = MatchedPredictor(b, r, p, mode="rpy_free_space")
    rows = []
    state_rows = {}
    all_candidates = []
    rng = np.random.default_rng(proto["search"]["seed"])
    optim_total = 0.0
    for state, x in states.items():
        spatial = support(x)
        cf = spatial["flow"]
        uf = cf["velocity_m_s"] if state != "mature" else np.zeros_like(x)
        assert np.isfinite(uf).all()
        dx = x - x.mean(0)
        rg = np.sqrt(np.mean(np.sum(dx**2, axis=1)))
        dist = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
        np.fill_diagonal(dist, np.inf)
        predictor.prepare(x, uf)
        rpy.prepare(x, uf)
        samples = predictor.samples.numpy()
        samples = samples[..., 0] + 1j * samples[..., 1]
        gradients = predictor.gradients.numpy()
        gradients = gradients[..., 0] + 1j * gradients[..., 1]
        grad = cf["gradient_s_inv"]
        strain = 0.5 * (grad + grad.swapaxes(1, 2))
        residence = np.min(
            np.where(
                abs(uf.mean(0)) > 1e-15,
                np.where(uf.mean(0) > 0, upper - x.mean(0), x.mean(0) - b.origin)
                / np.maximum(abs(uf.mean(0)), 1e-30),
                np.inf,
            )
        )
        state_rows[state] = dict(
            positions_m=x,
            particle_ids=np.arange(32),
            translation_m=x.mean(0) - x0.mean(0) if state != "mature" else None,
            minimum_pair_distance_m=dist.min(),
            initial_max_overlap_m=max(0.0, 2 * radius - dist.min()),
            minimum_wall_gap_m=spatial["gap"].min(),
            wall_gap_m=spatial["gap"],
            support_valid=spatial["valid"],
            per_bubble_flow_m_s=uf,
            mean_flow_m_s=uf.mean(0),
            flow_range_m_s=np.ptp(uf, axis=0),
            flow_relative_rms_m_s=np.sqrt(np.mean(np.sum((uf - uf.mean(0)) ** 2, axis=1))),
            shear_rate_s_inv=np.sqrt(2 * np.sum(strain**2, axis=(1, 2))),
            estimated_stored_field_residence_s=residence,
            pressure_original_Pa=stats(abs(samples @ original_u)),
            source="F8 environment_A mature initial state (no flow/no wall)"
            if state == "mature"
            else "F9/F8 environment_Y same32 relative coordinates; geometry-only rigid translation",
            actual_CFD_gradient_s_inv=grad,
            CFD_tetra_id=cf["tetra_id"],
        )
        seeds = [np.r_[np.ones(4), np.angle(original_u[1:] / original_u[0])]] + [
            np.r_[rng.uniform(0.1, 1, 4), rng.uniform(-np.pi, np.pi, 3)]
            for _ in range(proto["search"]["random_candidates_per_state"])
        ]
        for pressure in proto["pressures_Pa"]:
            scaled = original_u * pressure / np.median(abs(samples @ original_u))
            budget = float(np.sum(abs(scaled) ** 2))
            maxP = float(abs(samples @ scaled).max())

            def control(z):
                v = z[:4] * np.exp(1j * np.r_[0, z[4:]])
                median = np.median(abs(samples @ v))
                if median <= 1e-15:
                    raise ValueError("zero candidate")
                return v * pressure / median

            def margins(u):
                pp = abs(samples @ u)
                return np.r_[
                    1 - abs(u) / proto["source_cap_Pa"],
                    1 - np.sum(abs(u) ** 2) / budget,
                    1.01 - pp.max() / maxP,
                    pp.min() / 0.1 - 1,
                    1 - pp.max() / 40000,
                ]

            best = [-np.inf, None]
            calls = 0
            invalid = 0
            cached = {}
            t0 = time.perf_counter()

            def evaluate(z):
                nonlocal calls, invalid
                key = np.asarray(z).tobytes()
                if key in cached:
                    return cached[key]
                if calls >= proto["search"]["max_optimizer_evaluations_per_state_pressure"] + len(seeds):
                    return 1e10
                calls += 1
                try:
                    u = control(z)
                    m = margins(u)
                    # All objective calls query actual harmonics/pair closure on the GPU.
                    v = predictor.evaluate(abs(u), np.angle(u))
                    score = float((v - uf).mean(0)[2])
                    legal = bool(m.min() >= -1e-8)
                    if legal and score > best[0]:
                        best[:] = [score, u.copy()]
                    cost = -score / 1e-5
                    all_candidates.append(
                        dict(
                            state=state,
                            pressure_Pa=pressure,
                            evaluation=calls,
                            source_amplitudes_Pa=abs(u),
                            source_phases_rad=np.angle(u),
                            constraint_min=m.min(),
                            feasible=legal,
                            mean_acoustic_Z_m_s=score,
                        )
                    )
                except ValueError as e:
                    invalid += 1
                    cost = 1e10
                    all_candidates.append(
                        dict(state=state, pressure_Pa=pressure, evaluation=calls, error=str(e))
                    )
                cached[key] = cost
                return cost

            costs = [evaluate(z) for z in seeds]
            feasible_seeds = [
                (cost, z) for cost, z in zip(costs, seeds) if margins(control(z)).min() >= -1e-8
            ]
            initial = min(feasible_seeds, key=lambda v: v[0])[1] if feasible_seeds else seeds[0]
            fit = minimize(
                evaluate,
                initial,
                method="SLSQP",
                bounds=[(0.001, 1)] * 4 + [(-np.pi, np.pi)] * 3,
                constraints=[dict(type="ineq", fun=lambda z: margins(control(z)))],
                options=dict(maxiter=proto["search"]["max_optimizer_iterations"], ftol=1e-10),
            )
            optim_s = time.perf_counter() - t0
            optim_total += optim_s
            controls = {"S": scaled, "O": best[1]}
            for mode, u in controls.items():
                if u is None:
                    rows.append(
                        dict(state=state, pressure_Pa=pressure, mode=mode, status="NO_FEASIBLE_CONTROL")
                    )
                    continue
                data = predictor.evaluate(abs(u), np.angle(u), details=True)
                pp = samples @ u
                gg = np.einsum("nka,a->nk", gradients, u)
                h = r.query(pp)
                fp = -0.5 * np.real(h[:, 0, None] * gg.conj())
                fs = pair_reference(x, h, pp, p, predictor.chi)
                np.testing.assert_allclose(data["primary_N"], fp, rtol=3e-10, atol=1e-23)
                np.testing.assert_allclose(data["secondary_N"], fs, rtol=3e-10, atol=1e-23)
                vel = local_numpy(x, data["velocity_m_s"], p)
                acoustic = vel - uf
                # The component balance makes the group M0 mean equal to its diagonal predictor; verify.
                np.testing.assert_allclose(vel.mean(0), data["velocity_m_s"].mean(0), rtol=1e-10, atol=1e-16)
                horizon = proto["static_horizon_s"]
                forecast = np.concatenate([x + t * vel for t in np.linspace(0, horizon, 9)])
                future = support(forecast)
                rate = np.mean(np.sum(dx * (vel - vel.mean(0)), axis=1)) / rg
                radial_ok = radial["pressure_admission"][str(pressure)] and all(
                    max(v.values()) <= 0.001 for v in radial["convergence"].values()
                )
                effect = acoustic.mean(0)[2] * horizon
                admitted = bool(
                    state != "mature"
                    and mode == "O"
                    and radial_ok
                    and spatial["valid"]
                    and future["valid"]
                    and effect >= proto["minimum_static_Z_increment_m"]
                    and np.max(abs(u)) <= proto["source_cap_Pa"]
                    and dist.min() >= 1.99 * radius
                )
                row = dict(
                    state=state,
                    pressure_Pa=pressure,
                    mode=mode,
                    status="VALID_STATIC",
                    event=event(u),
                    per_bubble_pressure_Pa=abs(pp),
                    pressure_statistics_Pa=stats(abs(pp)),
                    source_budget_Pa2=np.sum(abs(u) ** 2),
                    S_budget_Pa2=budget,
                    S_max_pressure_Pa=maxP,
                    primary_N=fp,
                    secondary_N=fs,
                    primary_mean_N=fp.mean(0),
                    primary_rms_N=np.sqrt(np.mean(np.sum(fp**2, axis=1))),
                    secondary_rms_N=np.sqrt(np.mean(np.sum(fs**2, axis=1))),
                    per_bubble_velocity_M0_m_s=vel,
                    mean_acoustic_M0_m_s=acoustic.mean(0),
                    mean_total_M0_m_s=vel.mean(0),
                    relative_velocity_rms_m_s=np.sqrt(np.mean(np.sum((vel - vel.mean(0)) ** 2, axis=1))),
                    Rg_rate_m_s=rate,
                    forecast_minimum_gap_m=future["gap"].min(),
                    forecast_support_valid=future["valid"],
                    predicted_Z_increment_10ms_m=effect,
                    high_harmonic_dotV_energy_fraction=np.sum((np.arange(2, 9) * abs(h[:, 1:])) ** 2)
                    / np.sum((np.arange(1, 9) * abs(h)) ** 2),
                    F2_near_stats=predictor.near.numpy(),
                    table_radial_gate=radial_ok,
                    DIRECTION_FEASIBLE=bool(effect > 0),
                    trajectory_admitted=admitted,
                    HIGH_DRIVE_EXPLORATORY=bool(np.max(abs(u)) > proto["historic_source_cap_Pa"]),
                    optimizer=dict(
                        evaluations=calls,
                        invalid=invalid,
                        wall_s=optim_s,
                        success=bool(fit.success),
                        message=str(fit.message),
                    )
                    if mode == "O"
                    else None,
                )
                if pressure == 3000 and state in ["near", "mature"]:
                    row["RPY_free_space_acoustic_mean_m_s"] = (rpy.evaluate(abs(u), np.angle(u)) - uf).mean(0)
                if pressure == 1000:
                    actual = probe(x, event(u), obj, "M0", state != "mature", state != "mature")
                    row["actual_M0_probe_max_error_m_s"] = np.max(abs(actual["velocity_m_s"] - vel))
                    np.testing.assert_allclose(actual["velocity_m_s"], vel, rtol=2e-8, atol=2e-15)
                rows.append(row)
            print(state, pressure, "O Z um/s", None if best[1] is None else best[0] * 1e6, flush=True)
    write(out / "CANDIDATES.json", all_candidates)
    write(
        out / "RESULTS.json",
        dict(
            states=state_rows,
            rows=rows,
            total_evaluations=len(all_candidates),
            optimizer_s=optim_total,
            process_wall_s=time.perf_counter() - begin,
            scope="Frozen-state GPU force queries, CPU reference M0 relative mobility checked against production; not integrated trajectories. RPY free space sensitivity only.",
            admitted=[
                dict(state=q["state"], pressure_Pa=q["pressure_Pa"])
                for q in rows
                if q.get("trajectory_admitted")
            ],
        ),
    )


if __name__ == "__main__":
    main()
