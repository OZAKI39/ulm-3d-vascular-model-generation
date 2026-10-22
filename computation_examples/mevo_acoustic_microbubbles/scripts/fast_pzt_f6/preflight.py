"""Freeze targets first, then screen only static drift before any navigation run."""

import json
import argparse
import copy
import h5py
import numpy as np
from scipy.optimize import minimize
from mevo_mb.fast_sim.consistent_field import query_numpy
from .common import ROOT, BASE, F5, clean, write, js, identity, objects


def field_samples(basis, x):
    vals = [query_numpy(v, basis.origin, basis.spacing, x) for v in basis.db.pressure_basis]
    return np.stack([v[0] for v in vals], axis=-1), np.stack([v[1] for v in vals], axis=-1)


def primary_from_control(response, p, P, G, amplitude, phase):
    control = np.asarray(amplitude) * np.exp(1j * np.asarray(phase))
    pressure = P @ control
    gradient = np.einsum("ndk,k->nd", G, control)
    V = response.query(pressure)[:, 0]
    force = -0.5 * np.real(V[:, None] * gradient.conj())
    return force / (6 * np.pi * p.viscosity_Pa_s * p.radius_m), pressure


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--checks", default="checks_v1")
    a = ap.parse_args()
    assert js(BASE / a.checks / "CHECKS.json")["status"] == "PASS"
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    p, basis, response, geometry, flow = objects()
    with h5py.File(ROOT / "results/fast_pzt/f2/mature_v1/unchanged/trajectory.h5") as h:
        mature = h["positions_m"][0, 0]
    with h5py.File(F5 / "science_v1/main/trajectory.h5") as h:
        switch = h["positions_m"][600, 0]
    cap = 291.237687721856
    old5 = js(F5 / "ACOUSTIC_RELEASE_PROTOCOL.json")["events"][0]
    proto = dict(
        schema="fast_pzt_f6_protocol_v1",
        identity=identity(),
        amplitude_cap_Pa=cap,
        cap_scope="Inherited campaign cap, not hardware calibration",
        dt_s=1e-5,
        sample_interval_s=0.0005,
        target_radius_m=2e-6,
        dwell_s=0.05,
        min_largest_fraction=0.9,
        min_resolvable_control_effect_m=0.5e-6,
        dt_centroid_error_limit_m=0.5e-6,
        dt_rg_relative_limit=0.05,
        dt_largest_fraction_limit=0.1,
        max_overlap_over_R=0.01,
        max_correction_over_R=0.1,
        evidence_level="EXPLORATORY",
        environments=[],
    )
    for name, x, offset, horizon in [
        ("zero_flow", mature, [0, -10e-6, 0], 1.0),
        ("slow_affine", mature, [0, -10e-6, 0], 1.0),
        ("Y_CFD", switch, [0, 0, 10e-6], 0.3),
    ]:
        center = x.mean(0)
        target = center + offset
        # Targets are fixed by this geometric rule, not by screening or trajectories.
        target_check = geometry.query(target[None])
        assert target_check["inside"].all() and target_check["distance_m"][0] > p.radius_m + 2e-6
        query_numpy(basis.db.pressure_basis[0], basis.origin, basis.spacing, np.vstack([x, target]))
        if name == "Y_CFD":
            assert flow.sample(np.vstack([x, target]))["valid"].all()
        cfg = dict(
            name=name,
            positions_m=x,
            target_m=target,
            initial_centroid_m=center,
            duration_s=horizon,
            initial_target_distance_m=np.linalg.norm(offset),
            initial_state_source="F5 exact 300 ms switch"
            if name == "Y_CFD"
            else "F2/F3 archived mature state",
            flow_enabled=name != "zero_flow",
            flow_velocity_m_s=[2e-6, 0, 0],
            flow_gradient_s_inv=[[0, 0.05, 0], [0, 0, 0], [0, 0, 0]],
            wall_enabled=name == "Y_CFD",
            split_group="water_2p5MHz_" + name + "_target_01",
            partition={"zero_flow": "train", "slow_affine": "validation", "Y_CFD": "test"}[name],
            events={},
        )
        proto["environments"].append(cfg)
    write(out / "TARGETS_BEFORE_SCREENING.json", proto)
    diagnostics = []
    for cfg in proto["environments"]:
        name = cfg["name"]
        center = cfg["initial_centroid_m"]
        target = cfg["target_m"]
        delta = target - center
        probes = np.vstack(
            [center, center + 0.5 * delta, target, target + np.eye(3) * 2e-6, target - np.eye(3) * 2e-6]
        )
        P, G = field_samples(basis, probes)
        if name == "zero_flow":
            uf = np.zeros_like(probes)
        elif name == "slow_affine":
            uf = (
                np.array(cfg["flow_velocity_m_s"])
                + (probes - center) @ np.array(cfg["flow_gradient_s_inv"]).T
            )
        else:
            fq = flow.sample(probes)
            assert fq["valid"].all()
            uf = fq["velocity_m_s"]
        phase = -np.angle(P[2])
        amp = np.full(4, cap)
        old = (
            copy.deepcopy(old5)
            if name == "Y_CFD"
            else dict(amplitude_Pa=[151.340514] * 4, phase_rad=[0.0] * 4, enabled=[True] * 4)
        )
        old.pop("mode", None)
        old.update(time_s=0.0, flow_enabled=cfg["flow_enabled"], secondary_enabled=True)
        events = dict(
            off=dict(old, enabled=[False] * 4, secondary_enabled=False),
            unchanged=old,
            target_phase=dict(old, amplitude_Pa=amp.tolist(), phase_rad=phase.tolist(), enabled=[True] * 4),
        )
        # One 7-variable bounded static fit. Desired local drift contracts toward
        # the fixed target; controls are not fitted to integrated positions.
        desired = (target - probes) * 2.0

        def unpack(z):
            return z[:4] * cap, np.r_[0.0, z[4:]]

        def objective(z):
            aa, pp = unpack(z)
            try:
                ua, pressure = primary_from_control(response, p, P, G, aa, pp)
            except ValueError:
                return 1e8
            weights = np.array([2.0, 1.0, 3.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
            err = np.sum((ua + uf - desired) ** 2, axis=1) / (10e-6) ** 2
            # A lower 100 Pa local pressure preference avoids solving drift by
            # merely switching OFF; this is an objective, not fabricated attraction.
            cohesion = np.maximum(100 - abs(pressure), 0) ** 2 / 100**2
            return float(np.mean(weights * (err + 0.1 * cohesion)))

        z0 = np.r_[np.ones(4), np.angle(np.exp(1j * (phase[1:] - phase[0])))]
        fit = minimize(
            objective,
            z0,
            method="L-BFGS-B",
            bounds=[(0.0, 1.0)] * 4 + [(-np.pi, np.pi)] * 3,
            options={"maxiter": 80, "maxfun": 800, "ftol": 1e-10, "gtol": 1e-8},
        )
        aa, pp = unpack(fit.x)
        rows = {}
        candidates = dict(
            events, optimized=dict(old, amplitude_Pa=aa.tolist(), phase_rad=pp.tolist(), enabled=[True] * 4)
        )
        for mode, event in candidates.items():
            amplitudes = np.array(event["amplitude_Pa"]) * np.array(event["enabled"])
            ua, pressure = primary_from_control(response, p, P, G, amplitudes, event["phase_rad"])
            u = ua + uf
            direction = delta / np.linalg.norm(delta)
            restoring = np.einsum("ij,ij->i", u[3:] - u[2], probes[3:] - target)
            rows[mode] = dict(
                probe_velocity_m_s=u,
                acoustic_velocity_m_s=ua,
                pressure_Pa=abs(pressure),
                initial_target_progress_m_s=u[0] @ direction,
                target_speed_m_s=np.linalg.norm(u[2]),
                all_six_restoring=bool((restoring < 0).all()),
                restoring_dot_m2_s=restoring,
            )
        eligible = (
            objective(fit.x) < objective(z0) - 1e-6 and rows["optimized"]["initial_target_progress_m_s"] > 0
        )
        if eligible:
            events["optimized"] = candidates["optimized"]
        cfg["events"] = events
        diagnostics.append(
            dict(
                environment=name,
                controls=rows,
                optimizer=dict(
                    success=fit.success,
                    message=str(fit.message),
                    iterations=fit.nit,
                    evaluations=fit.nfev,
                    initial_objective=objective(z0),
                    final_objective=objective(fit.x),
                    trajectory_eligible=eligible,
                ),
                probes_m=probes,
                flow_m_s=uf,
            )
        )
    write(out / "STATIC_PREFLIGHT.json", diagnostics)
    write(out / "PROTOCOL.json", proto)
    write(out / "IDENTITY.json", identity())
    print(
        json.dumps(
            clean(
                [
                    dict(
                        environment=d["environment"],
                        controls={
                            k: dict(
                                progress_um_s=v["initial_target_progress_m_s"] * 1e6,
                                target_speed_um_s=v["target_speed_m_s"] * 1e6,
                                restoring=v["all_six_restoring"],
                            )
                            for k, v in d["controls"].items()
                        },
                        optimizer=d["optimizer"],
                    )
                    for d in diagnostics
                ]
            ),
            indent=2,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
