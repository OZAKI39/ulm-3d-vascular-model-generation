"""Bounded static original-Y authority and explicit virtual-neighborhood probes."""

import argparse
import time
import h5py
import numpy as np
from scipy.optimize import minimize
from mevo_mb.fast_sim.f7_predictor import CollectivePredictor
from mevo_mb.fast_sim.f7_navigation import MobilityFeedback
from .common import ROOT, BASE, write, js, identity, objects
from .fixed import local_numpy, independent_rpy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    assert js(BASE / "predictor_v1/CHECKS.json")["status"] == "PASS"
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    obj = objects()
    p, b, r, g, f = obj
    env = proto["environments"][2]
    x = np.array(env["positions_m"])
    target = np.array(env["target_m"])
    center = x.mean(axis=0)
    direction = target - center
    direction /= np.linalg.norm(direction)
    uf = f.sample(x)
    assert uf["valid"].all()
    pred = CollectivePredictor(b, r, p)
    pred.prepare(x, uf["velocity_m_s"])
    cap = proto["amplitude_cap_Pa"]
    rows = []
    candidates = []

    def save_candidate(name, amp, phase, enabled=None, detail=False):
        v = pred.evaluate(amp, phase, enabled, details=True)
        force = v["primary_N"] + v["secondary_N"]
        free = uf["velocity_m_s"] + force / (6 * np.pi * p.viscosity_Pa_s * p.radius_m)
        modes = {"M0": local_numpy(x, free, p), "M1": free, "M2": v["velocity_m_s"]}
        row = dict(
            name=name,
            amplitude_Pa=amp,
            phase_rad=phase,
            enabled=enabled if enabled is not None else [True] * 4,
            models={
                m: dict(
                    mean_velocity_m_s=z.mean(axis=0),
                    acoustic_mean_m_s=(z - uf["velocity_m_s"]).mean(axis=0),
                    target_progress_m_s=z.mean(axis=0) @ direction,
                )
                for m, z in modes.items()
            },
            min_pressure_Pa=np.linalg.norm(v["pressure_complex_Pa"], axis=1).min(),
            max_pressure_Pa=np.linalg.norm(v["pressure_complex_Pa"], axis=1).max(),
        )
        rows.append(row)
        if detail:
            return row, v
        return row

    for name in ["off", "unchanged", "target_phase"]:
        e = env["events"][name]
        save_candidate(name, e["amplitude_Pa"], e["phase_rad"], e["enabled"])
    start = time.perf_counter()
    # Six directional probes and one bounded 3-D hold fit, each a single warm start.
    for axis, sign in [(i, s) for i in range(3) for s in [-1, 1]] + [(3, 0)]:
        best = [np.inf, None]
        invalid = [0]

        def objective(z):
            try:
                v = pred.evaluate(z[:4] * cap, np.r_[0.0, z[4:]]).mean(axis=0)
                cost = float(-sign * v[axis] / 1e-4) if axis < 3 else float(np.sum((v / 1e-4) ** 2))
                candidates.append(
                    dict(
                        axis=axis,
                        sign=sign,
                        amplitude_Pa=z[:4] * cap,
                        phase_rad=np.r_[0.0, z[4:]],
                        mean_velocity_m_s=v,
                        cost=cost,
                    )
                )
                if cost < best[0]:
                    best[:] = [cost, z.copy()]
                return cost
            except ValueError:
                invalid[0] += 1
                return 1e12

        minimize(
            objective,
            np.r_[np.ones(4), np.zeros(3)],
            method="L-BFGS-B",
            bounds=[(0.0, 1.0)] * 4 + [(-np.pi, np.pi)] * 3,
            options=dict(maxiter=20, maxfun=160, ftol=1e-10, gtol=1e-8),
        )
        if best[1] is not None:
            z = best[1]
            row = save_candidate(
                ("XYZ"[axis] + ("+" if sign > 0 else "-")) if axis < 3 else "hold_fit",
                z[:4] * cap,
                np.r_[0.0, z[4:]],
            )
            row["invalid_candidates"] = invalid[0]
    positions = np.vstack([center, target, target + np.array(proto["stability_offsets_m"])])
    geom = g.query(positions)
    flow = f.sample(positions)
    lower = b.origin
    upper = lower + b.spacing * (np.array(b.shape) - 1)
    fieldvalid = np.all((positions >= lower) & (positions <= upper), axis=1)
    values = np.array([c["mean_velocity_m_s"] for c in candidates])
    tested = np.vstack([values, [v["models"]["M2"]["mean_velocity_m_s"] for v in rows]])
    admitted = bool(np.any((tested @ direction > 0) & (np.abs(tested[:, 0]) < 40e-6)))
    report = dict(
        identity=identity(),
        status="COMPLETE",
        positions_m=x,
        target_m=target,
        initial_centroid_m=center,
        target_direction=direction,
        flow_mean_m_s=uf["velocity_m_s"].mean(axis=0),
        rows=rows,
        candidates=candidates,
        tested_total_velocity_min_m_s=tested.min(axis=0),
        tested_total_velocity_max_m_s=tested.max(axis=0),
        neighborhood=dict(
            positions_m=positions,
            geometry=geom,
            flow=flow,
            field_valid=fieldvalid,
            field_lower_m=lower,
            field_upper_m=upper,
        ),
        full_Y_admitted=admitted,
        criterion=proto["Y_full_admission"],
        optimizer_wall_s=time.perf_counter() - start,
        scope="Finite candidate set at exact F5 300 ms state. Free-space RPY, no wall hydrodynamics. Directional extrema are observed candidates, not global bounds or proof of reachability.",
    )
    np.savez_compressed(
        out / "Y_state.npz",
        positions_m=x,
        flow_velocity_m_s=uf["velocity_m_s"],
        flow_gradient_s_inv=uf["gradient_s_inv"],
        wall_gap_m=g.query(x)["distance_m"] - p.radius_m,
    )
    write(out / "REACHABILITY.json", report)
    fixed = np.load(BASE / "fixed_v1/fixed_states.npz")
    decomposition = {}
    for name in ["mature", "dispersed", "near_wall"]:
        xx = fixed[name + "_positions_m"]
        mm = independent_rpy(xx, p.radius_m, p.viscosity_Pa_s)
        decomposition[name] = {
            key: (mm @ fixed[name + "_M2_" + key].ravel()).reshape(-1, 3).mean(axis=0)
            for key in ["primary_N", "secondary_N"]
        }
    write(
        out / "RPY_FORCE_CONTRIBUTIONS.json",
        dict(
            mean_velocity_m_s=decomposition,
            scope="Linear free RPY contributions at the fixed actual states. Secondary forces sum to zero, but configuration-dependent mobility row sums need not produce zero group velocity. PBD does not use the RPY constraint metric.",
        ),
    )
    # Virtual probes explicitly separate from actual integrated coordinates.
    formal = BASE / "controller_v2/run"
    if (formal / "SUMMARY.json").exists():
        summary = js(formal / "SUMMARY.json")
        with h5py.File(formal / "trajectory.h5") as h:
            final = h["positions_m"][-1, 0]
        target0 = np.array(proto["environments"][0]["target_m"])
        last = summary["controls"][-1]
        policy = dict(proto["f7_controller"], flow_enabled=False)
        stability = []
        for offset in proto["stability_offsets_m"]:
            virtual = final - final.mean(axis=0) + target0 + offset
            pred.prepare(virtual, np.zeros_like(virtual))
            frozen = pred.evaluate(last["amplitude_Pa"], last["phase_rad"]).mean(axis=0)
            controller = MobilityFeedback(b, r, p, target0, lambda q: np.zeros_like(q), policy)
            phase = np.array(last["phase_rad"])
            controller.previous = np.r_[
                np.array(last["amplitude_Pa"]) / cap, np.angle(np.exp(1j * (phase[1:] - phase[0])))
            ]
            event = controller.update(virtual, 0.0)
            updated = np.array(controller.records[-1]["predicted_mean_velocity_m_s"])
            stability.append(
                dict(
                    offset_m=offset,
                    frozen_velocity_m_s=frozen,
                    updated_velocity_m_s=updated,
                    frozen_restoring_dot_m2_s=-np.array(offset) @ frozen,
                    updated_restoring_dot_m2_s=-np.array(offset) @ updated,
                    event=event,
                    optimization=controller.records[-1],
                )
            )
        write(
            out / "STABILITY.json",
            dict(
                rows=stability,
                scope="Virtual translations of saved final relative morphology to fixed target +/-2 um; instantaneous free RPY drift, frozen and recomputed controls. No integrated recovery or Lyapunov proof.",
            ),
        )
    print("Y precheck full admission:", admitted, flush=True)


if __name__ == "__main__":
    main()
