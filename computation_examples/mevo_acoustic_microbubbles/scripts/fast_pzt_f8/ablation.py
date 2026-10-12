"""One-factor controller decisions and a bounded paired RPY trajectory study."""

import argparse
import time
import h5py
import numpy as np
from mevo_mb.fast_sim.f8_navigation import MatchedFeedback
from scripts.fast_pzt_f7.fixed import probe
from .common import ROOT, BASE, objects, write, js, identity, sha
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["static", "pilot", "formal"], required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    write(out / "IDENTITY.json", identity())
    write(out / "REGISTERED_PROTOCOL.json", proto)
    begin = time.perf_counter()
    obj = objects()
    p, b, r, _, _ = obj
    env = proto["environment_A"]
    spec = proto["A"]
    if a.phase == "static":
        rows = []
        arrays = {}
        source = ROOT / spec["static_source"]
        with h5py.File(source) as h:
            for frame in spec["static_frames"]:
                x = h["positions_m"][frame, 0]
                t = float(h["time_s"][frame])
                pair = {}
                arrays[f"x_{frame}"] = x
                for mode in ["A0", "A1"]:
                    policy = dict(
                        proto["feedback_policy"],
                        flow_enabled=False,
                        prediction_mobility="diagonal_stokes" if mode == "A0" else "rpy_free_space",
                    )
                    c = MatchedFeedback(b, r, p, env["target_m"], lambda q: np.zeros_like(q), policy)
                    event = c.update(x, t)
                    predicted = c.predictor.evaluate(event["amplitude_Pa"], event["phase_rad"], details=True)
                    actual = probe(x, event, obj, "M2")
                    diagonal = (actual["primary_N"] + actual["secondary_N"]) / (
                        6 * np.pi * p.viscosity_Pa_s * p.radius_m
                    )
                    reference = diagonal if mode == "A0" else actual["velocity_m_s"]
                    np.testing.assert_allclose(predicted["velocity_m_s"], reference, rtol=1e-9, atol=1e-16)
                    for key in ["primary_N", "secondary_N"]:
                        np.testing.assert_allclose(predicted[key], actual[key], rtol=1e-9, atol=1e-23)
                    arrays[f"{mode}_{frame}_predicted_m_s"] = predicted["velocity_m_s"]
                    arrays[f"{mode}_{frame}_actual_m_s"] = actual["velocity_m_s"]
                    desired = policy["gain_s_inv"] * (np.array(env["target_m"]) - x.mean(0))
                    direction = np.array(env["target_m"]) - x.mean(0)
                    direction /= np.linalg.norm(direction)
                    real = actual["velocity_m_s"].mean(0)
                    predicted_mean = predicted["velocity_m_s"].mean(0)
                    pair[mode] = dict(
                        policy=policy,
                        event=event,
                        optimization=c.records[-1],
                        predicted_mean_m_s=predicted_mean,
                        actual_RPY_mean_m_s=real,
                        mean_prediction_residual_m_s=np.linalg.norm(predicted_mean - real),
                        actual_desired_residual_m_s=np.linalg.norm(real - desired),
                        actual_target_progress_m_s=real @ direction,
                        actual_transverse_m_s=np.linalg.norm(real - (real @ direction) * direction),
                    )
                delta = np.linalg.norm(
                    np.array(pair["A0"]["actual_RPY_mean_m_s"]) - pair["A1"]["actual_RPY_mean_m_s"]
                )
                rows.append(
                    dict(frame=frame, time_s=t, modes=pair, selected_actual_mean_difference_m_s=delta)
                )
        admitted = (
            max(row["selected_actual_mean_difference_m_s"] for row in rows)
            > spec["full_admission_mean_velocity_difference_m_s"]
        )
        np.savez_compressed(out / "decisions.npz", **arrays)
        write(
            out / "CHECKS.json",
            dict(
                status="PASS",
                rows=rows,
                full_pair_admitted=admitted,
                source=str(source.relative_to(ROOT)),
                source_sha256=sha(source),
                initialization="Identical all-one normalized amplitudes and zero relative phases for each independent state decision; no trajectory warm start.",
                scope="Instantaneous unconstrained velocities; identical all-bubble forces, optimizer, objective, bounds and repairs. Only candidate mobility differs.",
                process_wall_s=time.perf_counter() - begin,
            ),
        )
        print("A static PASS; pair admitted:", admitted, flush=True)
        return
    checks = js(BASE / "a_static_v1/CHECKS.json")
    assert checks["status"] == "PASS" and checks["full_pair_admitted"]
    if a.phase == "formal":
        for mode in ["A0", "A1"]:
            assert js(BASE / f"a_pilot_v1/{mode}/SUMMARY.json")["numerical_quality"] == "PASS"
    rows = {}
    event = js(ROOT / "configs/fast_pzt/f7/protocol.json")["frozen_F6_events"][0]
    for mode in ["A0", "A1"]:
        rows[mode] = run_case(
            out / mode,
            env,
            event,
            proto,
            obj,
            feedback="stokes" if mode == "A0" else "rpy",
            mobility="rpy_free_space",
            duration=spec["pilot_duration_s"] if a.phase == "pilot" else spec["duration_s"],
            dt=spec["dt_s"],
        )
    write(
        out / "CAMPAIGN.json",
        dict(status="COMPLETE", phase=a.phase, runs=rows, process_wall_s=time.perf_counter() - begin),
    )


if __name__ == "__main__":
    main()
