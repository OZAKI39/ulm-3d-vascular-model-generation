"""Preregister an original-Y downstream gate, then run admitted matched controls."""

import argparse
from datetime import datetime, timezone
import time
import numpy as np
from mevo_mb.fast_sim.f8_navigation import MatchedFeedback
from mevo_mb.fast_sim.f8_gate import GateSimulator
from mevo_mb.fast_sim.f6_engine import F6Simulator
from scripts.fast_pzt_f7.fixed import probe
from .common import ROOT, BASE, objects, write, js, identity, sha
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["static", "pilot", "formal", "half", "rpy"], required=True)
    ap.add_argument("--dt", type=float, default=1e-5)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    obj = objects()
    p, b, r, g, flow = obj
    start = time.perf_counter()
    ident = identity()
    write(out / "IDENTITY.json", ident)
    env = dict(proto["environment_Y"])
    design = proto["Y_design"]
    x = np.asarray(env["positions_m"])
    center = x.mean(0)
    target = center + np.r_[design["gate_offset_x_m"], design["target_offset_yz_m"]]
    env.update(
        target_m=target,
        x_gate_m=target[0],
        duration_s=design["duration_s"],
        task_type="flow_guided_transverse",
        target_tolerance_m=design["target_tolerance_m"],
    )
    # Dedicated transverse task normalization, identical for every Y2 evaluation.
    proto["feedback_policy"] = dict(proto["feedback_policy"], velocity_scale_m_s=1e-6)
    if a.phase == "static":
        assert js(BASE / "a_static_v1/CHECKS.json")["status"] == "PASS"
        initial_flow = flow.sample(x)
        assert initial_flow["valid"].all()
        flight = design["gate_offset_x_m"] / initial_flow["velocity_m_s"].mean(0)[0]
        # These are geometry probes, not actual or plotted integrated bubble positions.
        corridor = np.concatenate([x + u * (target - center) for u in np.linspace(0, 1, 9)])
        offsets = (
            np.array([[0, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]) * design["target_tolerance_m"]
        )
        test_points = np.vstack([corridor, target + offsets])
        geom = g.query(test_points)
        field_upper = b.origin + b.spacing * (np.array(b.shape) - 1)
        margin = p.radius_m + design["safe_radius_margin_m"]
        field_valid = np.all(
            (test_points >= b.origin + margin) & (test_points <= field_upper - margin), axis=1
        )
        cf = flow.sample(test_points)
        geometry_valid = np.all(geom["inside"]) and np.min(geom["distance_m"] - p.radius_m) >= 0
        support_valid = bool(geometry_valid and field_valid.all() and cf["valid"].all())
        registration = dict(
            registered_UTC=datetime.now(timezone.utc).isoformat(),
            identity=ident,
            task_type="flow_guided_transverse",
            environment=env,
            policy=proto["feedback_policy"],
            design=design,
            particle_ids=list(range(32)),
            initial_positions_m=x,
            target_m=target,
            estimated_CFD_flight_s=flight,
            initial_wall_gap_m=g.query(x)["distance_m"] - p.radius_m,
            target_wall_gap_m=g.query(target[None])["distance_m"][0] - p.radius_m,
            supported_geometry_probes=support_valid,
            min_corridor_gap_m=np.min(geom["distance_m"] - p.radius_m),
            source="Exact F5 frame600 at300ms copied into F6/F7 registered environment; original CFD and open Y lumen",
            rationale="20um downstream uses ~0.248s at F7 measured80.79um/s, within0.3s. +Z0.25um target and0.15um tolerance fixed from F7 directional probes before any F8 trajectory. Minimum paired effect50nm; additionally >3x new paired discretization bound. This is model-discretization resolution, not physical calibration or acoustic mesh accuracy.",
            plans=[
                "Y0 CFD sound OFF",
                "Y1 CFD original fixed PZT",
                "Y2 CFD target feedback; production M0, diagonal predictor with same force sampling",
            ],
            crossing="Fixed32-ID centroid crossing, linear intersection of consecutive accepted physical steps; no50ms dwell criterion",
        )
        write(out / "TARGETS_BEFORE_TRAJECTORIES.json", registration)
        rows = []
        for mode in ["M0", "M2"]:
            controls = {"Y0": env["events"]["off"], "Y1": env["events"]["unchanged"]}
            policy = dict(
                proto["feedback_policy"],
                flow_enabled=True,
                task_type=env["task_type"],
                prediction_mobility="diagonal_stokes" if mode == "M0" else "rpy_free_space",
            )
            controller = MatchedFeedback(b, r, p, target, lambda q: flow.sample(q)["velocity_m_s"], policy)
            controls["Y2"] = controller.update(x, 0.0)
            means = {}
            for label, event in controls.items():
                actual = probe(x, event, obj, mode, True, True)
                means[label] = actual["velocity_m_s"].mean(0)
                rows.append(
                    dict(
                        production_mode=mode,
                        condition=label,
                        event=event,
                        actual_mean_m_s=means[label],
                        acoustic_mean_m_s=means[label] - initial_flow["velocity_m_s"].mean(0),
                        per_bubble_velocity_m_s=actual["velocity_m_s"],
                        optimization=controller.records[-1] if label == "Y2" else None,
                    )
                )
            for row in rows[-3:]:
                row["estimated_Z_gain_vs_Y0_m"] = flight * (row["actual_mean_m_s"][2] - means["Y0"][2])
                row["estimated_Z_gain_vs_Y1_m"] = flight * (row["actual_mean_m_s"][2] - means["Y1"][2])
        admission = {}
        for mode in ["M0", "M2"]:
            row = next(q for q in rows if q["production_mode"] == mode and q["condition"] == "Y2")
            admission[mode] = bool(
                support_valid
                and flight < design["duration_s"]
                and min(row["estimated_Z_gain_vs_Y0_m"], row["estimated_Z_gain_vs_Y1_m"])
                > design["static_admission_min_effect_m"]
            )
        # Actual tiny OFF integration verifies observer transparency and dt bracket.
        gate = center[0] + 4e-9
        opts = dict(
            flow=flow,
            flow_enabled=True,
            vessel=g,
            vessel_enabled=True,
            mobility_mode="f5_local",
            acoustic_enabled=False,
        )
        plain = F6Simulator(x, p, r, b, **opts)
        watched = GateSimulator(x, p, r, b, x_gate=gate, **opts)
        off = dict(env["events"]["off"], flow_enabled=True)
        a0 = plain.run(0.0001, 1e-5, [off], 1e-5)
        a1 = watched.run(0.0001, 1e-5, [off], 1e-5)
        np.testing.assert_array_equal(a0["positions_m"], a1["positions_m"])
        record = watched.gate_record()
        assert record["status"] == "CROSSED" and record["bracket_width_s"] <= 1.00001e-5
        idx = int(round(record["bracket_s"][0] / 1e-5))
        np.testing.assert_array_equal(record["accepted_before_m"], a1["positions_m"][idx, 0])
        np.testing.assert_array_equal(record["accepted_after_m"], a1["positions_m"][idx + 1, 0])
        write(
            out / "GATE_OBSERVER_CHECK.json",
            dict(
                status="PASS",
                gate=record,
                untouched_positions_bitwise_equal=True,
                actual_accepted_endpoint_indices=[idx, idx + 1],
            ),
        )
        np.savez_compressed(
            out / "static_data.npz",
            positions_m=x,
            flow_velocity_m_s=initial_flow["velocity_m_s"],
            geometry_probe_positions_m=test_points,
            geometry_probe_gap_m=geom["distance_m"] - p.radius_m,
        )
        write(
            out / "CHECKS.json",
            dict(
                status="PASS",
                admitted=admission,
                rows=rows,
                registration_sha256=sha(out / "TARGETS_BEFORE_TRAJECTORIES.json"),
                process_wall_s=time.perf_counter() - start,
                scope="Finite initial-state controls and geometry probes; not global reachability, wall-hydrodynamic validation or integrated future positions.",
            ),
        )
        print("Y static admission:", admission, flush=True)
        return
    regpath = BASE / "y_static_v1/TARGETS_BEFORE_TRAJECTORIES.json"
    registration = js(regpath)
    env = registration["environment"]
    proto["feedback_policy"] = registration["policy"]
    checks = js(BASE / "y_static_v1/CHECKS.json")
    required = "M2" if a.phase == "rpy" else "M0"
    assert checks["status"] == "PASS" and checks["admitted"][required]
    if a.phase in ["formal", "half"]:
        for label in ["Y0", "Y1", "Y2"]:
            assert js(BASE / f"y_pilot_v1/{label}/SUMMARY.json")["numerical_quality"] == "PASS"
    rows = {}
    labels = ["Y1", "Y2"] if a.phase == "half" else ["Y2"] if a.phase == "rpy" else ["Y0", "Y1", "Y2"]
    write(
        out / "REGISTERED_PROTOCOL.json",
        dict(
            proto=proto, target_registration_sha256=sha(regpath), phase=a.phase, dt_s=a.dt, conditions=labels
        ),
    )
    for label in labels:
        event = env["events"]["off" if label == "Y0" else "unchanged"]
        rows[label] = run_case(
            out / label,
            env,
            event,
            proto,
            obj,
            dt=a.dt,
            duration=0.03 if a.phase == "pilot" else env["duration_s"],
            feedback=("rpy" if a.phase == "rpy" else "stokes") if label == "Y2" else False,
            mobility="rpy_free_space" if a.phase == "rpy" else "f5_local",
        )
    write(
        out / "CAMPAIGN.json", dict(status="COMPLETE", runs=rows, process_wall_s=time.perf_counter() - start)
    )


if __name__ == "__main__":
    main()
