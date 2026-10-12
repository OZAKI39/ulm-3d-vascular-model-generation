"""Closed-data F7 acceptance and reproducible plot-data extraction; CPU only."""

import json
import h5py
import numpy as np
from .common import ROOT, BASE, js, write, sha


def load(path):
    s = js(path / "SUMMARY.json")
    cfg = js(path / "CONFIG.json")
    with h5py.File(path / "trajectory.h5") as h:
        data = dict(
            t=h["time_s"][:],
            x=h["positions_m"][:, 0],
            c=h["metrics/global_centroid_m"][:, 0],
            rg=h["metrics/global_rg_m"][:, 0],
            fraction=h["metrics/largest_size"][:, 0] / 32,
            error=h["trajectory/target_error_m"][:],
            quality=h["quality/raw_F1_diagonal"][:, 0],
            contacts=h["metrics/hard_contact_pair_count"][:, 0],
        )
        assert np.isfinite(data["x"]).all()
        assert json.loads(h["control"].attrs["events_json"]) == s["controls"]
        np.testing.assert_allclose(
            data["error"], np.linalg.norm(data["c"] - cfg["environment"]["target_m"], axis=1), rtol=1e-12
        )
        assert h["active"][:].all()
        assert h["trajectory/physical_dt_s"][:].sum() > 0
    return s, cfg, data


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    paths = {m: BASE / "replay_v1" / m for m in ["M0", "M1", "M2"]}
    paths.update(
        old=ROOT / "results/fast_pzt/f6/navigation_v1/zero_flow_feedback_RPY",
        failed=BASE / "controller_v1/run",
        new=BASE / "controller_v2/run",
        half=BASE / "controller_half_v1/run",
    )
    paths["slow_failed"] = BASE / "slow_v2/run"
    paths["slow"] = BASE / "slow_v3/run"
    arrays = {}
    summaries = {}
    provenance = {}
    loaded = {}
    for key, path in paths.items():
        s, cfg, d = load(path)
        loaded[key] = (s, cfg, d)
        summaries[key] = dict(
            status=s["status"],
            reason=s.get("reason"),
            completed_physical_s=s["completed_physical_s"],
            navigation=s["navigation"],
            quality=s["quality"],
            numerical_quality=s["numerical_quality"],
            wall_s=s["wall_s"],
            process_body_wall_s=s["process_body_wall_s"],
        )
        for k, v in d.items():
            arrays[key + "_" + k] = v
        arrays[key + "_target"] = np.array(cfg["environment"]["target_m"])
        events = s["controls"]
        arrays[key + "_control_t"] = np.array([e["time_s"] for e in events])
        arrays[key + "_amplitudes"] = np.array([e["amplitude_Pa"] for e in events])
        arrays[key + "_phases"] = np.array([e["phase_rad"] for e in events])
        for field in ["time_s", "amplitude_Pa", "phase_rad", "enabled", "secondary_enabled", "flow_enabled"]:
            if key in ["M0", "M1", "M2"]:
                assert [e[field] for e in events] == [e[field] for e in proto["frozen_F6_events"]]
        provenance[key] = dict(
            path=str(path.relative_to(ROOT)),
            files={
                f: sha(path / f) for f in ["trajectory.h5", "SUMMARY.json", "CONFIG.json", "IDENTITY.json"]
            },
            execution_commit=js(path / "IDENTITY.json")["execution_commit"],
        )
        if "feedback_optimization" in s:
            opts = s["feedback_optimization"]
            costs = np.array([r["wall_s"] for r in opts])
            summaries[key]["optimizer"] = dict(
                total_s=costs.sum(),
                median_s=np.median(costs),
                max_s=costs.max(),
                evaluations=sum(r["evaluations"] for r in opts),
                max_evaluations=max(r["evaluations"] for r in opts),
                updates=len(opts),
                cost_scope="CPU optimizer + cached basis preparation + GPU candidates and transfers; excludes full-grid production control combination",
            )
        if key in ["new", "half", "slow"]:
            pred = []
            actual = []
            for i, rec in enumerate(s["feedback_optimization"]):
                t = rec["time_s"]
                end = events[i + 1]["time_s"] if i + 1 < len(events) else d["t"][-1]
                i0 = int(np.argmin(abs(d["t"] - t)))
                i1 = int(np.argmin(abs(d["t"] - end)))
                if i1 <= i0:
                    continue
                pred.append(rec["predicted_mean_velocity_m_s"])
                actual.append((d["c"][i1] - d["c"][i0]) / (d["t"][i1] - d["t"][i0]))
            pred = np.array(pred)
            actual = np.array(actual)
            arrays[key + "_predicted_velocity"] = pred
            arrays[key + "_actual_segment_velocity"] = actual
            err = np.linalg.norm(pred - actual, axis=1)
            summaries[key]["prediction"] = dict(
                rms_error_m_s=np.sqrt(np.mean(err**2)),
                max_error_m_s=err.max(),
                mean_error_m_s=err.mean(),
                scope="Free instantaneous RPY mean at event vs actual next 25ms centroid displacement, includes changing morphology and PBD",
            )
    a = loaded["new"][2]
    b = loaded["half"][2]
    np.testing.assert_allclose(a["t"], b["t"], atol=1e-12)
    dt = dict(
        max_centroid_m=np.linalg.norm(a["c"] - b["c"], axis=1).max(),
        max_target_error_m=np.abs(a["error"] - b["error"]).max(),
        max_Rg_relative=np.max(np.abs(a["rg"] - b["rg"]) / b["rg"]),
        max_fraction_difference=np.abs(a["fraction"] - b["fraction"]).max(),
        max_contact_count_difference=np.abs(a["contacts"] - b["contacts"]).max(),
        same_navigation_gate=loaded["new"][0]["navigation"]["arrived_and_held"]
        == loaded["half"][0]["navigation"]["arrived_and_held"],
        arrival_dwell={k: loaded[k][0]["navigation"] for k in ["new", "half"]},
    )
    limits = proto["f7_dt_comparison"]
    dt["passed"] = bool(
        dt["max_centroid_m"] <= limits["centroid_m"]
        and dt["max_target_error_m"] <= limits["target_error_m"]
        and dt["max_Rg_relative"] <= limits["Rg_relative"]
        and dt["max_fraction_difference"] <= limits["largest_fraction"]
        and dt["same_navigation_gate"]
        and all(loaded[k][0]["numerical_quality"] == "PASS" for k in ["new", "half"])
    )
    y = js(BASE / "reachability_v1/REACHABILITY.json")
    # The repair affects only the previously failing700ms control update.
    prefix = loaded["slow_failed"][2]["x"]
    np.testing.assert_array_equal(loaded["slow"][2]["x"][: len(prefix)], prefix)
    repair_events = [r for r in loaded["slow"][0]["feedback_optimization"] if r.get("feasibility_repair")]
    summaries["slow"]["feasibility_repairs"] = repair_events
    summaries["slow"]["earlier_failed_prefix_bitwise_equal"] = True
    new = loaded["new"][0]
    old = loaded["old"][0]
    gate = dict(
        software_checks=js(BASE / "fixed_v1/CHECKS.json")["status"],
        cached_predictor=js(BASE / "predictor_v1/CHECKS.json")["status"],
        frozen_replay=js(BASE / "replay_v1/REPLAY.json")["status"],
        all_full_numerical_pass=all(s["numerical_quality"] == "PASS" for s in summaries.values()),
        all_completed_trajectories_numerical_pass=all(
            s["numerical_quality"] == "PASS" for s in summaries.values() if s["status"] == "COMPLETE"
        ),
        incomplete_cases=[k for k, s in summaries.items() if s["status"] != "COMPLETE"],
        failure_retention_check=js(BASE / "failure_record_v1/CHECKS.json")["status"],
        timestep="PASS" if dt["passed"] else "FAIL",
        new_controller_navigation_pass=bool(
            new["navigation"]["arrived_and_held"] and new["numerical_quality"] == "PASS"
        ),
        old_controller_navigation_pass=bool(
            old["navigation"]["arrived_and_held"] and old["numerical_quality"] == "PASS"
        ),
        new_final_error_lower_than_old=new["navigation"]["final_error_m"]
        < old["navigation"]["final_error_m"],
        Y_full_admitted=y["full_Y_admitted"],
        slow_validation=summaries["slow"]["navigation"]["classification"]
        if summaries["slow"]["status"] == "COMPLETE"
        else "INCOMPLETE_BOUNDED_CONTROL_SEARCH_FAILED",
        evidence_level="EXPLORATORY",
        translation_calibration="UNCALIBRATED",
        physical_validation=False,
        figures="PENDING",
        archive="PENDING",
    )
    gate["software_numerical_acceptance"] = (
        "PASS"
        if all(
            gate[k] == "PASS" for k in ["software_checks", "cached_predictor", "frozen_replay", "timestep"]
        )
        and gate["all_completed_trajectories_numerical_pass"]
        and gate["failure_retention_check"] == "PASS"
        else "PARTIAL"
    )
    gate["acceptance_scope"] = (
        "Completed trajectories and explicit safe failure retention; incomplete slow-flow case is not a successful or complete trajectory. No physical validation."
    )
    np.savez_compressed(out / "plot_data.npz", **arrays)
    write(
        out / "RESULTS.json",
        dict(
            cases=summaries,
            timestep=dt,
            provenance=provenance,
            Y=y,
            gradient=js(BASE / "fixed_v1/CHECKS.json")["states"],
        ),
    )
    write(out / "FINAL_GATE.json", gate)
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
