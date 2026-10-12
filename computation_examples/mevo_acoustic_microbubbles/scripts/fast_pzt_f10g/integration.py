"""Static full-environment gates, then at most two original-state short replays."""

import argparse
from time import perf_counter
import numpy as np
from mevo_mb.fast_sim.f10g_diagnostic import SoftContactDiagnostic
from mevo_mb.fast_sim.f10d_metrics import process
from scripts.fast_pzt_f9.common import objects
from scripts.fast_pzt_f10c.runner import run
from .common import BASE, ROOT, CASES, state, js, write, identity, sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    begin = perf_counter()
    assert js(BASE / "prototype_v1/SMALL_SYSTEMS.json")["SOFT_CONTACT_FOUR_BUBBLE_PASS"]
    selected = [
        r
        for r in js(BASE / "prototype_v1/DENSE_STATIC.json")["rows"]
        if r["mode"] == "B2" and r["multiplier"] == 2
    ]
    assert len(selected) == 2 and all(r["contact_budget_pass"] for r in selected)
    source = identity()
    write(out / "IDENTITY.json", source)
    obj = objects()
    p, b, r, g, f = obj
    k = 2 * js(ROOT / "configs/fast_pzt/f10g/stiffness.json")["kc_N_m"]
    selection = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")
    static = []
    for case, chosen in zip(CASES, [selection["branch"], selection["preparation"]]):
        _, t, _ = state(case)
        x = np.array(t["initial_positions_m"])
        cfg = js(ROOT / chosen["relative_path"] / "CONFIG.json")
        ev = [e for e in cfg["events"] if e["time_s"] <= t["start_s"]][-1] | {"time_s": 0}
        s = SoftContactDiagnostic(x, p, r, b, vessel=g, flow=f, kc_N_m=k)
        s.configure_schedule([ev], 1e-5)
        s.control(ev)
        s.advance(1e-5, 1, graphs=False)
        result = s.soft_steps[-1]
        assert np.max(np.abs(s.trace["total_force_N"] - t["total_force_N"])) < 1e-20
        assert np.max(np.abs(s.trace["free_velocity_m_s"] - t["free_velocity_m_s"])) < 1e-12
        ref = [v for v in selected if v["case"] == case][0]["result"]
        delta = np.max(np.abs(s.x.numpy() - ref["positions_m"]))
        assert delta < 1e-13 and s.quality()["error_bits"] == 0
        assert all(np.max(np.abs(w["correction_m"])) == 0 for w in s.trace["wall_passes"])
        static.append(
            dict(
                case=case,
                full_environment_pass=True,
                quality=s.quality(),
                cpu_position_max_error_m=delta,
                step=result,
                controls=ev,
            )
        )
    write(
        out / "STATIC_FULL_GATE.json",
        dict(
            rows=static,
            SOFT_CONTACT_32_STATIC_PASS=True,
            scope="Independent saved-input reconstruction plus original GPU force/CFD/wall/domain acceptance; no new production model",
        ),
    )
    # Deliberately invalid two-bubble proposal verifies transactional fail-stop.
    center = np.array(state(CASES[0])[1]["initial_positions_m"]).mean(axis=0)
    x = np.array([center + [-p.radius_m + 200e-9, 0, 0], center + [p.radius_m - 200e-9, 0, 0]])
    s = SoftContactDiagnostic(x, p, r, b, vessel=g, flow=f, kc_N_m=k)
    ev = static[0]["controls"]
    s.configure_schedule([ev], 1e-5)
    s.control(ev)
    rejected = False
    try:
        s.advance(1e-5, 1, graphs=False)
    except RuntimeError:
        rejected = True
    assert rejected and s.accepted_steps == 0 and np.array_equal(s.x.numpy(), x)
    try:
        s.advance(1e-5, 1, graphs=False)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Rejected state advanced")
    write(
        out / "FAIL_STOP_CHECK.json",
        dict(pass_=True, first_rejection=s.first_rejection, last_legal_state_m=x),
    )
    # All preregistered necessary gates passed. CPU hybrid is bounded evidence,
    # intentionally not CUDA production acceleration or batch integration.
    rows = []
    for chosen in [selection["branch"], selection["preparation"]]:
        parent = ROOT / chosen["relative_path"]
        assert sha(parent / "trajectory.h5") == chosen["sha256"]
        cfg = js(parent / "CONFIG.json")
        summ = js(parent / "SUMMARY.json")
        path = out / (chosen["case_id"] + "_B2_K2")
        timing = run(
            path,
            np.array(cfg["initial_positions_m"]),
            cfg["events"],
            obj,
            duration=chosen["target_end_s"],
            graphs=False,
            simulator_class=SoftContactDiagnostic,
            substep_protocol=True,
            export_observer=True,
            kc_N_m=k,
            sample_ids=summ["sample_id"],
            postprocess=process,
            identity_provider=lambda: source,
            dataset_metadata=dict(
                initial_state_family=summ["initial_state_family"],
                split_groups=summ["split_groups"],
                split_group=summ["split_group"],
                labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED", "CPU_CONTACT_DIAGNOSTIC"],
                independent_scientific_initial_states=1,
                initial_state_sources=[dict(parent_sha256=chosen["sha256"], original_initial_state=True)],
                original_physics_mode="F10G_B2_K2_DT10",
                postprocessing_mode="F10D_optimized",
            ),
        )
        ob = js(path / "OBSERVER.json")
        summary = js(path / "SUMMARY.json")
        rows.append(
            dict(
                case=chosen["case_id"],
                status=summary["status"],
                attempted_steps=ob["attempted_steps"],
                accepted_steps=ob["accepted_steps"],
                first_rejection=ob["first_rejection_step"],
                timing=timing,
            )
        )
        write(
            out / "RESULTS.json",
            dict(
                rows=rows,
                process_body_s=perf_counter() - begin,
                full_10ms_admitted=False,
                training_ready=False,
            ),
        )
        print(rows[-1]["case"], summary["status"], ob["attempted_steps"], flush=True)
    write(
        out / "STOP_DECISION.json",
        dict(
            no_more_trajectories=True,
            reason="Two preregistered short windows exhausted; any failure forbids extension. CPU dense reference has no production speed or morphology admission.",
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
