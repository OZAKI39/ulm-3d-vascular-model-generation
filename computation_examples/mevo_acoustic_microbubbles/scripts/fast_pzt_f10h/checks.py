"""GPU adapter gates: existing static references, controls, graph and fail-stop."""

import argparse
import numpy as np
from dataclasses import asdict
from mevo_mb.fast_sim.f10h_lubrication import LubricationDiagnostic, MODES
from mevo_mb.fast_sim.f10f_smallsteps import SmallStepsBatch
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.f10g_contact import resistance, jacobi, balance, pbd
from scripts.fast_pzt_f9.common import objects
from .common import BASE, ROOT, CASES, state, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    obj = objects()
    p, b, r, g, f = obj
    sel = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")
    rows = []
    for case, chosen in zip(CASES, [sel["branch"], sel["preparation"]]):
        _, t, _ = state(case)
        x = np.array(t["initial_positions_m"])
        cfg = js(ROOT / chosen["relative_path"] / "CONFIG.json")
        ev = [e for e in cfg["events"] if e["time_s"] <= t["start_s"]][-1] | {"time_s": 0}
        for mode in ["M0_LUB2_PBD12_DT10", "M0_LUB4_PBD12_DT10"]:
            s = LubricationDiagnostic(x, p, r, b, vessel=g, flow=f, numerical_mode=mode)
            s.configure_schedule([ev], 1e-5)
            s.control(ev)
            try:
                s.advance(1e-5, 1, graphs=False)
            except RuntimeError:
                pass
            assert len(s.trace["lubrication"]) == MODES[mode][1]
            np.testing.assert_allclose(
                s.trace["free_velocity_m_s"], t["free_velocity_m_s"], rtol=0, atol=1e-12
            )
            op = resistance(x, p.radius_m)
            free = np.array(t["free_velocity_m_s"])
            u = jacobi(op, free, MODES[mode][1])
            ub = balance(op, free, u)
            ref, rounds = pbd(x + 1e-5 * ub, p.radius_m)
            err = float(np.abs(s.final_proposal.numpy() - ref).max())
            assert err < 1e-13
            if "LUB2" in mode:
                assert s.first_rejection and s.accepted_steps == 0
                np.testing.assert_array_equal(s.x.numpy(), x)
                try:
                    s.advance(1e-5, 1, graphs=False)
                except RuntimeError:
                    pass
                else:
                    raise AssertionError("rejected state advanced")
            else:
                assert s.accepted_steps == 1 and not s.first_rejection
            rows.append(
                dict(
                    case=case,
                    mode=mode,
                    cpu_position_error_m=err,
                    accepted_steps=s.accepted_steps,
                    first_rejection=s.first_rejection,
                    quality=s.quality(),
                )
            )
    # Same original start, a genuine control change at physical 20us; graph/observer parity.
    cfg = js(ROOT / sel["branch"]["relative_path"] / "CONFIG.json")
    x = np.array(cfg["initial_positions_m"])
    e = {k: v for k, v in cfg["events"][0].items() if k != "mode"}
    schedule = [e | {"time_s": 0}, e | {"time_s": 2e-5, "secondary_enabled": False}]
    for mode in MODES:
        dt = MODES[mode][0]
        states = []
        for diagnostic in [True, False]:
            s = LubricationDiagnostic(
                x, p, r, b, vessel=g, flow=f, numerical_mode=mode, diagnostic=diagnostic
            )
            s.configure_schedule(schedule, dt)
            s.control(schedule[0])
            s.advance(dt, round(4e-5 / dt), graphs=not diagnostic)
            assert len(s.executed_controls) == 2 and s.executed_controls[1]["time_s"] == 2e-5
            states.append(s.x.numpy().copy())
        np.testing.assert_array_equal(*states)
    # Original mode remains byte-identical to inherited F10F B0 on real initial state.
    s = SmallStepsBatch(x, p, r, b, vessel=g, flow=f)
    s.configure_schedule([e], 1e-5)
    s.control(e)
    s.advance(1e-5, 2, graphs=False)
    q = LubricationDiagnostic(x, p, r, b, vessel=g, flow=f)
    q.configure_schedule([e], 1e-5)
    q.control(e)
    q.advance(1e-5, 2, graphs=False)
    np.testing.assert_array_equal(s.x.numpy(), q.x.numpy())
    assert asdict(p) == asdict(F2Physics())
    try:
        LubricationDiagnostic(np.stack([x.reshape(-1, 3)] * 2), p, r, b, vessel=g, flow=f)
    except ValueError as exc:
        assert "batch not admitted" in str(exc)
    else:
        raise AssertionError("batch silently enabled")
    write(
        out / "RESULTS.json",
        dict(
            status="PASS",
            static=rows,
            controls_and_graph_parity=True,
            original_two_rounds_unchanged=True,
            batch_explicitly_rejected=True,
        ),
    )
    print("F10H ENGINEERING PASS", flush=True)


if __name__ == "__main__":
    main()
