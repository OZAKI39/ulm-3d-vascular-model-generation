"""Geometry CPU gates -> original full-environment static gate -> one bounded case."""

import argparse
import numpy as np
from mevo_mb.fast_sim.f10h_diagnostic import GeometricContactDiagnostic
from mevo_mb.fast_sim.f10g_contact import resistance
from mevo_mb.fast_sim.f10h_contact import geometric_contact
from scripts.fast_pzt_f9.common import objects
from .common import ROOT, BASE, CASES, state, js, write, identity
from .trials import record


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["short", "dt5"], default="short")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    source = identity()
    write(out / "IDENTITY.json", source)
    small = js(BASE / "geometry_v1/SMALL_SYSTEMS.json")
    assert small["two_pass"] and small["four_pass"]
    static = js(BASE / "geometry_v1/DENSE_STATIC.json")
    assert static["both_pass"]
    kc = js(ROOT / "configs/fast_pzt/f10h/protocol.json")["B2_kc_N_m"]
    obj = objects()
    p, b, r, g, f = obj
    sel = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")
    if a.phase == "short":
        gates = []
        for case, chosen in zip(CASES, [sel["branch"], sel["preparation"]]):
            _, t, _ = state(case)
            x = np.array(t["initial_positions_m"])
            cfg = js(ROOT / chosen["relative_path"] / "CONFIG.json")
            ev = [e for e in cfg["events"] if e["time_s"] <= t["start_s"]][-1] | {"time_s": 0}
            s = GeometricContactDiagnostic(x, p, r, b, vessel=g, flow=f, kc_N_m=kc)
            s.configure_schedule([ev], 1e-5)
            s.control(ev)
            s.advance(1e-5, 1, graphs=False)
            ref = next(v for v in static["rows"] if v["name"] == case)["result"]
            np.testing.assert_allclose(s.trace["total_force_N"], t["total_force_N"], atol=1e-20, rtol=0)
            np.testing.assert_allclose(s.x.numpy(), ref["positions_m"], atol=1e-13, rtol=0)
            assert s.quality()["error_bits"] == 0
            gates.append(dict(case=case, quality=s.quality(), step=s.soft_steps[-1]))
        # Full environment rejection keeps coordinates, active mask and clock unchanged.
        center = np.array(state(CASES[0])[1]["initial_positions_m"]).mean(0)
        x = np.array([center + [-p.radius_m + 200e-9, 0, 0], center + [p.radius_m - 200e-9, 0, 0]])
        s = GeometricContactDiagnostic(x, p, r, b, vessel=g, flow=f, kc_N_m=kc)
        s.configure_schedule([ev], 1e-5)
        s.control(ev)
        try:
            s.advance(1e-5, 1, graphs=False)
        except RuntimeError:
            pass
        else:
            raise AssertionError("invalid geometry accepted")
        assert s.accepted_steps == 0
        np.testing.assert_array_equal(s.x.numpy(), x)
        write(
            out / "FAIL_STOP_CHECK.json",
            dict(first_rejection=s.first_rejection, step=s.soft_steps[-1], last_legal=x),
        )
        # Two dt values, force-changing controls at 20/40us and independent CPU equations.
        checks = []
        for dt, mode in [(1e-5, "B2_GEOM_K2_DT10"), (5e-6, "B2_GEOM_K2_DT5")]:
            x = np.array([center + [-10e-6, 0, 0], center + [10e-6, 0, 0]])
            e = {key: value for key, value in ev.items() if key != "mode"}
            schedule = [
                e | {"time_s": 0},
                e | {"time_s": 2e-5, "secondary_enabled": False},
                e | {"time_s": 4e-5, "secondary_enabled": True},
            ]
            s = GeometricContactDiagnostic(x, p, r, b, vessel=g, flow=f, kc_N_m=kc, numerical_mode=mode)
            s.configure_schedule(schedule, dt)
            s.control(schedule[0])
            for _ in range(round(6e-5 / dt)):
                s.advance(dt, 1, graphs=False)
                row = s.soft_steps[-1]
                ref = geometric_contact(
                    row["initial_positions_m"],
                    row["result"]["free_m_s"],
                    p.radius_m,
                    6 * np.pi * p.viscosity_Pa_s * p.radius_m,
                    dt,
                    kc,
                    resistance(row["initial_positions_m"], p.radius_m)["A"],
                    atol_m=1e-15,
                )
                np.testing.assert_allclose(s.x.numpy(), ref["positions_m"], atol=1e-13, rtol=0)
            assert len(s.executed_controls) == 3
            checks.append(dict(mode=mode, accepted_steps=s.accepted_steps, controls=s.executed_controls))
        write(out / "STATIC_FULL_GATE.json", dict(status="PASS", rows=gates, control_checks=checks))
    else:
        previous = js(BASE / "b2_short_v1/RESULTS.json")["rows"][0]
        assert previous["status"] == "COMPLETE"
        cost = previous["timing"]["integration_quality_wall_s"] / previous["accepted_steps"]
        old = js(ROOT / "results/fast_pzt/f10g/integration_v1/RESULTS.json")["rows"][0]
        assert cost <= 10 * old["timing"]["integration_quality_wall_s"] / old["accepted_steps"]
    mode = "B2_GEOM_K2_DT10" if a.phase == "short" else "B2_GEOM_K2_DT5"
    chosen = sel["branch"]
    row = record(
        out / (chosen["case_id"] + "_" + mode),
        chosen,
        obj,
        mode,
        0.005,
        source,
        simulator=GeometricContactDiagnostic,
        kc_N_m=kc,
    )
    write(out / "RESULTS.json", dict(rows=[row], training_ready=False, full_10ms_admitted=False))
    print(row["case"], mode, row["status"], row["accepted_steps"], flush=True)


if __name__ == "__main__":
    main()
