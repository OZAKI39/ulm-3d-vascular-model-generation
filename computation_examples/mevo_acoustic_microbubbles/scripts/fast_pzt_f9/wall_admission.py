"""Check true-wall instrumentation through onset, plus rejected-interval isolation."""

import argparse
import numpy as np
import h5py
from mevo_mb.fast_sim.f9_ledger import AuditSimulator
from mevo_mb.fast_sim.flow_gpu import AffineFlow
from .checks import box_geometry
from .common import ROOT, BASE, objects, identity, js, write
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    assert js(BASE / "checks_v1/CHECKS.json")["status"] == "PASS"
    write(out / "IDENTITY.json", identity())
    write(
        out / "REGISTERED_PROTOCOL.json",
        dict(
            duration_s=0.03,
            reason="F8 first Y1 wall corrections begin around13ms; prior10ms check only tested the precontact Y segment. Cover actual wall events before full audits.",
            thresholds=js(ROOT / "configs/fast_pzt/f9/protocol.json"),
        ),
    )
    obj = objects()
    cfg = js(ROOT / "results/fast_pzt/f8/y_formal_v1/Y1/CONFIG.json")
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    summaries = {}
    for audit in [False, True]:
        name = "Y_wall_on" if audit else "Y_wall_off"
        summaries[name] = run_case(
            out / name, cfg["environment"], cfg["events"], proto, obj, duration=0.03, audit=audit
        )
        assert summaries[name]["numerical_quality"] == "PASS"
    with (
        h5py.File(out / "Y_wall_on/trajectory.h5") as on,
        h5py.File(out / "Y_wall_off/trajectory.h5") as off,
        h5py.File(ROOT / "results/fast_pzt/f8/y_formal_v1/Y1/trajectory.h5") as old,
    ):
        x = on["positions_m"][:]
        assert np.array_equal(x, off["positions_m"][:])
        assert np.array_equal(x, old["positions_m"][: len(x)])
        ne = len(on["motion_ledger/wall_contact_events"])
        assert ne > 0
    p, b, r, g, flow = obj
    g = box_geometry(p.radius_m)
    x = np.array([[0.0, 0.0, p.radius_m + 0.5e-9]])
    sim = AuditSimulator(
        x,
        p,
        r,
        flow=AffineFlow([0.0, 0.0, -0.1], np.zeros((3, 3)), np.zeros(3)),
        acoustic_enabled=False,
        vessel=g,
        audit_enabled=True,
    )
    event = dict(
        time_s=0.0,
        amplitude_Pa=[0.0] * 4,
        phase_rad=[0.0] * 4,
        enabled=[False] * 4,
        flow_enabled=True,
        secondary_enabled=False,
    )
    result = sim.run(1e-5, 1e-5, [event], sample_interval_s=1e-5)
    assert result["summary"]["status"] == "FAILED"
    assert len(result["motion_ledger"]["displacement_m"]) == 0
    assert len(result["motion_ledger"]["rejected_intervals"]) == 1
    with h5py.File(out / "expected_rejection.h5", "w") as h:
        h["positions_m"] = result["positions_m"]
        h["time_s"] = result["time_s"]
        for i, rec in enumerate(result["motion_ledger"]["rejected_intervals"]):
            for k, v in rec.items():
                h.create_dataset(f"rejected/{i}/{k}", data=v)
    write(out / "EXPECTED_REJECTION.json", result["summary"])
    write(
        out / "CHECKS.json",
        dict(
            status="PASS",
            actual_wall_event_count=ne,
            audit_on_off_bitwise_equal=True,
            F8_first_61_frames_bitwise_equal=True,
            expected_wall_budget_failure_isolated=True,
            performance=summaries,
        ),
    )
    print("Actual-wall and rejected-interval admission PASS", flush=True)


if __name__ == "__main__":
    main()
