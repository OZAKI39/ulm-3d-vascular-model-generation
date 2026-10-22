"""Necessary F8 audit reproductions and the missing Y0 half-step control."""

import argparse
import h5py
import numpy as np
from .common import ROOT, BASE, js, write, sha, identity, objects
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--case", choices=["Y0_5", "Y1_10", "Y1_5", "Y2_10", "Y2_5"], required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    assert js(BASE / "wall_admission_v1/CHECKS.json")["status"] == "PASS"
    condition, dtname = a.case.split("_")
    dt = 1e-5 if dtname == "10" else 5e-6
    historic = (
        ROOT
        / "results/fast_pzt/f8"
        / ("y_half_v1" if dtname == "5" and condition != "Y0" else "y_formal_v1")
        / condition
    )
    cfg = js(historic / "CONFIG.json")
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    reg = ROOT / "results/fast_pzt/f8/y_static_v1/TARGETS_BEFORE_TRAJECTORIES.json"
    proto["feedback_policy"] = js(reg)["policy"]
    if condition == "Y2":
        admission = js(BASE / "analysis/ACTIVE_AUDIT_ADMISSION.json")
        assert a.case in admission["admitted_cases"]
    write(out / "IDENTITY.json", identity())
    write(
        out / "REGISTERED_PROTOCOL.json",
        dict(
            case=a.case,
            dt_s=dt,
            duration_s=0.3,
            target_registration_sha256=sha(reg),
            f8_source=str(historic.relative_to(ROOT)),
            source_hashes={
                p.name: sha(p)
                for p in historic.iterdir()
                if p.name in ["CONFIG.json", "SUMMARY.json", "trajectory.h5", "GATE_CROSSING.json"]
            },
            new_independent_sample=False,
            role="new_missing_half_step_control" if a.case == "Y0_5" else "instrumented_audit_reproduction",
            protocol=js(ROOT / "configs/fast_pzt/f9/protocol.json"),
        ),
    )
    s = run_case(
        out / "run",
        cfg["environment"],
        cfg["events"],
        proto,
        objects(),
        dt=dt,
        duration=0.3,
        feedback="stokes" if condition == "Y2" else False,
        audit=True,
    )
    assert (
        s["status"] == "COMPLETE"
        and s["numerical_quality"] == "PASS"
        and s["ledger_validation"]["status"] == "PASS"
    )
    if a.case != "Y0_5":
        with h5py.File(out / "run/trajectory.h5") as new, h5py.File(historic / "trajectory.h5") as old:
            bitwise = np.array_equal(new["positions_m"][:], old["positions_m"][:])
            assert bitwise
            np.testing.assert_array_equal(new["time_s"][:], old["time_s"][:])
        oldgate = js(historic / "GATE_CROSSING.json")
        newgate = js(out / "run/GATE_CROSSING.json")
        assert oldgate == newgate, "Gate observer must remain identical"
        report = dict(
            status="PASS",
            positions_bitwise_equal=True,
            gate_exact_equal=True,
            actual_controls_equal=js(historic / "SUMMARY.json")["controls"] == s["controls"],
            scope="Audit reproduction, not an independent sample",
        )
        assert report["actual_controls_equal"]
    else:
        report = dict(
            status="PASS",
            scope="New Y0 5us numerical control; differs from F8 Y0 only in dt plus passive audit",
            positions_bitwise_equal=None,
        )
    write(out / "REPRODUCTION.json", report)
    write(out / "COMPLETE.json", dict(status="COMPLETE", case=a.case, summary=s))
    print(report, flush=True)


if __name__ == "__main__":
    main()
