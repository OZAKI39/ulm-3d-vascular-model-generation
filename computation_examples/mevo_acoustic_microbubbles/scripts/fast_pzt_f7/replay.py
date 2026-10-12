"""Frozen F6 control replay; pilots precede three formal trajectories."""

import argparse
import time
import h5py
import numpy as np
from .common import ROOT, BASE, MODES, write, js, identity, objects
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--full", action="store_true")
    a = ap.parse_args()
    assert js(BASE / "fixed_v1/CHECKS.json")["status"] == "PASS"
    if a.full:
        assert js(BASE / "pilot_v1/REPLAY.json")["status"] == "PASS"
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f7/protocol.json")
    env = proto["environments"][0]
    duration = 1.0 if a.full else proto["pilot_duration_s"]
    events = [e for e in proto["frozen_F6_events"] if e["time_s"] < duration]
    write(out / "REGISTERED_PROTOCOL.json", proto)
    write(out / "IDENTITY.json", identity())
    start = time.perf_counter()
    obj = objects()
    rows = {}
    for mode, mobility in MODES.items():
        rows[mode] = run_case(out / mode, env, events, proto, obj, duration=duration, mobility=mobility)
    with (
        h5py.File(out / "M0/trajectory.h5") as h,
        h5py.File(ROOT / "results/fast_pzt/f6/navigation_v1/zero_flow_feedback/trajectory.h5") as ref,
    ):
        x = h["positions_m"][:]
        old = ref["positions_m"][: len(x)]
        diff = float(np.max(np.linalg.norm(x - old, axis=-1)))
        cent = float(np.max(np.linalg.norm(x.mean(axis=2) - old.mean(axis=2), axis=-1)))
    # A stronger reproduction threshold than scientific dt tolerance; preserve discrepancy.
    good = all(s["numerical_quality"] == "PASS" for s in rows.values()) and cent < 1e-9
    write(
        out / "REPLAY.json",
        dict(
            status="PASS" if good else "FAIL",
            runs=rows,
            process_wall_s=time.perf_counter() - start,
            M0_F6_max_same_ID_difference_m=diff,
            M0_F6_max_centroid_difference_m=cent,
            frozen_control_source="F6 zero_flow_feedback actual events; no controller invoked",
        ),
    )


if __name__ == "__main__":
    main()
