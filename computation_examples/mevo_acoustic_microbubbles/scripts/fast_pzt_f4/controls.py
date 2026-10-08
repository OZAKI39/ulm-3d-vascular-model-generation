"""Matched wall/control factorial continuation from one actual near-wall cluster."""

import argparse
import hashlib
import numpy as np
import h5py
from mevo_mb.fast_sim.f4_engine import GPUVessel
from scripts.fast_pzt_f3.repeatability import pair_comparison
from scripts.fast_pzt_f2.analyze import read
from .common import BASE, ROOT, js, write, identity
from .experiment import run_case


def main(label, aggregation):
    formation = js(BASE / aggregation / "RESULT.json")
    assert formation["status"] == "COMPLETE"
    region = "parent" if formation["rows"]["parent_wall_on"]["final_largest_size"][0] >= 16 else "junction"
    source = BASE / aggregation / (region + "_wall_on") / "trajectory.h5"
    with h5py.File(source) as h:
        x = h["positions_m"][-1]
        offset = float(h["time_s"][-1])
        assert h["active"][-1].all()
    protocol = js(ROOT / "configs/fast_pzt/f4/protocol.json")
    events = js(ROOT / "configs/fast_pzt/f4/controls.json")["regions"][region]
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(
        out / "CONTINUATION_STATE.npz",
        positions_m=x,
        particle_id=np.arange(x.shape[1]),
        absolute_time_s=offset,
    )
    write(
        out / "STATE_IDENTITY.json",
        dict(
            source=str(source.relative_to(BASE)),
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            absolute_time_s=offset,
            region=region,
            selection="Predeclared parent if largest >=16; otherwise junction. All 32 IDs retained; no new maturity claim",
        ),
    )
    g = GPUVessel(BASE / "inputs/prepared_geometry.npz", BASE / "inputs/labelled_lumen_m.vtp")
    rows = {}
    for enabled in [True, False]:
        for control in ["baseline", "switched"]:
            name = ("wall_on_" if enabled else "wall_off_") + control
            cfg = dict(
                positions_m=x.tolist(),
                duration_s=protocol["continuation_duration_s"],
                dt_s=protocol["dt_s"],
                sample_interval_s=protocol["sample_interval_s"],
                wall_budget_s=protocol["wall_budget_s"],
                wall_enabled=enabled,
                basis=True,
                flow_enabled=False,
                events=[events[control]["event"]],
                initial_region=region,
                absolute_time_offset_s=offset,
                execution_mode="stable_reduction",
            )
            rows[name] = run_case(out / name, cfg, g)
            write(out / "PARTIAL.json", dict(identity=identity(), rows=rows))
            if rows[name]["numerical_quality"] != "PASS":
                write(
                    out / "RESULT.json",
                    dict(status="FAILED_REPRESENTATIVE_REQUIRES_DIAGNOSIS", identity=identity(), rows=rows),
                )
                return
    contrasts = {}
    for a, b in [
        ("wall_on_baseline", "wall_on_switched"),
        ("wall_off_baseline", "wall_off_switched"),
        ("wall_off_baseline", "wall_on_baseline"),
        ("wall_off_switched", "wall_on_switched"),
    ]:
        d = read(out / a)
        e = read(out / b)
        np.testing.assert_array_equal(d["positions_m"][0], e["positions_m"][0])
        c = pair_comparison(d, e)
        name = b + "_vs_" + a
        np.savez_compressed(out / (name + ".npz"), time_s=d["time_s"], **c)
        contrasts[name] = {k: np.asarray(v[-1]).tolist() for k, v in c.items()}
    write(
        out / "RESULT.json",
        dict(
            status="COMPLETE",
            identity=identity(),
            rows=rows,
            region=region,
            endpoint_contrasts=contrasts,
            control="Y+ phase +30 degrees; amplitude/enable unchanged",
            interpretation="Matched 2x2 control/wall comparison; not an independent physical calibration or robustness proof. Inherited F3 sensitivity remains.",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--aggregation", required=True)
    a = ap.parse_args()
    main(a.label, a.aggregation)
