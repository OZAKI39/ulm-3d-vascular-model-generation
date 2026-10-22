"""Two preregistered patches; paired wall OFF/ON acoustics at identical inputs."""

import argparse
import numpy as np
from .common import BASE, ROOT, js, write, identity
from .experiment import run_case
from mevo_mb.fast_sim.f4_engine import GPUVessel


def main(label, checks):
    assert js(BASE / checks / "RESULT.json")["status"] == "PASS"
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    protocol = js(ROOT / "configs/fast_pzt/f4/protocol.json")
    controls = js(ROOT / "configs/fast_pzt/f4/controls.json")["regions"]
    g = GPUVessel(BASE / "inputs/prepared_geometry.npz", BASE / "inputs/labelled_lumen_m.vtp")
    rows = {}
    for region in ["parent", "junction"]:
        x = np.load(BASE / "inputs" / (region + "_initial.npz"))["positions_m"]
        for enabled in [False, True]:
            name = region + ("_wall_on" if enabled else "_wall_off")
            cfg = dict(
                positions_m=x.tolist(),
                duration_s=protocol["aggregation_duration_s"],
                dt_s=protocol["dt_s"],
                sample_interval_s=protocol["sample_interval_s"],
                wall_budget_s=protocol["wall_budget_s"],
                wall_enabled=enabled,
                basis=True,
                flow_enabled=False,
                events=[controls[region]["baseline"]["event"]],
                initial_region=region,
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
    write(
        out / "RESULT.json",
        dict(
            status="COMPLETE",
            identity=identity(),
            rows=rows,
            scope="Actual original Y, actual PZT, zero flow; wall OFF retains original F3 physics",
        ),
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--checks", required=True)
    a = ap.parse_args()
    main(a.label, a.checks)
