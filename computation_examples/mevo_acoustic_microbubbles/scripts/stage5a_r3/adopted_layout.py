"""Reassess layout at the actual accepted half-step horizon, not only 50 ns."""

import subprocess
from time import perf_counter
import numpy as np
from scripts.stage5a_r3.common import BASE, write, identity
from scripts.stage5a_r3.layout import benchmark


def main():
    # A stopped scheduler alone is insufficient: the current science child must
    # have safely completed before this diagnostic competes for resources.
    commands = subprocess.check_output(["ps", "-eo", "args"], text=True).splitlines()
    assert not any("-m scripts.stage5a_r3.run " in c for c in commands), "WAIT_FOR_ACTIVE_TRAJECTORY"
    begin = perf_counter()
    rows = []
    layouts = ["global", "serial", "parallel2", "parallel4"]
    windows = [
        ("B2", "initial", 5e-9),
        ("B2", "stable", 2.5e-6),
        ("B2", "event", 2.5e-6),
        ("B3", "late", 2.5e-6),
    ]
    comparisons = []
    for case, window, dt in windows:
        for layout in layouts:
            rows.extend(benchmark((case, window, layout, True, "adopted_horizon", dt)))
        paths = [
            BASE / "benchmarks/layout_v2/adopted_horizon" / f"{case}_{window}_{layout}_optimized.npz"
            for layout in layouts
        ]
        arrays = [np.load(p) for p in paths]
        for layout, a in zip(layouts[1:], arrays[1:]):
            b = arrays[0]
            same = (
                np.array_equal(a["particle_id"], b["particle_id"])
                and a["positions_m"].shape == b["positions_m"].shape
            )
            err = float(np.linalg.norm(a["positions_m"] - b["positions_m"], axis=2).max()) if same else None
            comparisons.append(
                dict(
                    case=case,
                    window=window,
                    dt_s=dt,
                    layout=layout,
                    same_IDs_and_steps=same,
                    max_position_difference_m=err,
                    status="PASS" if same and err <= 1e-11 else "FAIL",
                    comparison_scope="Actual sequential substeps, global versus named layout; unchanged 1e-11 m numerical resolution check",
                )
            )
    result = dict(
        rows=rows,
        comparisons=comparisons,
        wall_s=perf_counter() - begin,
        identity=identity(),
        reason="Formal B3 demonstrated that 50 ns layout tests did not predict adopted 2.5 us half-step cost",
        physical_model_unchanged=True,
        status="PASS" if all(r["status"] == "PASS" for r in rows + comparisons) else "FAIL",
    )
    write(BASE / "benchmarks/LAYOUT_ADOPTED_HORIZON.json", result)
    print(dict(status=result["status"], wall_s=result["wall_s"], comparisons=comparisons), flush=True)


if __name__ == "__main__":
    main()
