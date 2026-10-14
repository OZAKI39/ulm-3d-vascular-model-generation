"""Freeze at most two full-window controls from actual successful short runs."""

import numpy as np
import h5py
from .common import BASE, read, write, sha, identity


def main():
    rows = []
    controls = read(BASE / "CANDIDATES.json")["candidates"]
    for c in controls:
        folder = BASE / "b_screen_v1/runs" / c["name"]
        r = read(folder / "RESULT.json")
        with h5py.File(folder / "trajectory.h5") as h:
            pressure = float(h["extra_metrics"][:, 4].max())
        displacement = np.array(r["displacement_m"])
        target = np.array(c["target_m"])
        direction = float(displacement @ target / np.linalg.norm(target))
        # Scaling observed short displacement is a ranking heuristic, not a new trajectory.
        score = float(
            np.sum((displacement * 4 - target) ** 2) / np.dot(target, target) + 0.02 * (pressure / 1800) ** 2
        )
        rows.append(
            dict(
                name=c["name"],
                axis=c["axis"],
                score=score,
                short_result=r,
                short_sha256=sha(folder / "trajectory.h5"),
                measured_pressure_max_Pa=pressure,
                forward_displacement_m=direction,
                eligible=bool(r["status"] == "PASS" and direction > 0),
                control=c["control"],
            )
        )
    selected = [
        min([r for r in rows if r["axis"] == a and r["eligible"]], key=lambda r: r["score"])
        for a in ["X", "Y"]
        if any(r["axis"] == a and r["eligible"] for r in rows)
    ]
    write(
        BASE / "NAVIGATION_SELECTION.json",
        dict(
            identity=identity(),
            rule="Complete 0.5ms, all original gates and forward motion. Rank squared target error from observed displacement times four plus 0.02*(observed max pressure/1800Pa)^2 as pressure-risk cost. Heuristic only, no success inferred. One X and one Y; fresh full integration and independent dt/cache reruns.",
            selected=selected,
            all_rankings=rows,
        ),
    )
    print([(r["name"], r["score"], r["short_result"]["max_local_qps"]) for r in rows])
    print("selected", [r["name"] for r in selected])


if __name__ == "__main__":
    main()
