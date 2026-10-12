"""Independent all-pair accepted-edge reconstruction and worst B2 timing audit."""

from collections import defaultdict
import h5py
import numpy as np
from .common import ROOT, BASE, js, write
from .analyze import read_steps


def main():
    summary = (
        js(BASE / "analysis/LUBRICATION_TRAJECTORIES.json")["trials"]
        + js(BASE / "analysis/B2_GEOMETRY_TRAJECTORIES.json")["trials"]
    )
    records = []
    series = {}
    for row in summary:
        path = ROOT / row["path"]
        times, x, _ = read_steps(path)
        i, j = np.triu_indices(32, 1)
        gaps = np.linalg.norm(x[:, j] - x[:, i], axis=-1) - 3.56e-6
        linked = gaps <= 0.178e-6
        changes = np.argwhere(linked[1:] != linked[:-1])
        expected = []
        for frame, pair in changes:
            expected.append(
                (round(times[frame + 1], 12), int(i[pair]), int(j[pair]), int(linked[frame + 1, pair]))
            )
        with h5py.File(path / "trajectory.h5") as h:
            actual = [
                (round(float(t), 12), int(a), int(b), int(on))
                for t, (_, a, b, on) in zip(h["edge_events/time_s"][:], h["edge_events/batch_i_j_linked"][:])
                if t <= times[-1] + 1e-14
            ]
        assert sorted(expected) == sorted(actual), (path, len(expected), len(actual))
        records.append(
            dict(
                path=row["path"],
                accepted_times=len(times),
                reconstructed_events=len(expected),
                exact_GPU_edge_reconstruction=True,
            )
        )
        if row["mode"].startswith("B2_GEOM"):
            series[row["mode"]] = dict(times=times, x=x, gaps=gaps, events=actual, i=i, j=j)
    a, b = series["B2_GEOM_K2_DT10"], series["B2_GEOM_K2_DT5"]
    d = []
    for s in [a, b]:
        v = defaultdict(list)
        for time, i, j, on in s["events"]:
            v[i, j].append((on, time))
        d.append(v)
    matches = []
    for pair in sorted(set(d[0]) | set(d[1])):
        assert [e[0] for e in d[0][pair]] == [e[0] for e in d[1][pair]]
        for u, v in zip(d[0][pair], d[1][pair]):
            matches.append(
                dict(pair=pair, on=u[0], time_DT10_s=u[1], time_DT5_s=v[1], difference_s=abs(u[1] - v[1]))
            )
    worst = max(matches, key=lambda e: e["difference_s"])
    i, j = worst["pair"]
    index = np.flatnonzero((a["i"] == i) & (a["j"] == j))[0]
    lo = max(0, min(worst["time_DT10_s"], worst["time_DT5_s"]) - 0.0003)
    hi = max(worst["time_DT10_s"], worst["time_DT5_s"]) + 0.0003
    window = []
    for mode, s in [("DT10", a), ("DT5", b)]:
        keep = (s["times"] >= lo) & (s["times"] <= hi)
        window.append(dict(mode=mode, time_s=s["times"][keep], surface_gap_m=s["gaps"][keep, index]))
    same = [(k, int(np.argmin(abs(b["times"] - t)))) for k, t in enumerate(a["times"])]
    difference = max(abs(a["gaps"][k, index] - b["gaps"][other, index]) for k, other in same)
    write(
        BASE / "analysis/EDGE_AUDIT.json",
        dict(
            status="PASS_EXACT_GPU_EVENT_RECONSTRUCTION",
            trials=records,
            B2_matched_events=matches,
            worst=worst,
            worst_pair_max_gap_difference_m=difference,
            worst_pair_window=window,
            edge_threshold_surface_gap_m=0.178e-6,
            interpretation="Connectivity threshold is d<=2.1R, distinct from physical contact at d<2R. Same on/off pair histories but event-time and cross-pair-order gates fail. No time interpolation or changed threshold.",
        ),
    )
    print(
        "Exact accepted-edge audit PASS",
        len(records),
        "trials; worst",
        worst,
        "max same-time gap difference",
        difference,
    )


if __name__ == "__main__":
    main()
