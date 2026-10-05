"""Graph changes and censored lifetimes from every accepted dense-probe state."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np

from scripts.stage4l.analyze_dense import edges, largest
from scripts.stage4l.dense_probe import write, sha


def graph_statistics(frames, times, cutoff):
    """Streaming lifetime accounting; no interpolation or fixed-clock resampling."""
    times = np.asarray(times, float)
    if len(times) < 2 or np.any(np.diff(times) <= 0):
        raise ValueError("At least two strictly increasing accepted times required")
    previous = edges(frames[0], cutoff)
    live = {edge: times[0] for edge in previous}
    initial = set(previous)
    changes, lifetimes, left, duration_sum = 0, 0, 0, 0.0
    maximum_turnover = 0.0
    for frame, time in zip(frames[1:], times[1:]):
        graph = edges(frame, cutoff)
        removed, added = previous - graph, graph - previous
        changes += len(removed) + len(added)
        maximum_turnover = max(maximum_turnover, (len(removed) + len(added)) / max(1, len(previous | graph)))
        for edge in removed:
            duration_sum += time - live.pop(edge)
            lifetimes += 1
            left += int(edge in initial)
            initial.discard(edge)
        for edge in added:
            live[edge] = time
        previous = graph
    for edge, start in live.items():
        duration_sum += times[-1] - start
        lifetimes += 1
        left += int(edge in initial)
    return dict(
        accepted_intervals=len(times) - 1,
        edge_changes=changes,
        edge_changes_per_s=changes / (times[-1] - times[0]),
        maximum_single_interval_turnover_fraction=maximum_turnover,
        observed_lifetimes=lifetimes,
        observed_lifetime_mean_s=duration_sum / lifetimes if lifetimes else None,
        left_censored_lifetimes=left,
        right_censored_lifetimes=len(live),
        final_edges=len(previous),
        observation_duration_s=float(times[-1] - times[0]),
    )


def hard_threshold_motion(frames, radius, tolerance=1e-11):
    """Magnitude of real gap changes at hard-contact membership switches."""
    previous = edges(frames[0], 2 * radius + tolerance)
    gap_changes, endpoint_gaps, relative_motion = [], [], []
    for old, new in zip(frames[:-1], frames[1:]):
        current = edges(new, 2 * radius + tolerance)
        changed = np.asarray(list(previous ^ current), dtype=int).reshape(-1, 2)
        if len(changed):
            i, j = changed.T
            before, after = old[j] - old[i], new[j] - new[i]
            old_gap = np.linalg.norm(before, axis=1) - 2 * radius
            new_gap = np.linalg.norm(after, axis=1) - 2 * radius
            gap_changes.append(abs(new_gap - old_gap))
            endpoint_gaps.append(np.maximum(abs(old_gap), abs(new_gap)))
            relative_motion.append(np.linalg.norm(after - before, axis=1))
        previous = current
    quantiles = [0, 0.5, 0.9, 0.99, 1]
    return dict(
        hard_contact_tolerance_m=tolerance,
        switched_edge_events=sum(map(len, gap_changes)),
        quantiles=quantiles,
        absolute_gap_change_m=np.quantile(np.concatenate(gap_changes), quantiles).tolist()
        if gap_changes
        else [],
        endpoint_absolute_gap_m=np.quantile(np.concatenate(endpoint_gaps), quantiles).tolist()
        if endpoint_gaps
        else [],
        relative_pair_vector_change_m=np.quantile(np.concatenate(relative_motion), quantiles).tolist()
        if relative_motion
        else [],
        scope="Descriptive diagnostic only. No hysteresis, coordinate filtering, tolerance change or modified contact handling.",
    )


def analyze(run):
    summary = json.loads((run / "SUMMARY.json").read_text())
    if summary["status"] != "COMPLETE":
        raise ValueError("Native graph comparison requires completed probes")
    with h5py.File(run / "trajectories.h5") as h:
        if not h.attrs["all_accepted_substeps_recorded"]:
            raise ValueError("Probe does not contain all accepted states")
        t = h["time_s"][:]
        r = float(h["radii_m"][0])
        active = h["active"][0]
        indices = np.flatnonzero(active)[largest(h["positions_m"][0, active], 2.1 * r)]
        x = h["positions_m"][:][:, indices]
        ids = h["particle_id"][:][indices]
    result = dict(
        case=summary["case"],
        backend=summary["backend"],
        same_ID_particle_ids=ids.tolist(),
        all_accepted_states=True,
        near_contact=graph_statistics(x, t, 2.1 * r),
        hard_contact=graph_statistics(x, t, 2 * r + 1e-11),
        hard_threshold_motion=hard_threshold_motion(x, r),
        script_sha256=sha(Path(__file__)),
        semantics="Actual accepted states; no position interpolation or clock resampling. Edge lifetimes are censored by the observation window, not material persistence estimates.",
    )
    write(run / "NATIVE_GRAPH_METRICS.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("runs", type=Path, nargs="+")
    args = parser.parse_args()
    for run in args.runs:
        result = analyze(run)
        print(json.dumps(dict(run=str(run), **result)), flush=True)


if __name__ == "__main__":
    main()
