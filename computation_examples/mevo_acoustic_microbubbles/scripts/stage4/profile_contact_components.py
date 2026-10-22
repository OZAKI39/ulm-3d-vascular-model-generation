"""Diagnostic only: compare component-wise and global contact NNLS solves.

Disconnected contact graphs have independent mobility-metric projections. This
trial retains the original global scaling and the same NNLS/feasibility gates.
It does not install the alternative into production or relax a physical guard.
"""
import argparse
import cProfile
import io
import json
import pstats
import time

import h5py
import numpy as np
from scipy.optimize import nnls
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from mevo_mb.bubbles.contact import contact_candidates, constrain_pair_velocities
from mevo_mb.bubbles import tracker
from scripts.stage4.common import ROOT, new_run, require, finish
from scripts.stage4.experiments import setup, run_case


def component_projection(x, radii, velocity, tolerance_m, max_iterations,
                         candidate_backend="brute", timings=None):
    pairs = contact_candidates(x, radii, tolerance_m, candidate_backend, timings)
    if not len(pairs):
        return velocity.copy()
    matrix = np.zeros((len(pairs), 3*len(x)))
    for row, (i, j) in enumerate(pairs):
        normal = (x[j]-x[i])/np.linalg.norm(x[j]-x[i])
        matrix[row, 3*i:3*i+3] = -normal
        matrix[row, 3*j:3*j+3] = normal
    sqrt_r = np.repeat(np.sqrt(radii/np.median(radii)), 3)
    scale = max(float(np.linalg.norm(velocity, axis=1).max()), np.finfo(float).tiny)
    free = velocity.ravel()/scale
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                       shape=(len(x), len(x))).tocsr()
    _, labels = connected_components(graph, directed=False)
    edge_labels = labels[pairs[:, 0]]
    multiplier = np.zeros(len(pairs))
    for component in np.unique(edge_labels):
        edges = np.flatnonzero(edge_labels == component)
        particles = np.flatnonzero(labels == component)
        columns = (3*particles[:, None]+np.arange(3)).ravel()
        block = matrix[np.ix_(edges, columns)]
        multiplier[edges], _ = nnls((block/sqrt_r[columns]).T,
                                   -sqrt_r[columns]*free[columns],
                                   maxiter=max_iterations*max(1, len(pairs)))
    constrained = (free+matrix.T@multiplier/(sqrt_r*sqrt_r))*scale
    if (matrix@constrained).min() < -max(1e-12, 1e-10*scale):
        raise RuntimeError("Component contact projection violates normal feasibility")
    return constrained.reshape(velocity.shape)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shadow-global", action="store_true",
                        help="Also check the original global solve at every candidate query")
    args = parser.parse_args()
    baseline = require("performance/profile_B2_secondary")
    baseline_folder = ROOT/baseline["row"]["run"]
    checks = []
    with h5py.File(baseline_folder/"particles/trajectories.h5") as h:
        assert h.attrs["complete"]
        for frame, t in enumerate(h["time_s"][:]):
            x, r, v = h["positions_m"][frame], h["radii_m"][:], h["free_total_velocity_m_s"][frame]
            values, timings = [], []
            for method in [constrain_pair_velocities, component_projection]:
                durations = []
                for _ in range(20):
                    start = time.perf_counter()
                    value = method(x, r, v, 1e-11, 3000, "ckdtree")
                    durations.append(time.perf_counter()-start)
                values.append(value)
                timings.append(float(np.median(durations)))
            error = float(np.linalg.norm(values[0]-values[1])/max(np.linalg.norm(values[0]), 1e-30))
            checks.append(dict(time_s=float(t), global_component_median_s=timings,
                               maximum_velocity_difference_m_s=float(np.max(abs(values[0]-values[1]))),
                               relative_L2=error, bitwise_equal=bool(np.array_equal(*values))))
            if error > 1e-12:
                raise RuntimeError("Component projection differs beyond the diagnostic roundoff gate")
    out = new_run("performance/component_contact_shadow" if args.shadow_global else "performance/component_contact_trial",
                  dict(secondary=True, duration_s=float(t), projection_checks=checks,
                       shadow_global=args.shadow_global, production_deployed=False), ["performance/profile_B2_secondary"])
    provided = setup()
    shadow = dict(queries=0, maximum_absolute_velocity_difference_m_s=0.,
                  maximum_relative_L2=0., queries_above_relative_1e_minus12=0)

    def checked_projection(*arguments, **keywords):
        candidate = component_projection(*arguments, **keywords)
        reference = constrain_pair_velocities(*arguments, **keywords)
        error = float(np.linalg.norm(candidate-reference)/max(np.linalg.norm(reference), 1e-30))
        shadow["queries"] += 1
        shadow["maximum_relative_L2"] = max(shadow["maximum_relative_L2"], error)
        shadow["maximum_absolute_velocity_difference_m_s"] = max(
            shadow["maximum_absolute_velocity_difference_m_s"], float(np.max(abs(candidate-reference))))
        shadow["queries_above_relative_1e_minus12"] += int(error > 1e-12)
        return candidate

    tracker.constrain_pair_velocities = checked_projection if args.shadow_global else component_projection
    profiler = cProfile.Profile()
    profiler.enable()
    row = run_case(out, "B2", provided[1]["particles"]["seed"], secondary=True,
                   duration=float(t), provided=provided)
    profiler.disable()
    profiler.dump_stats(str(out/"tracker.prof"))
    stream = io.StringIO()
    pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(35)
    (out/"PROFILE.txt").write_text(stream.getvalue())
    candidate_folder = ROOT/row["run"]
    equality, errors = {}, {}
    with h5py.File(baseline_folder/"particles/trajectories.h5") as a, h5py.File(candidate_folder/"particles/trajectories.h5") as b:
        for key in a:
            av, bv = a[key][:], b[key][:]
            equality[key] = bool(np.array_equal(av, bv, equal_nan=True))
            if av.shape == bv.shape and av.dtype.kind in "fc":
                defined = np.isfinite(av) & np.isfinite(bv)
                errors[key] = float(np.max(abs(av[defined]-bv[defined]))) if defined.any() else None
    events_equal = (baseline_folder/"particles/events.jsonl").read_bytes() == (candidate_folder/"particles/events.jsonl").read_bytes()
    validity_equal = json.loads((baseline_folder/"physics_validity.json").read_text()) == json.loads((candidate_folder/"physics_validity.json").read_text())
    equivalent = all(equality.values()) and events_equal and validity_equal
    finish(out, dict(status="PASS" if equivalent else "WARNING", baseline=baseline["row"], candidate=row,
                     projection_checks=checks, datasets_bitwise_equal=equality,
                     dataset_maximum_absolute_differences=errors, events_bitwise_equal=events_equal,
                     validity_equal=validity_equal, production_deployed=False,
                     shadow_global=shadow if args.shadow_global else None,
                     scope="Short instrumented diagnostic only. Mathematically separable projections can differ in floating-point arithmetic. No long-trajectory equivalence or production acceptance follows from this trial."))


if __name__ == "__main__":
    main()
