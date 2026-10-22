"""Same-array regression, real release events, and independent limiting cases."""

import argparse
from time import perf_counter
import h5py
import numpy as np
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.f10d_metrics import process, components, topology_fast
from mevo_mb.fast_sim.morphodynamics import topology, rigid_fit
from .postprocess_reference import process as reference
from .common import BASE, ROOT, F10C, F5, identity, js, sha, write, clean


def arrays_equal(a, b, key):
    if a.shape != b.shape:
        return dict(status="FAIL", key=key, shapes=[a.shape, b.shape])
    exact = np.array_equal(a, b, equal_nan=True)
    if a.dtype.kind in "biu":
        return dict(
            key=key,
            status="PASS" if exact else "FAIL",
            bitwise=exact,
            max_absolute_error=0 if exact else None,
        )
    category = (
        "tensor_m2"
        if key.endswith("_m2")
        else "length_m"
        if key.endswith("_m")
        else "angle_deg"
        if key.endswith("_deg")
        else "dimensionless"
    )
    gate = js(ROOT / "configs/fast_pzt/f10d/protocol.json")["float_gates"][category]
    good = np.allclose(a, b, **gate, equal_nan=True)
    finite = np.isfinite(a) & np.isfinite(b)
    error = float(np.max(np.abs(a[finite] - b[finite]))) if finite.any() else 0.0
    return dict(
        key=key,
        status="PASS" if good else "FAIL",
        bitwise=bool(exact),
        max_absolute_error=error,
        category=category,
        gate=gate,
    )


def compare_results(a, b):
    checks = []
    for index, group in [(0, "metrics"), (2, "morphology")]:
        assert set(a[index]) == set(b[index])
        for key in a[index]:
            checks.append(arrays_equal(a[index][key], b[index][key], group + "/" + key))
    for index, key in [(1, "sampled_lineage"), (3, "accepted_step_topology")]:
        checks.append(
            dict(
                key=key,
                status="PASS" if clean(a[index]) == clean(b[index]) else "FAIL",
                bitwise=clean(a[index]) == clean(b[index]),
            )
        )
    return dict(status="PASS" if all(c["status"] == "PASS" for c in checks) else "FAIL", checks=checks)


def independent_checks():
    rng = np.random.default_rng(570)
    graphs = rng.uniform(size=(40, 32, 32)) < 0.04
    graphs |= graphs.swapaxes(1, 2)
    graphs[0] = False
    graphs[1] = True
    labels = components(graphs)
    for a, b in zip(graphs, labels):
        np.testing.assert_array_equal(connected_components(a, directed=False)[1], b)
    # Tetrahedral repeated eigenvalues, linear rank deficiency, zero cloud,
    # single point and an ordinary proper rotation with fixed IDs.
    clouds = [
        np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) * 1e-6,
        np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0], [3, 0, 0]]) * 1e-6,
        np.zeros((4, 3)),
        np.zeros((1, 3)),
        rng.normal(size=(8, 3)) * 1e-6,
    ]
    angle = 0.3
    rot = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    cases = []
    for i, x in enumerate(clouds):
        xx = np.stack([x, x @ rot + np.array([1e-5, -2e-6, 1e-6])])[:, None]
        # These are analytic software checks, not physical trajectories.
        # Avoid topology because cloud separation/contact may change under translation rounding.
        from mevo_mb.fast_sim.metrics import trajectory_metrics as old_metrics
        from mevo_mb.fast_sim.morphodynamics import morphology_series
        from mevo_mb.fast_sim.f10d_metrics import trajectory_metrics, morphology_batch

        with np.errstate(divide="ignore", invalid="ignore"):
            ma, la = old_metrics(xx, 1.78e-6, [0, 1])
            mb, lb = trajectory_metrics(xx, 1.78e-6, [0, 1])
            mo = {
                k: v[:, None]
                for k, v in morphology_series(xx[:, 0], np.arange(len(x)), ma["labels"][:, 0]).items()
            }
            mn = morphology_batch(xx, mb["labels"], mb)
        report = compare_results((ma, la, mo, []), (mb, lb, mn, []))
        rotation, rmsd = rigid_fit(xx[0, 0], xx[1, 0])
        assert rmsd < 1e-18 and np.linalg.det(rotation) > 0
        cases.append(dict(case=i, **report))
    # Explicit original persistence boundary: split is right-censored; then merge.
    x = np.array([[0.0, 0, 0], [3.6e-6, 0, 0]])
    ft = np.array([0.0, 0.001, 0.004])
    ei = np.array([[0, 0, 1, 0], [0, 0, 1, 1]])
    et = np.array([0.001, 0.003])
    fl = np.array([[0, 0], [0, 1], [0, 0]])
    old = topology(x, ft, fl, et, ei, 0.004)
    new = topology_fast(x, ft, fl, et, ei, 0.004)
    assert old[0] == new[0]
    np.testing.assert_array_equal(old[1], new[1])
    np.testing.assert_array_equal(old[2], new[2])
    try:
        from mevo_mb.fast_sim.f10d_metrics import trajectory_metrics

        trajectory_metrics(np.empty((2, 1, 0, 3)), 1.78e-6, [0, 1])
    except ValueError:
        pass
    else:
        raise AssertionError("Empty cloud must be explicitly rejected")
    return dict(
        status="PASS" if all(r["status"] == "PASS" for r in cases) else "FAIL",
        cases=cases,
        graph_cases=40,
        persistence=old[0],
        scope="Analytic software checks only; not manufactured scientific split/merge",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    checks = independent_checks()
    write(out / "INDEPENDENT_CHECKS.json", checks)
    rows = []
    for name, path in [
        ("B64", F10C / "exclusive_v1_batch/B64_batch/trajectory.h5"),
        ("F5_release", F5 / "science_v1/main/trajectory.h5"),
    ]:
        started = perf_counter()
        with h5py.File(path) as h:
            x = h["positions_m"][:]
            times = h["time_s"][:]
            radius = float(h["radii_m"][0])
            ei = h["edge_events/batch_i_j_linked"][:]
            et = h["edge_events/time_s"][:]
        read_s = perf_counter() - started
        before = sha(path)
        for repeat in ["first_after_import", "warm"] if name == "B64" else ["once"]:
            old = reference(x, radius, times, ei, et)
            new = process(x, radius, times, ei, et)
            report = compare_results(old, new)
            row = dict(
                name=name,
                repeat=repeat,
                input_sha256=before,
                read_s=read_s,
                reference_timing=old[4],
                optimized_timing=new[4],
                speedup=old[4]["total_s"] / new[4]["total_s"],
                **report,
            )
            rows.append(row)
            print(
                dict(
                    name=name,
                    repeat=repeat,
                    status=report["status"],
                    old_s=old[4]["total_s"],
                    new_s=new[4]["total_s"],
                ),
                flush=True,
            )
            if name == "F5_release":
                write(out / "F5_TOPOLOGY.json", new[3])
                np.savez_compressed(
                    out / "F5_RELEASE_METRICS.npz",
                    time_s=times,
                    positions_m=x,
                    labels=new[0]["labels"],
                    largest_size=new[0]["largest_size"],
                    rg_m=new[0]["global_rg_m"],
                    aspect=new[2]["fixed_aspect_ratio"],
                    edge_times_s=et,
                    edge_ids=ei,
                )
        assert before == sha(path)
    status = "PASS" if checks["status"] == "PASS" and all(r["status"] == "PASS" for r in rows) else "FAIL"
    write(out / "RESULTS.json", dict(status=status, rows=rows, training_ready=False))
    assert status == "PASS", "Retain failed regression; do not promote optimized path"


if __name__ == "__main__":
    main()
