"""Three matched fixed-state pairs only after independent numerical admission."""

import argparse
import json
from time import perf_counter
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f10i_gpu import FixedStateB2
from mevo_mb.fast_sim.f10j_gpu import CorrectedB2, OptimizedB2
from mevo_mb.fast_sim.f10h_contact import geometric_contact
from mevo_mb.fast_sim.f10g_contact import resistance, normals
from .common import ROOT, BASE, js, write, identity, clean, sha


def common_output(result, A, x, free, ref):
    start = perf_counter()
    r, dt, kc = 1.78e-6, 1e-5, 0.04641319518378658
    zeta = 6 * np.pi * 0.001 * r
    y = np.asarray(result["positions_m"])
    f = np.asarray(result["contact_force_N"])
    i, j, d, n = normals(y)
    delta = np.maximum(2 * r - d, 0)
    actual = np.zeros_like(y)
    np.add.at(actual, i, -kc * delta[:, None] * n)
    np.add.at(actual, j, kc * delta[:, None] * n)
    res = y - x - dt * np.linalg.solve(A, free.ravel() + actual.ravel() / zeta).reshape(-1, 3)
    metrics = dict(
        equation_residual_m=float(np.linalg.norm(res, axis=1).max()),
        position_vs_independent_m=float(np.linalg.norm(y - ref, axis=1).max()),
        force_geometry_mismatch_N=float(np.max(abs(f - actual))),
        net_force_N=float(np.linalg.norm(f.sum(0))),
        overlap_m=float(delta.max()),
        single_correction_m=result["max_single_iteration_m"],
    )
    passed = all(
        [
            metrics["equation_residual_m"] <= 1e-14,
            metrics["position_vs_independent_m"] <= 1e-14,
            metrics["overlap_m"] <= 17.8e-9,
            metrics["single_correction_m"] <= 178e-9,
            metrics["net_force_N"] <= 1e-20,
            metrics["force_geometry_mismatch_N"] < 1e-18,
        ]
    )
    quality = perf_counter() - start
    start = perf_counter()
    payload = json.dumps(
        clean(
            dict(
                positions_m=y,
                contact_force_N=f,
                active_pairs=np.column_stack([i[delta > 0], j[delta > 0]]),
                metrics=metrics,
                converged=result["converged"],
            )
        ),
        sort_keys=True,
        allow_nan=False,
    ).encode()
    output = perf_counter() - start
    return dict(
        metrics=metrics,
        quality_pass=passed,
        quality_wall_s=quality,
        matched_output_commit_wall_s=output,
        matched_output_bytes=len(payload),
    )


def profile(g, x, free):
    """One explicitly instrumented diagnostic, separate from formal paired timings."""
    graphs = []
    evaluations = []
    launch = wp.capture_launch
    evaluate = g.evaluate

    def captured(graph, *args, **kwargs):
        ev0, ev1 = wp.Event(enable_timing=True), wp.Event(enable_timing=True)
        wp.record_event(ev0)
        start = perf_counter()
        value = launch(graph, *args, **kwargs)
        enqueue = perf_counter() - start
        wp.record_event(ev1)
        graphs.append((ev0, ev1, enqueue, "A" if graph is g.graphA else "H"))
        return value

    def assessed(q):
        start = perf_counter()
        reads = g.scalar_reads
        v = evaluate(q)
        evaluations.append(
            dict(trial=q is g.trial, wall_s=perf_counter() - start, reads=g.scalar_reads - reads)
        )
        return v

    wp.capture_launch = captured
    g.evaluate = assessed
    try:
        result = g.solve(x, free)
    finally:
        wp.capture_launch = launch
        g.evaluate = evaluate
    return dict(
        result=result,
        pcg_graphs=[
            dict(kind=k, graph_cuda_span_s=wp.get_event_elapsed_time(a, b) / 1000, host_enqueue_s=t)
            for a, b, t, k in graphs
        ],
        evaluations=evaluations,
        interpretation="PCG spans enclose only captured graphs; evaluate wall includes pending-device waits. Categories overlap and must not be summed as independent kernel costs.",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    assert js(BASE / "corrected_v1/RESULTS.json")["numerical_admitted"]
    assert js(BASE / "optimized_v2/RESULTS.json")["numerical_admitted"]
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    write(
        out / "IDENTITY.json",
        dict(**identity(), performance_protocol_sha256=sha(ROOT / "configs/fast_pzt/f10j/performance.json")),
    )
    saved = js(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json")["rows"][-1]
    x, free = np.array(saved["initial_positions_m"]), np.array(saved["free_velocity_m_s"])
    ref = np.array(js(BASE / "corrected_v1/RESULTS.json")["rows"][-1]["independent_reference"]["positions_m"])
    r, dt, kc = 1.78e-6, 1e-5, 0.04641319518378658
    Aref = resistance(x, r)["A"]
    wp.init()
    solvers = {}
    construction = []
    for name, cls in [("original", FixedStateB2), ("corrected", CorrectedB2), ("optimized", OptimizedB2)]:
        t = perf_counter()
        solvers[name] = cls(32)
        g = solvers[name]
        construction.append(
            dict(
                name=name,
                whole_construction_s=perf_counter() - t,
                allocation_s=g.allocation_s,
                original_graph_build_s=g.graph_build_s,
                corrected_graph_build_s=getattr(g, "correction_graph_build_s", 0.0),
            )
        )

    def run(name):
        start = perf_counter()
        if name == "cpu":
            A = resistance(x, r)["A"]
            assembly = perf_counter() - start
            result = geometric_contact(x, free, r, 6 * np.pi * 0.001 * r, dt, kc, A, atol_m=1e-15)
            bound = float(
                np.linalg.norm(
                    A @ result["contact_displacement_m"].ravel()
                    - dt * result["contact_force_N"].ravel() / (6 * np.pi * 0.001 * r)
                )
            )
            assert result["converged"] and bound <= 1e-14
            timing = dict(cpu_assembly_s=assembly, cpu_contact_s=result["total_s"], gradient_bound_m=bound)
        else:
            result = solvers[name].solve(x, free)
            A = Aref
            timing = dict(result["timing"])
        solve = perf_counter() - start
        common = common_output(result, A, x, free, ref)
        total = perf_counter() - start
        if name != "original":
            assert result["converged"] and common["quality_pass"]
        return dict(
            name=name,
            solver_converged=result["converged"],
            solver_reason=result.get("reason"),
            solve_and_preparation_s=solve,
            matched_full_wall_s=total,
            components=timing,
            **common,
        )

    # Exactly one warmup per implementation, excluded and explicitly stored.
    names = ["original", "corrected", "optimized", "cpu"]
    warm = [run(n) for n in names]
    samples = []
    for repeat in range(3):
        order = names[repeat:] + names[:repeat]
        samples.append(dict(repeat=repeat, order=order, measurements=[run(n) for n in order]))
    diagnostics = {name: profile(g, x, free) for name, g in solvers.items()}
    write(
        out / "RESULTS.json",
        dict(
            warmup=warm,
            samples=samples,
            construction=construction,
            diagnostics=diagnostics,
            original_failed_not_production_baseline=True,
            production=False,
            training_ready=False,
            fairness="Same state, equation, input preparation, original quality budgets, independent reference and matched serialized output. CPU original implementation atol=1e-15 additionally meets the unchanged GPU gradient bound 1e-14 m; all numerical targets verified. No disk output in either per-solve timer. Resident reused GPU graphs; construction and independent process costs separate.",
        ),
    )
    print(
        {
            name: float(
                np.median(
                    [
                        next(v for v in s["measurements"] if v["name"] == name)["matched_full_wall_s"]
                        for s in samples
                    ]
                )
            )
            for name in names
        },
        flush=True,
    )


if __name__ == "__main__":
    main()
