"""Single connected fixed-state benchmarks; process watchdog enforces all budgets.

No trajectories. Each N worker owns one CUDA context; CPU and GPU share identical
saved x/free/r/kc. Reused graphs are measured separately from complete solves.
"""

import argparse
import json
import resource
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
from scipy import linalg
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from threadpoolctl import threadpool_info
from .common import ROOT, BASE, js, sha, write, identity, clean

R = 1.78e-6
DT = 1e-5
KC = 0.04641319518378658
ZETA = 6 * np.pi * 0.001 * R


def estimate(n):
    pairs = n * (n - 1) // 2
    # Unchanged CPU contact_geometry explicitly allocates pairs x3N, including
    # distant pairs. A live current and trial B coexist during Armijo evaluation.
    return dict(
        N=n,
        CPU_peak_array_estimate_bytes=2 * pairs * 3 * n * 8
        + 4 * (3 * n) ** 2 * 8
        + 3 * n * n * 9 * 8
        + 64 * pairs,
        GPU_peak_array_estimate_bytes=192 * n * n + 4096 * n + 16 * 1024**2,
        formula="CPU2 simultaneous all-pair B +4 dense3Nx3N +3 pair3x3 blocks +pair scratch; GPU conservative192N^2+4096N+16MiB; excludes interpreter/CUDA context",
    )


def state(n):
    if n == 32:
        p = ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json"
        old = js(p)["rows"][-1]
        x = np.array(old["initial_positions_m"])
        free = np.array(old["free_velocity_m_s"])
        label = "REAL_ARCHIVED_F10J_FIXED_STATE"
        source = {str(p.relative_to(ROOT)): sha(p)}
    else:
        side = int(np.ceil(n ** (1 / 3)))
        ijk = np.indices((side, side, side)).reshape(3, -1).T[:n]
        x = ijk.astype(float) * (2 * R - 0.001 * R)
        x -= x.mean(0)
        free = -1.0 * x
        label = "BENCHMARK_ONLY_SYNTHETIC_STATE"
        source = {}
    d = squareform(pdist(x))
    np.fill_diagonal(d, np.inf)
    graph = d < 2 * R
    count, labels = connected_components(graph, directed=False)
    overlap = np.maximum(2 * R - d[np.triu_indices(n, 1)], 0)
    active = overlap[overlap > 0]
    assert np.isfinite(x + free).all() and active.size > 0 and overlap.max() <= 17.8e-9
    if n > 32:
        assert count == 1
    return (
        x,
        free,
        dict(
            label=label,
            source=source,
            N=n,
            B=1,
            components=count,
            contact_pairs=int(graph.sum() // 2),
            per_particle_contact_degree=graph.sum(1),
            initial_positive_overlap_quantiles_m=np.quantile(active, [0, 0.25, 0.5, 0.75, 1]),
            radius_m=R,
            dt_s=DT,
            kc_N_m=KC,
            free_force_N=ZETA * free,
            synthetic_construction=None
            if n == 32
            else "One contiguous nearly touching simple-cubic cluster; spacing1.999R; free=-1s^-1*(x-centroid), external force=zeta*free. Not CFD or acoustic force. Regular lattice is simpler than disordered scientific contacts.",
            training_ready=False,
        ),
    )


def quality(result, A, x, free):
    from mevo_mb.fast_sim.f10g_contact import normals

    y = np.asarray(result["positions_m"])
    f = np.asarray(result["contact_force_N"])
    i, j, d, n = normals(y)
    overlap = np.maximum(2 * R - d, 0)
    force = np.zeros_like(y)
    pair = KC * overlap[:, None] * n
    np.add.at(force, i, -pair)
    np.add.at(force, j, pair)
    v = linalg.solve(A, free.ravel() + force.ravel() / ZETA, assume_a="pos").reshape(-1, 3)
    eq = y - x - DT * v
    metric = dict(
        equation_residual_m=float(np.linalg.norm(eq, axis=1).max()),
        max_overlap_m=float(overlap.max()),
        single_correction_m=result["max_single_iteration_m"],
        force_geometry_mismatch_N=float(abs(f - force).max()),
        net_force_N=float(np.linalg.norm(force.sum(0))),
        nonnegative_scalar_reaction=bool(np.all(KC * overlap >= 0)),
        separated_pair_force_exact_zero=bool(np.all((KC * overlap)[d >= 2 * R] == 0)),
        active_contacts=int(np.count_nonzero(overlap)),
        gradient_bound_m=result.get("conservative_equation_bound_m", result.get("gradient_bound_m")),
    )
    passed = bool(
        result["converged"]
        and metric["equation_residual_m"] <= 1e-14
        and metric["max_overlap_m"] <= 17.8e-9
        and metric["single_correction_m"] <= 178e-9
        and metric["force_geometry_mismatch_N"] < 1e-18
        and metric["net_force_N"] < 1e-20
        and metric["gradient_bound_m"] <= 1e-14
        and np.isfinite(y).all()
    )
    # Exactly matched scientific output fields on both implementations; disk outside timer.
    payload = dict(
        positions_m=y,
        contact_force_N=f,
        active_pairs=np.c_[i[overlap > 0], j[overlap > 0]],
        metrics=metric,
        converged=passed,
    )
    return passed, metric, payload


def worker(out, n):
    from mevo_mb.fast_sim.f10g_contact import resistance
    from mevo_mb.fast_sim.f10h_contact import geometric_contact
    from mevo_mb.fast_sim.f10j_gpu import OptimizedB2
    from mevo_mb.fast_sim.f10k_scaling import ScalingB2
    import warp as wp

    start = time.perf_counter()
    out.mkdir(exist_ok=False)

    def phase(name, **kwargs):
        value = dict(phase=name, start_monotonic=time.monotonic(), **kwargs)
        temp = out / "PROGRESS.tmp"
        write(temp, value)
        temp.replace(out / "PROGRESS.json")

    phase("prepare")
    x, free, meta = state(n)
    np.savez_compressed(out / "INPUT.npz", positions_m=x, free_velocity_m_s=free, radius_m=R)
    write(out / "INPUT.json", dict(**meta, npz_sha256=sha(out / "INPUT.npz"), static_memory=estimate(n)))
    t = time.perf_counter()
    Aref = resistance(x, R)["A"]
    reference_setup = time.perf_counter() - t
    phase("construct")
    wp.init()
    t = time.perf_counter()
    g = ScalingB2(n)
    construction = time.perf_counter() - t
    # Default constructors must continue rejecting sizes outside their original domain.
    try:
        OptimizedB2(128)
    except ValueError:
        default_guard = True
    else:
        raise AssertionError("Original default limit changed")
    write(
        out / "CONSTRUCTION.json",
        dict(
            wall_s=construction,
            allocation_s=g.allocation_s,
            inherited_graph_s=g.graph_build_s,
            replacement_graph_s=g.correction_graph_build_s,
            explicit_array_capacity_bytes=sum(
                v.capacity for v in g.__dict__.values() if isinstance(v, wp.array)
            ),
            buffer_scope="Explicit class-owned arrays, excluding internal CG work/context; measured process memory separately",
            reference_setup_s=reference_setup,
            default_limit_guard=default_guard,
            actual_BLAS=threadpool_info(),
            static_memory=estimate(n),
        ),
    )
    assert all(v["num_threads"] == 1 for v in threadpool_info() if v["user_api"] == "blas")
    measurements = []
    last = {}
    sum_s = 0.0
    for repeat in [-1, 0, 1]:
        order = ["cpu", "gpu"] if repeat % 2 else ["gpu", "cpu"]
        for method in order:
            phase("solve", method=method, repeat=repeat, completed_solve_s=sum_s)
            t = time.perf_counter()
            if method == "cpu":
                A = resistance(x, R)["A"]
                prep = time.perf_counter() - t
                linear = []
                original = linalg.solve

                def audit_solve(a, b, *args, **kwargs):
                    sol = original(a, b, *args, **kwargs)
                    residual = float(np.linalg.norm(b - a @ sol))
                    target = max(1e-14, 1e-12 * float(np.linalg.norm(b)))
                    linear.append(dict(true_residual=residual, target=target, pass_=residual <= target))
                    return sol

                linalg.solve = audit_solve
                try:
                    result = geometric_contact(x, free, R, ZETA, DT, KC, A, atol_m=1e-15)
                finally:
                    linalg.solve = original
                bound = float(
                    np.linalg.norm(
                        A @ result["contact_displacement_m"].ravel()
                        - DT * result["contact_force_N"].ravel() / ZETA
                    )
                )
                rhs = free.ravel() * DT / R
                u = result["mobility_free_m_s"].ravel() * DT / R
                base_res = float(np.linalg.norm(rhs - A @ u))
                target = max(1e-14, 1e-12 * float(np.linalg.norm(rhs)))
                result["gradient_bound_m"] = bound
                result["converged"] = bool(
                    result["converged"]
                    and bound <= 1e-14
                    and base_res <= target
                    and all(v["pass_"] for v in linear)
                )
                detail = dict(
                    preparation_wall_s=prep,
                    solve_wall_s=result["total_s"],
                    nonlinear_updates=len(result["iterations"]),
                    base_true_residual=base_res,
                    base_target=target,
                    direction_linear_audit=linear,
                    PCG="NOT_APPLICABLE direct Cholesky; true residual uses same norm/targets",
                )
            else:
                result = g.solve(x, free)
                A = Aref
                detail = dict(
                    result["timing"],
                    nonlinear_updates=len(result["iterations"]),
                    linear_audit=result["linear_audit"],
                    linear_solves=result["linear_solves"],
                    solver_reason=result["reason"],
                )
            solve = time.perf_counter() - t
            tq = time.perf_counter()
            passed, metrics, payload = quality(result, A, x, free)
            qtime = time.perf_counter() - tq
            to = time.perf_counter()
            serialized = json.dumps(clean(payload), sort_keys=True, allow_nan=False).encode()
            outputtime = time.perf_counter() - to
            whole = time.perf_counter() - t
            sum_s += whole
            item = dict(
                method=method,
                repeat=repeat,
                kind="warmup" if repeat < 0 else "formal",
                converged=bool(result["converged"]),
                quality_pass=passed,
                metrics=metrics,
                full_wall_s=whole,
                preparation_solve_readback_s=solve,
                quality_s=qtime,
                matched_serialization_s=outputtime,
                matched_payload_bytes=len(serialized),
                components=detail,
            )
            measurements.append(item)
            last[method] = payload
            if method == "gpu":
                # Compare actual returned operator coefficients/normals to the same CPU A,
                # outside performance timers but inside bounded stage budget.
                coeff = result["resistance_coefficients"]
                norm = result["resistance_normals"]
                ag = np.eye(3 * n)
                for i, j in np.argwhere(np.triu(coeff > 0, 1)):
                    block = coeff[i, j] * np.outer(norm[i, j], norm[i, j])
                    s = slice(3 * i, 3 * i + 3)
                    v = slice(3 * j, 3 * j + 3)
                    ag[s, s] += block
                    ag[v, v] += block
                    ag[s, v] -= block
                    ag[v, s] -= block
                item["operator_max_abs_difference"] = float(np.max(abs(ag - Aref)))
                assert item["operator_max_abs_difference"] < 1e-8
            if n == 32 and method == "gpu":
                archived = js(ROOT / "results/fast_pzt/f10j/optimized_v2/RESULTS.json")["rows"][-1]["gpu"]
                item["F10J_archived_bitwise_positions"] = bool(
                    np.array_equal(result["positions_m"], archived["positions_m"])
                )
                item["F10J_archived_bitwise_forces"] = bool(
                    np.array_equal(result["contact_force_N"], archived["contact_force_N"])
                )
                assert item["F10J_archived_bitwise_positions"] and item["F10J_archived_bitwise_forces"]
            write(out / "MEASUREMENTS.json", measurements)
            np.savez_compressed(
                out / f"{method}_latest.npz",
                positions_m=payload["positions_m"],
                contact_force_N=payload["contact_force_N"],
                active_pairs=payload["active_pairs"],
            )
            phase("between", completed_solve_s=sum_s)
            if not passed or whole >= 8:
                write(
                    out / "STOP.json",
                    dict(
                        reason="NUMERICAL_FAILURE" if not passed else "SINGLE_SOLVE_NEAR10S",
                        method=method,
                        repeat=repeat,
                    ),
                )
                break
        if (out / "STOP.json").exists():
            break
    comparison = None
    if set(last) == {"cpu", "gpu"}:
        c, gpu = last["cpu"], last["gpu"]
        comparison = dict(
            position_max_difference_m=float(
                np.linalg.norm(c["positions_m"] - gpu["positions_m"], axis=1).max()
            ),
            force_max_difference_N=float(np.max(abs(c["contact_force_N"] - gpu["contact_force_N"]))),
            active_pairs_exact=bool(np.array_equal(c["active_pairs"], gpu["active_pairs"])),
        )
        comparison["pass"] = bool(
            comparison["position_max_difference_m"] <= 1e-14
            and comparison["force_max_difference_N"] <= 1e-15
            and comparison["active_pairs_exact"]
        )
    formal = {
        m: [
            v["full_wall_s"]
            for v in measurements
            if v["method"] == m and v["kind"] == "formal" and v["quality_pass"]
        ]
        for m in ["cpu", "gpu"]
    }
    admitted = bool(
        comparison
        and comparison["pass"]
        and all(v["quality_pass"] for v in measurements)
        and all(len(v) == 2 for v in formal.values())
    )
    write(
        out / "RESULTS.json",
        dict(
            N=n,
            B=1,
            status="PASS" if admitted else "NOT_ADMITTED",
            measurements=measurements,
            comparison=comparison,
            formal_median_s={m: float(np.median(v)) if v else None for m, v in formal.items()},
            CPU_over_GPU_speedup=float(np.median(formal["cpu"]) / np.median(formal["gpu"]))
            if admitted
            else None,
            solve_and_warmup_wall_s=sum_s,
            worker_wall_s=time.perf_counter() - start,
            peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
            training_ready=False,
        ),
    )
    phase("done", completed_solve_s=sum_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--worker", type=int)
    args = ap.parse_args()
    out = BASE / args.label
    if args.worker:
        worker(out / f"N{args.worker}", args.worker)
        return
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    protocol = js(ROOT / "configs/fast_pzt/f10k/protocol.json")["scaling"]
    write(out / "PROTOCOL.json", protocol)
    start = time.monotonic()
    spent = 0.0
    rows = []
    stop_reason = None
    previous = None
    for n in protocol["N"]:
        memory = estimate(n)
        if stop_reason:
            rows.append(dict(N=n, status="NOT_RUN", reason=stop_reason, static_memory=memory))
            continue
        if (
            memory["CPU_peak_array_estimate_bytes"] > protocol["CPU_memory_limit_GiB"] * 1024**3
            or memory["GPU_peak_array_estimate_bytes"] > protocol["GPU_memory_limit_GiB"] * 1024**3
        ):
            rows.append(dict(N=n, status="NOT_RUN", reason="STATIC_MEMORY_BUDGET", static_memory=memory))
            stop_reason = "STATIC_MEMORY_BUDGET"
            continue
        # Cost forecast only decides admission. Never plot it as measured timing.
        projected = None
        if previous and previous["N"] >= 128:
            projected = max(previous["formal_median_s"].values()) * (n / previous["N"]) ** 3
            if projected >= 10 or projected * 6 > 120 - spent:
                rows.append(
                    dict(
                        N=n,
                        status="NOT_RUN",
                        reason="CONSERVATIVE_COST_ADMISSION",
                        projected_single_s=projected,
                        static_memory=memory,
                    )
                )
                stop_reason = "CONSERVATIVE_COST_ADMISSION"
                continue
        if time.monotonic() - start >= 290 or spent >= 110:
            rows.append(dict(N=n, status="NOT_RUN", reason="STAGE_BUDGET_ADMISSION", static_memory=memory))
            stop_reason = "STAGE_BUDGET_ADMISSION"
            continue
        t = time.monotonic()
        timeout = None
        peakrss = 0
        with (out / f"N{n}.log").open("w") as log:
            proc = subprocess.Popen(
                [
                    sys.executable,
                    "-m",
                    "scripts.fast_pzt_f10k.scaling",
                    "--label",
                    args.label,
                    "--worker",
                    str(n),
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            latest = {}
            while proc.poll() is None:
                p = out / f"N{n}/PROGRESS.json"
                if p.exists():
                    latest = js(p)
                elapsed = time.monotonic() - latest.get("start_monotonic", time.monotonic())
                active = elapsed if latest.get("phase") == "solve" else 0.0
                partial = latest.get("completed_solve_s", 0.0) + active
                try:
                    fields = Path(f"/proc/{proc.pid}/status").read_text().splitlines()
                    peakrss = max(
                        peakrss, next(int(s.split()[1]) * 1024 for s in fields if s.startswith("VmRSS:"))
                    )
                except (FileNotFoundError, StopIteration):
                    pass
                if active > 10 or spent + partial > 120 or time.monotonic() - start > 299:
                    timeout = (
                        "SINGLE_SOLVE10S"
                        if active > 10
                        else "TOTAL_SOLVE120S"
                        if spent + partial > 120
                        else "STAGE300S"
                    )
                    proc.kill()
                    proc.wait(timeout=10)
                    break
                time.sleep(0.05)
            wall = time.monotonic() - t
        result = out / f"N{n}/RESULTS.json"
        if timeout:
            charged = latest.get("completed_solve_s", 0.0) + (
                time.monotonic() - latest["start_monotonic"] if latest.get("phase") == "solve" else 0.0
            )
            row = dict(
                N=n,
                status="TIME_BUDGET_EXCEEDED",
                reason=timeout,
                partial_state=latest,
                charged_solve_s=charged,
            )
            spent += charged
            stop_reason = timeout
        elif proc.returncode or not result.exists():
            row = dict(N=n, status="WORKER_FAILURE", exit_code=proc.returncode, latest=latest)
            spent += latest.get("completed_solve_s", 0.0)
            stop_reason = "WORKER_FAILURE"
        else:
            row = js(result)
            spent += row["solve_and_warmup_wall_s"]
            previous = row
            if row["status"] != "PASS":
                stop_reason = "PRIOR_NUMERICAL_OR_SINGLE_COST_FAILURE"
        row.update(process_wall_s=wall, controller_sampled_peak_RSS_bytes=peakrss, static_memory=memory)
        rows.append(row)
        write(out / "PARTIAL.json", dict(rows=rows, cumulative_solve_s=spent))
        print(n, row["status"], row.get("formal_median_s"), flush=True)
    write(
        out / "RESULTS.json",
        dict(
            rows=rows,
            solve_and_warmup_s=spent,
            stage_wall_s=time.monotonic() - start,
            solve_budget_s=120,
            stage_budget_s=300,
            all_cases_B=1,
            training_ready=False,
            timing_scope="Complete input preparation, solve, actual residual quality and identical JSON payload. Input load, reusable graph build, common Aref and disk artifacts separately charged to stage; no projected times as measurements.",
        ),
    )


if __name__ == "__main__":
    main()
