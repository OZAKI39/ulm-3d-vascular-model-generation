"""Conditional fixed-state Warp/CPU comparison; stop scientific expansion at failure."""

import argparse
from time import perf_counter
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f10i_gpu import FixedStateB2
from mevo_mb.fast_sim.f10h_contact import geometric_contact, geometric_reference
from mevo_mb.fast_sim.f10g_contact import resistance, normals
from .common import ROOT, BASE, js, write, identity


def fixtures(r):
    for name, gap, speed, tangent in [
        ("separated", 200e-9, 0, 0),
        ("approaching", 5e-9, 0.002, 0),
        ("overlapping", -10e-9, 0, 0),
        ("separating", -10e-9, -0.003, 0),
        ("tangential", -10e-9, 0.001, 0.0008),
        ("tangent_separates", -1e-10, 0, 0.03),
    ]:
        x = np.array([[-(2 * r + gap) / 2, 0, 0], [(2 * r + gap) / 2, 0, 0]])
        free = np.array([[speed, tangent, 0], [-speed, -tangent, 0]])
        for dt in [1e-5, 5e-6]:
            yield name, x, free, dt
    configs = [
        np.array([[0, 0, 0], [1.99, 0.08, 0.04], [-0.45, 1.94, 0.02], [-0.4, -0.55, 1.90]]),
        np.array([[-2.0, -0.06, 0], [0, 0, 0], [1.99, 0.10, 0.04], [0.1, 1.99, -0.1]]),
        np.array([[0, 0, 0], [1.99, 0, 0], [0.90, 1.75, 0.03], [0.85, 0.53, 1.72]]),
    ]
    for i, v in enumerate(configs):
        for dt in [1e-5, 5e-6]:
            yield "four_" + str(i), v * r, -0.001 * (v - v.mean(0)) + [0.0003, -0.0002, 0.0001], dt


def audit(x, free, dt, kc, A, gpu, cpu, reference, r):
    y = np.array(gpu["positions_m"])
    i, j, d, n = normals(y)
    gap = d - 2 * r
    pair = kc * np.maximum(-gap, 0)
    f = np.zeros_like(y)
    np.add.at(f, i, -pair[:, None] * n)
    np.add.at(f, j, pair[:, None] * n)
    zeta = 6 * np.pi * 0.001 * r
    residual = y - x - dt * np.linalg.solve(A, free.ravel() + f.ravel() / zeta).reshape(-1, 3)
    gc = np.array(cpu["geometric_gap_m"])
    D = np.array(gpu["resistance_coefficients"])
    N = np.array(gpu["resistance_normals"])
    reconstructed = np.eye(x.size)
    for a in range(len(x)):
        for b in range(a + 1, len(x)):
            C = D[a, b] * np.outer(N[a, b], N[a, b])
            s = slice(3 * a, 3 * a + 3)
            t = slice(3 * b, 3 * b + 3)
            reconstructed[s, s] += C
            reconstructed[t, t] += C
            reconstructed[s, t] -= C
            reconstructed[t, s] -= C
    return dict(
        max_position_vs_cpu_m=float(np.linalg.norm(y - cpu["positions_m"], axis=1).max()),
        max_position_vs_independent_m=float(np.linalg.norm(y - reference["positions_m"], axis=1).max()),
        max_force_vs_cpu_N=float(
            np.linalg.norm(np.array(gpu["contact_force_N"]) - cpu["contact_force_N"], axis=1).max()
        ),
        max_force_vs_independent_N=float(
            np.linalg.norm(np.array(gpu["contact_force_N"]) - reference["contact_force_N"], axis=1).max()
        ),
        force_actual_geometry_error_N=float(np.max(abs(f - gpu["contact_force_N"]))),
        active_sets_equal=bool(np.array_equal(gap < 0, gc < 0)),
        edge_sets_equal=bool(np.array_equal(gap <= 0.1 * r, gc <= 0.1 * r)),
        max_overlap_m=float(max(0, -gap.min())),
        max_equation_residual_m=float(np.linalg.norm(residual, axis=1).max()),
        net_force_N=float(np.linalg.norm(np.sum(gpu["contact_force_N"], axis=0))),
        torque_N_m=float(np.linalg.norm(np.sum(np.cross(y - y.mean(0), gpu["contact_force_N"]), axis=0))),
        resistance_operator_max_error=float(np.max(abs(A - reconstructed))),
        minimum_A_eigenvalue=float(np.linalg.eigvalsh(A)[0]),
        max_single_iteration_m=gpu["max_single_iteration_m"],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    assert js(BASE / "analysis/GPU_START_GATE.json")["admitted"]
    out = BASE / a.label
    out.mkdir(exist_ok=False, parents=True)
    write(out / "IDENTITY.json", identity())
    protocol = js(ROOT / "configs/fast_pzt/f10i/protocol.json")
    r = protocol["radius_m"]
    kc = protocol["kc_N_m"]
    q = protocol["gpu_comparison"]
    wp.init()
    cache = {}
    rows = []
    constructions = []

    def test(name, x, free, dt):
        key = (len(x), dt)
        if key not in cache:
            cache[key] = FixedStateB2(len(x), r, dt, kc)
            g = cache[key]
            constructions.append(
                dict(
                    n=len(x),
                    dt_s=dt,
                    allocation_s=g.allocation_s,
                    graph_build_s=g.graph_build_s,
                    explicit_buffer_bytes=g.buffer_bytes,
                )
            )
        g = cache[key]
        op = resistance(x, r)
        A = op["A"]
        zeta = 6 * np.pi * 0.001 * r
        cpu = geometric_contact(x, free, r, zeta, dt, kc, A, atol_m=1e-14)
        ref = geometric_reference(x, free, r, zeta, dt, kc, A)
        result = g.solve(x, free)
        m = audit(x, free, dt, kc, A, result, cpu, ref, r)
        passed = bool(
            result["converged"]
            and cpu["converged"]
            and ref["reference_pass"]
            and m["max_position_vs_cpu_m"] <= q["position_m"]
            and m["max_position_vs_independent_m"] <= q["position_m"]
            and m["max_force_vs_cpu_N"] <= q["contact_force_N"]
            and m["max_force_vs_independent_N"] <= q["contact_force_N"]
            and m["max_equation_residual_m"] <= q["equation_residual_m"]
            and m["net_force_N"] <= q["net_force_N"]
            and m["active_sets_equal"]
            and m["edge_sets_equal"]
            and (m["max_overlap_m"] <= q["overlap_m"]) == (max(cpu["geometric_overlap_m"]) <= q["overlap_m"])
            and m["max_single_iteration_m"] <= q["single_iteration_m"]
            and m["force_actual_geometry_error_N"] < 1e-18
        )
        row = dict(
            name=name,
            n=len(x),
            dt_s=dt,
            initial_positions_m=x,
            free_velocity_m_s=free,
            cpu=cpu,
            independent_reference=ref,
            gpu=result,
            metrics=m,
            contact_quality="PASS" if m["max_overlap_m"] <= q["overlap_m"] else "EXPECTED_BUDGET_REJECTION",
            pass_=passed,
        )
        rows.append(row)
        write(
            out / "RESULTS.json",
            dict(
                rows=rows,
                constructions=constructions,
                status="IN_PROGRESS" if passed else "GPU_B2_NOT_ADMITTED",
                training_ready=False,
            ),
        )
        print(
            name,
            len(x),
            dt,
            passed,
            m["max_position_vs_cpu_m"],
            m["max_equation_residual_m"],
            result["reason"],
            flush=True,
        )
        return passed

    try:
        for name, x, free, dt in fixtures(r):
            if not test(name, x, free, dt):
                return
        # Exact stored pre-step state, not a rejected proposal or a restart.
        from scripts.fast_pzt_f10h.common import state

        _, s, _ = state("O_1000_1102")
        x = np.array(s["initial_positions_m"])
        free = np.array(s["free_velocity_m_s"])
        if not test("saved32_O_1102", x, free, 1e-5):
            return
        g = cache[(32, 1e-5)]
        samples = []
        # One explicit warmup, then exactly five matched static timing samples.
        g.solve(x, free)
        for _ in range(5):
            t = perf_counter()
            A = resistance(x, r)["A"]
            assembly = perf_counter() - t
            cpu = geometric_contact(x, free, r, 6 * np.pi * 0.001 * r, 1e-5, kc, A, atol_m=1e-14)
            elapsed = perf_counter() - t
            candidate = g.solve(x, free)
            assert candidate["converged"]
            samples.append(
                dict(
                    cpu_assembly_s=assembly,
                    cpu_contact_s=cpu["total_s"],
                    cpu_whole_s=elapsed,
                    gpu=candidate["timing"],
                )
            )
        cpu_med = float(np.median([s["cpu_whole_s"] for s in samples]))
        gpu_med = float(np.median([s["gpu"]["whole_step_wall_s"] for s in samples]))
        status = "GPU_FIXED_STATE_NUMERICAL_PASS" if gpu_med < cpu_med else "GPU_B2_NOT_ADMITTED"
        write(
            out / "RESULTS.json",
            dict(
                rows=rows,
                constructions=constructions,
                numerical_status="PASS",
                status=status,
                performance_samples=samples,
                cpu_median_whole_s=cpu_med,
                gpu_median_whole_s=gpu_med,
                performance_reason="Fixed-state feasibility, no production claim"
                if gpu_med < cpu_med
                else "Full GPU step is slower than CPU; stop expansion under preregistered rule",
                production=False,
                training_ready=False,
            ),
        )
        print("Numerical PASS; cost CPU/GPU", cpu_med, gpu_med, status, flush=True)
    except Exception as e:
        write(
            out / "ENGINEERING_FAILURE.json",
            dict(exception_type=type(e).__name__, message=str(e), completed_scientific_cases=len(rows)),
        )
        raise


if __name__ == "__main__":
    main()
