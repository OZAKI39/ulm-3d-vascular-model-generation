"""Bounded original 18 fixtures + exact original 32; no trajectory or timing sweep."""

import argparse
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f10j_gpu import CorrectedB2, OptimizedB2
from mevo_mb.fast_sim.f10h_contact import geometric_contact, geometric_reference
from mevo_mb.fast_sim.f10g_contact import resistance
from scripts.fast_pzt_f10i.gpu import fixtures, audit
from .common import ROOT, BASE, js, write, identity, sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--optimized", action="store_true")
    args = ap.parse_args()
    if args.optimized:
        assert js(BASE / "corrected_v1/RESULTS.json")["numerical_admitted"]
    cls = OptimizedB2 if args.optimized else CorrectedB2
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    write(
        out / "IDENTITY.json",
        dict(**identity(), candidate_protocol_sha256=sha(ROOT / "configs/fast_pzt/f10j/candidate_v1.json")),
    )
    wp.init()
    r, kc = 1.78e-6, 0.04641319518378658
    saved = js(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json")["rows"][-1]
    cases = list(fixtures(r)) + [
        ("saved32_O_1102", np.array(saved["initial_positions_m"]), np.array(saved["free_velocity_m_s"]), 1e-5)
    ]
    rows = []
    cache = {}
    for name, x, free, dt in cases:
        key = (len(x), dt)
        if key not in cache:
            cache[key] = cls(len(x), r, dt, kc)
        g = cache[key]
        A = resistance(x, r)["A"]
        cpu = geometric_contact(x, free, r, 6 * np.pi * 0.001 * r, dt, kc, A, atol_m=1e-14)
        ref = geometric_reference(x, free, r, 6 * np.pi * 0.001 * r, dt, kc, A)
        result = g.solve(x, free)
        identical = True
        if args.optimized:
            previous = js(BASE / "corrected_v1/RESULTS.json")["rows"][len(rows)]["gpu"]
            identical = bool(
                np.array_equal(result["positions_m"], previous["positions_m"])
                and result["reason"] == previous["reason"]
            )
        m = audit(x, free, dt, kc, A, result, cpu, ref, r)
        reject = max(cpu["geometric_overlap_m"]) > 17.8e-9
        acceptable = all(
            [
                identical,
                cpu["converged"],
                ref["reference_pass"],
                m["max_position_vs_cpu_m"] <= 1e-14,
                m["max_position_vs_independent_m"] <= 1e-14,
                m["max_force_vs_cpu_N"] <= 1e-15,
                m["max_force_vs_independent_N"] <= 1e-15,
                m["max_equation_residual_m"] <= 1e-14,
                m["net_force_N"] <= 1e-20,
                m["active_sets_equal"],
                m["edge_sets_equal"],
                m["force_actual_geometry_error_N"] < 1e-18,
                result["conservative_equation_bound_m"] <= 1e-14,
                np.all(
                    np.array(result["linear_audit"])[:, :2] <= np.array(result["linear_audit"])[:, 2, None]
                ),
                (result["reason"] == "OVERLAP_BUDGET" and not result["converged"])
                if reject
                else result["converged"],
            ]
        )
        rows.append(
            dict(
                name=name,
                n=len(x),
                dt_s=dt,
                initial_positions_m=x,
                free_velocity_m_s=free,
                cpu=cpu,
                independent_reference=ref,
                gpu=result,
                metrics=m,
                expected_quality_rejection=reject,
                test_pass=bool(acceptable),
                optimized_identical_to_corrected=identical if args.optimized else None,
            )
        )
        write(out / "RESULTS.json", dict(rows=rows, status="IN_PROGRESS", training_ready=False))
        print(
            name,
            len(x),
            dt,
            acceptable,
            result["reason"],
            result["conservative_equation_bound_m"],
            flush=True,
        )
    invalid = []
    g = cache[(2, 1e-5)]
    for name, x, free in [
        ("nonfinite", np.array([[np.nan, 0, 0], [1e-6, 0, 0]]), np.zeros((2, 3))),
        ("coincident", np.zeros((2, 3)), np.zeros((2, 3))),
    ]:
        try:
            g.solve(x, free)
            invalid.append(dict(name=name, rejected=False))
        except ValueError as exc:
            invalid.append(dict(name=name, rejected=True, message=str(exc)))
    # Deliberately severe, finite overlap is rejected by the unchanged correction budget.
    bad = g.solve(np.array([[0.0, 0, 0], [0.1 * r, 0, 0]]), np.zeros((2, 3)))
    invalid.append(dict(name="severe_overlap", rejected=not bad["converged"], result=bad))
    passed = all(z["test_pass"] for z in rows) and all(z["rejected"] for z in invalid)
    write(
        out / "RESULTS.json",
        dict(
            rows=rows,
            invalid=invalid,
            constructions=[
                dict(
                    n=k[0],
                    dt_s=k[1],
                    allocation_s=v.allocation_s,
                    old_graph_build_s=v.graph_build_s,
                    correction_graph_build_s=v.correction_graph_build_s,
                )
                for k, v in cache.items()
            ],
            status="GPU_FIXED_STATE_NUMERICAL_PASS" if passed else "GPU_B2_NOT_ADMITTED",
            numerical_admitted=passed,
            production=False,
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
