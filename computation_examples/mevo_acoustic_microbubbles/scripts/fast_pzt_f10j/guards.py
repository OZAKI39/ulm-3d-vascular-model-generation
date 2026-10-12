"""Negative linear solve and independent stable-merit kernel checks, no trajectories."""

import argparse
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f10j_gpu import CorrectedB2, energy_increment, sum_increment, check_linear
from mevo_mb.fast_sim.f10i_gpu import product
from .analyze_audit import dense, stable_delta
from .common import ROOT, BASE, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    wp.init()
    data = js(BASE / "audit_v1/AUDIT.json")
    g = CorrectedB2(32)
    g.c.assign(np.array(data["result"]["resistance_coefficients"]))
    g.normal.assign(np.array(data["result"]["resistance_normals"]))
    A = dense(data["result"]["resistance_coefficients"], data["result"]["resistance_normals"])
    rows = []
    old = None
    for e in data["evaluations"]:
        if not e["trial"]:
            old = e
            continue
        for key in ["pred", "q", "direction"]:
            getattr(g, key).assign(np.array(old[key]))
        g.launch(
            energy_increment, g.n, [g.pred, g.q, g.direction, g.c, g.normal, e["rate"], g.alpha, g.increment]
        )
        g.launch(sum_increment, 1, [g.increment, g.stats])
        actual = float(g.stats.numpy()[0])
        reference = stable_delta(A, old["pred"], old["q"], e["rate"] * np.array(old["direction"]), g.alpha)
        rhs = 1e-4 * e["rate"] * old["stats"][2]
        rows.append(
            dict(
                outer=e["outer"],
                rate=e["rate"],
                gpu_delta_E=actual,
                extended_precision_delta_E=reference,
                absolute_error=abs(actual - reference),
                armijo_rhs=rhs,
                decisions_equal=bool((actual <= rhs) == (reference <= rhs)),
            )
        )
    # Deliberately insufficient A solve. All inherited nonlinear and quality checks stay enabled.
    g.cgA._maxiter = 1
    with wp.ScopedCapture(device=g.device) as cap:
        g.linearA = g.cgA()
        g.launch(product, g.n, [g.u, g.b, g.true_residual, g.c, g.normal, g.a, g.current, 0, -1.0, 1.0])
        g.launch(
            check_linear,
            1,
            [g.b, g.cgA._r_and_z_buf[0], g.true_residual, g.invA, g.linear_count, g.linear_records],
        )
    g.graphA = cap.graph
    saved = js(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json")["rows"][-1]
    result = g.solve(saved["initial_positions_m"], saved["free_velocity_m_s"])
    expected = not result["converged"] and result["reason"] == "TRUE_OR_RECURSIVE_LINEAR_RESIDUAL"
    write(
        out / "RESULTS.json",
        dict(
            merit_rows=rows,
            independent_merit_decisions_pass=all(r["decisions_equal"] for r in rows),
            maximum_merit_absolute_error=max(r["absolute_error"] for r in rows),
            deliberate_incomplete_linear_solve=result,
            linear_failure_rejected=expected,
            status="PASS" if expected and all(r["decisions_equal"] for r in rows) else "FAIL",
        ),
    )


if __name__ == "__main__":
    main()
