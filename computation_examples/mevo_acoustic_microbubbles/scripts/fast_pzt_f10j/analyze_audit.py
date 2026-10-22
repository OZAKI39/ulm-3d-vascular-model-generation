"""Independent dense and extended-precision analysis of unchanged GPU traces."""

import numpy as np
from scipy.linalg import block_diag
from mevo_mb.fast_sim.f10g_contact import resistance, normals
from .common import ROOT, BASE, js, write, sha


def dense(c, n):
    c, n = np.asarray(c), np.asarray(n)
    A = np.eye(3 * len(c))
    for i in range(len(c)):
        for j in range(len(c)):
            if i != j:
                C = c[i, j] * np.outer(n[i, j], n[i, j])
                A[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] += C
                A[3 * i : 3 * i + 3, 3 * j : 3 * j + 3] -= C
    return A


def stable_delta(A, pred, q, step, alpha):
    """Same quadratic + unilateral squared-distance energy; no subtraction of totals."""
    A, pred, q, s = [np.asarray(v, dtype=np.longdouble) for v in (A, pred, q, step)]
    change = q @ A @ s + (s @ A @ s) / 2
    pred, q, s = [v.reshape(-1, 3) for v in (pred, q, s)]
    for i in range(len(q)):
        for j in range(i + 1, len(q)):
            dx = pred[j] - pred[i] + q[j] - q[i]
            ds = s[j] - s[i]
            d = np.sqrt(dx @ dx)
            dn = np.sqrt((dx + ds) @ (dx + ds))
            dd = (2 * (dx @ ds) + ds @ ds) / (dn + d)
            old = max(2 - d, 0)
            new = max(2 - d - dd, 0)
            change += np.longdouble(alpha) * (new - old) * (new + old) / 2
    return float(change)


def main():
    (BASE / "analysis").mkdir(parents=True, exist_ok=True)
    source = BASE / "audit_v1/AUDIT.json"
    data = js(source)
    row = js(ROOT / "results/fast_pzt/f10i/gpu_v2/RESULTS.json")["rows"][-1]
    x, free = np.array(row["initial_positions_m"]), np.array(row["free_velocity_m_s"])
    r, dt, kc = 1.78e-6, 1e-5, 0.04641319518378658
    zeta = 6 * np.pi * 0.001 * r
    alpha = dt * kc / zeta
    A = dense(data["result"]["resistance_coefficients"], data["result"]["resistance_normals"])
    Acpu = resistance(x, r)["A"]
    origins = []
    evals = data["evaluations"]
    for index, e in enumerate(evals):
        if e["trial"] and (not origins or e["outer"] != origins[-1]["outer"]):
            previous = evals[index - 1]
            assert not previous["trial"] and previous["outer"] == e["outer"]
            origins.append(previous)
    matrices = [A] + [A + dense(e["active"], e["normals"]) - np.eye(len(A)) for e in origins]
    rhss = [(free * dt / r).ravel()] + [np.array(e["rhs"]) for e in origins]
    linear = []
    for k, (matrix, b) in enumerate(zip(matrices, rhss)):
        kind = "A" if k == 0 else "H"
        start = 0 if k == 0 else (k - 1) * 129
        trace = np.asarray(data["linear"][kind]["values"])[start : start + 129]
        sol = np.asarray(data["linear"][kind]["solutions"])[start : start + 129]
        exact = np.linalg.solve(matrix, b)
        M = block_diag(*[np.linalg.inv(matrix[i : i + 3, i : i + 3]) for i in range(0, len(matrix), 3)])
        tol = max(1e-14, 1e-12 * np.linalg.norm(b))
        # Only pAp from iterations with previous residual above tolerance is relevant.
        used = trace[1:, 5][trace[:-1, 0] > tol]
        linear.append(
            dict(
                kind=kind,
                outer=k - 1 if k else None,
                tolerance=tol,
                symmetry_max=float(abs(matrix - matrix.T).max()),
                eigen_min=float(np.linalg.eigvalsh(matrix)[0]),
                eigen_max=float(np.linalg.eigvalsh(matrix)[-1]),
                preconditioner_eigen_min=float(np.linalg.eigvalsh(M)[0]),
                nonpositive_curvature_count=int(np.sum(used <= 0)),
                min_used_curvature=float(used.min()) if len(used) else None,
                recursive_final=float(trace[-1, 0]),
                true_final=float(trace[-1, 1]),
                recursive_relative=float(trace[-1, 0] / np.linalg.norm(b)),
                true_relative=float(trace[-1, 1] / np.linalg.norm(b)),
                recursive_pass=bool(trace[-1, 0] <= tol),
                true_pass=bool(trace[-1, 1] <= tol),
                dense_residual_final=float(np.linalg.norm(b - matrix @ sol[-1])),
                max_position_equivalent_error_m=float(
                    np.linalg.norm((sol[-1] - exact).reshape(-1, 3), axis=1).max() * r
                ),
                position_equivalent_error_curve_m=np.linalg.norm(
                    (sol - exact).reshape(len(sol), -1, 3), axis=2
                ).max(axis=1)
                * r,
                trace=trace,
            )
        )
    records = []
    byouter = {e["outer"]: e for e in origins}
    for e in evals:
        y = (np.array(e["pred"]) + np.array(e["q"])).reshape(-1, 3) * r + x.mean(0)
        i, j, d, n = normals(y)
        gap = d - 2 * r
        pair = kc * np.maximum(-gap, 0)
        force = np.zeros_like(y)
        np.add.at(force, i, -pair[:, None] * n)
        np.add.at(force, j, pair[:, None] * n)
        res = y - x - dt * np.linalg.solve(Acpu, free.ravel() + force.ravel() / zeta).reshape(-1, 3)
        record = dict(
            outer=e["outer"],
            trial=e["trial"],
            rate=e["rate"],
            energy=e["stats"][0],
            gradient_bound_m=r * e["stats"][1],
            equation_residual_m=float(np.linalg.norm(res, axis=1).max()),
            active_pairs=np.column_stack([i[gap < 0], j[gap < 0]]),
            pair_force_N=pair,
            worst_pair=[int(i[np.argmin(gap)]), int(j[np.argmin(gap)])],
            max_overlap_m=max(0.0, -gap.min()),
        )
        if e["trial"]:
            old = byouter[e["outer"]]
            rate = e["rate"]
            change = stable_delta(A, old["pred"], old["q"], rate * np.array(old["direction"]), alpha)
            arm = 1e-4 * rate * old["stats"][2]
            record.update(
                direction_derivative=old["stats"][2],
                original_energy_lhs=e["stats"][0],
                original_energy_rhs=old["stats"][0] + arm,
                original_energy_change=e["stats"][0] - old["stats"][0],
                stable_energy_change=change,
                armijo_change_rhs=arm,
                original_accept=e["stats"][0] <= old["stats"][0] + arm,
                stable_accept=change <= arm,
                max_candidate_displacement_m=float(
                    rate * r * np.linalg.norm(np.array(old["direction"]).reshape(-1, 3), axis=1).max()
                ),
                max_force_change_N=float(
                    kc
                    * r
                    * np.linalg.norm((np.array(e["force"]) - old["force"]).reshape(-1, 3), axis=1).max()
                ),
            )
        records.append(record)
    final_trials = [e for e in records if e["trial"] and e["outer"] == 4]
    report = dict(
        source_sha256=sha(source),
        instrumentation_position_delta_m=data["instrumentation_position_delta_m"],
        original_saved_position_delta_m=data["original_saved_position_delta_m"],
        gpu_cpu_operator_max_diff=float(abs(A - Acpu).max()),
        linear=linear,
        evaluations=records,
        diagnosis=dict(
            base_pcg_original_target_failed=not linear[0]["true_pass"],
            spd_conditions_hold=all(
                e["eigen_min"] > 0 and e["preconditioner_eigen_min"] > 0 and e["symmetry_max"] < 1e-12
                for e in linear
            ),
            final_direction_true_linear_pass=linear[-1]["true_pass"],
            final_trial_stable_accept_count=sum(e["stable_accept"] for e in final_trials),
            final_trial_original_accept_count=sum(e["original_accept"] for e in final_trials),
            old_stopping_rule_not_changed=True,
        ),
    )
    write(BASE / "analysis/ORIGINAL_AUDIT.json", report)
    print(report["diagnosis"])
    print(
        "base",
        {k: v for k, v in linear[0].items() if k not in ["trace", "position_equivalent_error_curve_m"]},
    )
    for e in final_trials[:3]:
        print({k: v for k, v in e.items() if k not in ["active_pairs", "pair_force_N"]})


if __name__ == "__main__":
    main()
