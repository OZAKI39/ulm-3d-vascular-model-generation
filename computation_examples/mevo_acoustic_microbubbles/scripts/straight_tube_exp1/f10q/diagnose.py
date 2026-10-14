"""Read-only actual frozen failure audit, before sequential-solver implementation."""

import time
import numpy as np
from mevo_mb.fast_sim.f10g_contact import contact_geometry, resistance
from scripts.straight_tube_exp1.f10o.analyze import load
from .common import ROOT, BASE, CONFIG, read, write, sha, identity


def main():
    start = time.perf_counter()
    p = read(CONFIG)
    R = p["radius_m"]
    dt = 1e-6
    zeta = 6 * np.pi * 0.001 * R
    folder = ROOT / "results/straight_tube_exp1/f10p/validation_v1/runs/STRONG_DT1"
    t, x = load(folder)
    x = x[-1]
    saved = read(folder / "LAST_STEP.json")
    r = saved["result"]
    y = np.array(r["positions_m"])
    F = np.array(saved["fields"]["total_force_N"])
    B, pairs, g0 = contact_geometry(x, R)
    By, p2, gy = contact_geometry(y, R)
    assert np.array_equal(pairs, p2)
    pairs = np.asarray(pairs)
    rel = x[pairs[:, 1]] - x[pairs[:, 0]]
    dr = (y - x)[pairs[:, 1]] - (y - x)[pairs[:, 0]]
    d = np.linalg.norm(rel, axis=1)
    n = rel / d[:, None]
    normal = np.einsum("ij,ij->i", n, dr)
    tangent = dr - normal[:, None] * n
    delta = np.linalg.norm(tangent, axis=1)
    axial = d + normal
    exact = delta**2 / (np.sqrt(axial**2 + delta**2) + axial)
    approx = delta**2 / (2 * axial)
    special = delta**2 / (4 * R)
    linear = g0 + normal
    lam = np.array(r["pair_force_N"])
    scale = max(np.linalg.norm(F, axis=1).max(), abs(lam).max(), 1e-18)
    weighted = abs(lam * gy) / scale
    worst = int(np.argmax(weighted))
    active = lam > 1e-16
    A = resistance(x, R)["A"]
    oldbalance = zeta * A @ ((y - x) / dt).ravel() - F.ravel() - B.T @ lam
    truebalance = zeta * A @ ((y - x) / dt).ravel() - F.ravel() - By.T @ lam
    records = []
    for k in range(len(pairs)):
        if active[k] or k == worst:
            records.append(
                dict(
                    pair=pairs[k],
                    old_normal=n[k],
                    new_normal=By[k].reshape(-1, 3)[pairs[k, 1]],
                    normal_rotation_deg=float(
                        np.degrees(np.arccos(np.clip(n[k] @ By[k].reshape(-1, 3)[pairs[k, 1]], -1, 1)))
                    ),
                    relative_normal_displacement_m=normal[k],
                    relative_tangential_displacement_m=delta[k],
                    initial_gap_m=g0[k],
                    linear_gap_m=linear[k],
                    actual_gap_m=gy[k],
                    exact_curvature_m=exact[k],
                    second_order_m=approx[k],
                    contact_special_delta2_over_4R_m=special[k],
                    reaction_N=lam[k],
                    weighted_true_comp_m=weighted[k],
                )
            )
    w = records[[tuple(v["pair"]) for v in records].index(tuple(pairs[worst]))]
    result = dict(
        identity=identity(),
        input_sha256={
            str((folder / f).relative_to(ROOT)): sha(folder / f)
            for f in ["trajectory.h5", "LAST_STEP.json", "INPUT.json"]
        },
        last_legal_time_s=t[-1],
        rejected_time_s=t[-1] + dt,
        force_scale_N=scale,
        worst_pair=w,
        active_records=records,
        all_positions_before_m=x,
        all_positions_proposal_m=y,
        all_pairs=pairs,
        all_reactions_N=lam,
        all_initial_gaps_m=g0,
        all_actual_gaps_m=gy,
        active_pairs_count=int(active.sum()),
        contact_network_edge_changes=int(np.sum((g0 <= 0.1 * R) != (gy <= 0.1 * R))),
        old_normal_balance_N=float(np.linalg.norm(oldbalance.reshape(-1, 3), axis=1).max()),
        true_normal_balance_N=float(np.linalg.norm(truebalance.reshape(-1, 3), axis=1).max()),
        saved_math_residual_m=r["mobility_equation_residual_m"],
        recomputed_true_comp_m=float(weighted.max()),
        saved_true_comp_m=r["actual_complementarity_equivalent_gap_m"],
        max_exact_identity_error_m=float(abs(gy - linear - exact).max()),
        max_saved_pair_gap_error_m=float(abs(gy - np.array(r["actual_gap_m"])).max()),
        max_saved_linear_gap_error_m=float(abs(linear - np.array(r["linear_gap_m"])).max()),
        conclusion="Finite tangential motion produces a positive curvature gap while the frozen local QP retains compressive reactions. Actual normals rotate; force balance mapped through actual normals is separately inconsistent. Stored pair mapping, units and weighted residual reproduce; no unit error found. Only this frozen proposal is diagnosed.",
        wall_s=time.perf_counter() - start,
    )
    assert np.isclose(weighted.max(), r["actual_complementarity_equivalent_gap_m"], rtol=1e-10)
    assert result["max_exact_identity_error_m"] < 1e-20
    write(BASE / "analysis/TRUE_GAP_FAILURE_DIAGNOSIS.json", result)
    print(
        {
            k: result[k]
            for k in [
                "worst_pair",
                "old_normal_balance_N",
                "true_normal_balance_N",
                "max_exact_identity_error_m",
                "wall_s",
            ]
        }
    )


if __name__ == "__main__":
    main()
