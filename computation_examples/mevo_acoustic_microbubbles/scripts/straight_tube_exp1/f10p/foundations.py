"""Few deterministic analytic/independent checks; never a random stress test."""

import time
import numpy as np
from mevo_mb.fast_sim.hard_sphere import hard_step, enumerated_hard_reference
from mevo_mb.fast_sim.f10g_contact import resistance
from mevo_mb.fast_sim.f10h_contact import geometric_contact
from .common import CONFIG, read


def run():
    start = time.perf_counter()
    p = read(CONFIG)
    p.update(read(CONFIG.parent / "algorithm_revision_2.json"))
    R = p["radius_m"]
    zeta = 6 * np.pi * 0.001 * R
    dt = 1e-6
    pair = np.array([[-R, 0, 0], [R, 0, 0]])
    inward = np.array([[1e-9, 0, 0], [-1e-9, 0, 0]])
    square = np.array([[-R, -R, 0], [R, -R, 0], [-R, R, 0], [R, R, 0]])
    cases = [
        ("separating", pair, -inward),
        ("compressing", pair, inward),
        ("release", pair, -0.5 * inward),
        ("no_contact", pair * 4, inward),
        ("four_simultaneous", square, -square / R * 1e-9),
        ("four_network_translation", square, np.tile([1e-10, 2e-10, 0], (4, 1))),
    ]
    records = []
    for name, x, F in cases:
        A = resistance(x, R)["A"]
        r = hard_step(x, F, R, zeta, dt, A, p)
        ref = enumerated_hard_reference(x, F, R, zeta, dt, A)
        b2 = geometric_contact(x, F / zeta, R, zeta, dt, 0.04641319518378658, A, atol_m=1e-13)
        error = float(np.linalg.norm(r["positions_m"] - ref["positions_m"], axis=1).max())
        checks = dict(
            gates=r["converged"],
            independent_position=error < 1e-13,
            net_internal_force=np.linalg.norm(r["net_contact_force_N"]) < 1e-22,
            COM_force_balance=np.linalg.norm((r["positions_m"] - x).mean(0) - dt * F.mean(0) / zeta) < 1e-13,
        )
        if name in ["separating", "release", "no_contact", "four_network_translation"]:
            checks["no_spurious_reaction"] = bool(abs(r["pair_force_N"]).max(initial=0) < 1e-16)
        if name == "compressing":
            checks["analytic_static"] = bool(np.max(abs(r["positions_m"] - x)) < 1e-13)
            checks["analytic_reaction"] = bool(abs(r["pair_force_N"][0] - 1e-9) < 1e-16)
        records.append(
            dict(
                name=name,
                checks=checks,
                status="PASS" if all(checks.values()) else "FAIL",
                input=dict(positions_m=x, force_N=F, dt_s=dt, A=A),
                hard=r,
                reference=ref,
                B2=b2,
                independent_max_error_m=error,
            )
        )
    # Rejection fixtures are synthetic, distinct from physical trajectories.
    rejects = {}
    for name, A in [("indefinite", -np.eye(6)), ("nonsymmetric", np.eye(6) + np.triu(np.ones((6, 6)), 1))]:
        try:
            hard_step(pair, inward, R, zeta, dt, A, p)
            rejects[name] = False
        except (ValueError, np.linalg.LinAlgError):
            rejects[name] = True
    limited = dict(p, hard_solver="osqp", osqp=dict(p["osqp"], max_iter=1, polishing=False))
    limited_result = hard_step(pair, inward, R, zeta, dt, np.eye(6), limited)
    rejects["iteration_limit_not_accepted"] = not limited_result["converged"]
    return dict(
        status="PASS" if all(r["status"] == "PASS" for r in records) and all(rejects.values()) else "FAIL",
        cases=records,
        rejection_checks=rejects,
        limited_solver=limited_result,
        wall_s=time.perf_counter() - start,
    )
