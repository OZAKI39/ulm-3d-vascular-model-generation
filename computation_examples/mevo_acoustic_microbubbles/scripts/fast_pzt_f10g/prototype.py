"""Tasks B and conditional C: fixed-state CPU engineering/science evidence."""

import argparse
from time import perf_counter
import numpy as np
from mevo_mb.fast_sim.f10g_contact import (
    resistance,
    soft_contact,
    enumerated_reference,
    pbd,
    shape,
    compare,
    measure,
    contact_geometry,
)
from .common import BASE, ROOT, CASES, state, js, write, identity


def checks(r, zeta, k):
    two = []
    four = []
    setups = [
        ("separated", 200e-9, 0, 0),
        ("approaching", 5e-9, 0.002, 0),
        ("overlapping", -10e-9, 0, 0),
        ("separating", -10e-9, -0.003, 0),
        ("tangential", -10e-9, 0.001, 0.0008),
    ]
    for name, gap, speed, tangent in setups:
        x = np.array([[-(2 * r + gap) / 2, 0, 0], [(2 * r + gap) / 2, 0, 0]])
        free = np.array([[speed, tangent, 0], [-speed, -tangent, 0]])
        for mode in ["B1", "B2"]:
            A = None if mode == "B1" else resistance(x, r)["A"]
            for dt in [2.5e-6, 5e-6, 1e-5]:
                for mult in [0.5, 1, 2]:
                    s = soft_contact(x, free, r, zeta, dt, k * mult, A)
                    ref = enumerated_reference(x, free, r, zeta, dt, k * mult, A)
                    B, _, g = contact_geometry(x, r)
                    D = np.eye(6) if A is None else np.linalg.inv(A)
                    rel = float((B @ D @ B.T)[0, 0] / zeta)
                    b = float(-g[0] - dt * (B @ D @ free.ravel())[0])
                    exact = max(0, b) / (1 + dt * k * mult * rel)
                    assert s["converged"] and abs(s["frozen_normal_overlap_m"][0] - exact) < 1e-15
                    assert np.max(np.abs(s["positions_m"] - ref["positions_m"])) < 1e-12
                    assert s["force_position_consistency_m"] < 1e-13
                    assert s["pair_force_N"].min() >= 0 and np.linalg.norm(s["net_contact_force_N"]) < 1e-20
                    swapped = soft_contact(
                        x[::-1],
                        free[::-1],
                        r,
                        zeta,
                        dt,
                        k * mult,
                        None if A is None else resistance(x[::-1], r)["A"],
                    )
                    assert np.max(np.abs(s["positions_m"] - swapped["positions_m"][::-1])) < 1e-14
                    # Tangential velocity must not be changed by this fixed normal spring.
                    assert np.max(np.abs(s["corrected_velocity_m_s"][:, 1:] - free[:, 1:])) < 1e-12
                    explicit = (
                        max(0, -gap) - dt * k * mult * rel * max(0, -gap) - dt * (B @ D @ free.ravel())[0]
                    )
                    two.append(
                        dict(
                            name=name,
                            mode=mode,
                            multiplier=mult,
                            dt_s=dt,
                            gap_m=gap,
                            M_relative_m_N_s=rel,
                            delta_pred_m=b,
                            analytic_overlap_m=exact,
                            explicit_signed_overlap_m=explicit,
                            reference=ref,
                            result=s,
                        )
                    )
    x = np.array([[-r + 5e-9, 0, 0], [r - 5e-9, 0, 0]])
    free = np.array([[0.002, 0, 0], [-0.002, 0, 0]])
    hard = soft_contact(x, free, r, zeta, 1e-5, k * 1e4)
    assert hard["converged"] and hard["frozen_normal_overlap_m"].max() < 1e-12
    configurations = [
        np.array([[0, 0, 0], [1.99, 0.08, 0.04], [-0.45, 1.94, 0.02], [-0.4, -0.55, 1.90]]),
        np.array([[-2.0, -0.06, 0], [0, 0, 0], [1.99, 0.10, 0.04], [0.1, 1.99, -0.1]]),
        np.array([[0, 0, 0], [1.99, 0, 0], [0.90, 1.75, 0.03], [0.85, 0.53, 1.72]]),
    ]
    for index, x0 in enumerate(configurations):
        x = x0 * r
        free = -0.0015 * x0 + np.array(
            [
                [0, 0.00013, 0.00004],
                [0.00012, 0, 0.00003],
                [-0.00002, 0.00002, 0.0001],
                [0.00001, -0.00003, 0],
            ]
        )
        for mode in ["B1", "B2"]:
            A = None if mode == "B1" else resistance(x, r)["A"]
            for mult in [0.5, 1, 2]:
                s = soft_contact(x, free, r, zeta, 1e-5, k * mult, A)
                ref = enumerated_reference(x, free, r, zeta, 1e-5, k * mult, A)
                assert s["converged"] and np.max(np.abs(s["positions_m"] - ref["positions_m"])) < 1e-12
                assert s["force_position_consistency_m"] < 1e-13
                assert np.linalg.norm(s["net_contact_force_N"]) < 1e-20
                assert np.linalg.norm(s["centroid_contact_shift_m"]) < 1e-14
                p, rounds = pbd(s["prediction_m"], r)
                four.append(
                    dict(
                        configuration=index,
                        mode=mode,
                        multiplier=mult,
                        initial_positions_m=x,
                        free_m_s=free,
                        result=s,
                        reference=ref,
                        pbd_positions_m=p,
                        pbd_rounds=rounds,
                        comparison_to_pbd=compare(s["positions_m"], p),
                    )
                )
    return dict(
        two=two,
        four=four,
        hard_limit=hard,
        SOFT_CONTACT_TWO_BUBBLE_PASS=True,
        SOFT_CONTACT_FOUR_BUBBLE_PASS=True,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="prototype_local_v1")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    start = perf_counter()
    write(out / "IDENTITY.json", identity())
    proto = js(ROOT / "configs/fast_pzt/f10g/protocol.json")
    r = proto["radius_m"]
    zeta = 6 * np.pi * proto["viscosity_Pa_s"] * r
    stiffness = js(ROOT / "configs/fast_pzt/f10g/stiffness.json")
    k = stiffness["kc_N_m"]
    small = checks(r, zeta, k)
    write(out / "SMALL_SYSTEMS.json", small)
    print("Two/four independent references PASS", flush=True)
    rows = []
    for case in CASES:
        _, t, _ = state(case)
        x = np.array(t["initial_positions_m"])
        free = np.array(t["free_velocity_m_s"])
        A = resistance(x, r)["A"]
        original = np.array(t["proposal_positions_m"])
        for mode in ["B1", "B2"]:
            for mult in [0.5, 1, 2]:
                matrix = None if mode == "B1" else A
                s = soft_contact(x, free, r, zeta, 1e-5, k * mult, matrix)
                # Self convergence: 100x tighter equation tolerance, same physical model.
                tight = soft_contact(x, free, r, zeta, 1e-5, k * mult, matrix, atol_m=1e-15)
                good = (
                    s["converged"]
                    and s["geometric_overlap_m"].max() <= proto["residual_budget_m"]
                    and s["max_single_iteration_m"] <= proto["single_iteration_budget_m"]
                )
                row = dict(
                    case=case,
                    mode=mode,
                    multiplier=mult,
                    initial_positions_m=x,
                    initial_shape=shape(x, r),
                    result=s,
                    shape=shape(s["positions_m"], r),
                    relative_to_rejected_PBD=compare(s["positions_m"], original),
                    relative_to_initial=compare(s["positions_m"], x),
                    tight_reference=compare(s["positions_m"], tight["positions_m"]),
                    tight_converged=tight["converged"],
                    contact_budget_pass=bool(good),
                    full_domain_gate="NOT_EVALUATED_CPU_STATIC",
                    accepted_step=False,
                    edges_if_accepted=shape(s["positions_m"], r)["edges"],
                    original_accepted_edges=shape(x, r)["edges"],
                    timing=measure(lambda: soft_contact(x, free, r, zeta, 1e-5, k * mult, matrix)),
                )
                rows.append(row)
                print(
                    case,
                    mode,
                    mult,
                    "res nm",
                    s["geometric_overlap_m"].max() * 1e9,
                    "pass nm",
                    s["max_single_iteration_m"] * 1e9,
                    "budget",
                    good,
                    "iters",
                    len(s["iterations"]),
                    flush=True,
                )
    write(
        out / "DENSE_STATIC.json",
        dict(
            rows=rows,
            body_s=perf_counter() - start,
            training_ready=False,
            scope="Hypothetical one-step proposals, no time advanced; no accepted new trajectories",
        ),
    )


if __name__ == "__main__":
    main()
