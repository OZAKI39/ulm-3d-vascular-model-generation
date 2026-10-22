"""Two/four references first, then real32 and exact F10G spurious pair-step audit."""

import argparse
from time import perf_counter
import numpy as np
from mevo_mb.fast_sim.f10g_contact import resistance, soft_contact, compare
from mevo_mb.fast_sim.f10h_contact import geometric_contact, geometric_reference
from .common import ROOT, BASE, CASES, state, js, write, identity
from .analyze import gz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    protocol = js(ROOT / "configs/fast_pzt/f10h/protocol.json")
    r = 1.78e-6
    z = 6 * np.pi * 0.001 * r
    k = protocol["B2_kc_N_m"]
    start = perf_counter()
    assert k == 2 * js(ROOT / "configs/fast_pzt/f10g/stiffness.json")["kc_N_m"]
    rows = []

    def test(name, x, free, dt):
        A = resistance(x, r)["A"]
        result = geometric_contact(x, free, r, z, dt, k, A)
        ref = geometric_reference(x, free, r, z, dt, k, A)
        error = float(np.max(np.linalg.norm(result["positions_m"] - ref["positions_m"], axis=1)))
        passed = bool(
            result["converged"]
            and ref["reference_pass"]
            and error < 1e-12
            and np.min(result["pair_force_N"]) >= 0
            and np.linalg.norm(result["net_contact_force_N"]) < 1e-20
            and np.all(result["pair_force_N"][result["geometric_gap_m"] >= 0] == 0)
        )
        row = dict(
            name=name,
            bubble_count=len(x),
            dt_s=dt,
            reference_error_m=error,
            reference=ref,
            result=result,
            pass_=passed,
        )
        rows.append(row)
        return row

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
            row = test(name, x, free, dt)
            swapped = geometric_contact(x[::-1], free[::-1], r, z, dt, k, resistance(x[::-1], r)["A"])
            row["id_exchange_max_error_m"] = float(
                np.max(abs(swapped["positions_m"][::-1] - row["result"]["positions_m"]))
            )
            if tangent == 0:
                A = resistance(x, r)["A"]
                D = np.linalg.inv(A)
                B = np.array([[-1.0, 0, 0, 1, 0, 0]])
                exact = max(0, -gap - dt * float((B @ D @ free.ravel())[0])) / (
                    1 + dt * k / z * float((B @ D @ B.T)[0, 0])
                )
                row["analytic_overlap_error_m"] = abs(row["result"]["geometric_overlap_m"][0] - exact)
                row["pass_"] &= row["analytic_overlap_error_m"] < 1e-13
            row["pass_"] &= row["id_exchange_max_error_m"] < 1e-13
    write(out / "SMALL_SYSTEMS.json", dict(rows=rows, two_pass=all(v["pass_"] for v in rows)))
    if not all(v["pass_"] for v in rows):
        return
    configurations = [
        np.array([[0, 0, 0], [1.99, 0.08, 0.04], [-0.45, 1.94, 0.02], [-0.4, -0.55, 1.90]]),
        np.array([[-2.0, -0.06, 0], [0, 0, 0], [1.99, 0.10, 0.04], [0.1, 1.99, -0.1]]),
        np.array([[0, 0, 0], [1.99, 0, 0], [0.90, 1.75, 0.03], [0.85, 0.53, 1.72]]),
    ]
    for index, v in enumerate(configurations):
        x = v * r
        free = -0.001 * (v - v.mean(0)) + np.array([0.0003, -0.0002, 0.0001])
        for dt in [1e-5, 5e-6]:
            test("four_" + str(index), x, free, dt)
    smallpass = all(v["pass_"] for v in rows)
    write(out / "SMALL_SYSTEMS.json", dict(rows=rows, two_pass=True, four_pass=smallpass))
    if not smallpass:
        return
    statics = []
    for case in CASES:
        _, s, _ = state(case)
        x = np.array(s["initial_positions_m"])
        free = np.array(s["free_velocity_m_s"])
        v = test(case, x, free, 1e-5)
        v["budget_pass"] = (
            v["pass_"]
            and v["result"]["geometric_overlap_m"].max() <= 17.8e-9
            and v["result"]["max_single_iteration_m"] <= 178e-9
        )
        frozen = soft_contact(x, free, r, z, 1e-5, k, resistance(x, r)["A"])
        v["frozen_comparison"] = compare(v["result"]["positions_m"], frozen["positions_m"])
        statics.append(v)
    write(out / "DENSE_STATIC.json", dict(rows=statics, both_pass=all(v["budget_pass"] for v in statics)))
    # Exactly the archived false reactions, with true IDs and physical times.
    events = []
    replays = []
    for case in CASES:
        steps = gz(ROOT / f"results/fast_pzt/f10g/integration_v1/{case}_B2_K2/SOFT_STEPS.json.gz")
        for s in steps:
            old = s["result"]
            pairs = np.array(old["pairs"])
            fr = np.array(old["pair_force_N"])
            delta = np.array(old["geometric_overlap_m"])
            selected = np.flatnonzero((delta == 0) & (fr > 0))
            if not len(selected):
                continue
            x = np.array(s["initial_positions_m"])
            y = np.array(old["positions_m"])
            free = np.array(old["free_m_s"])
            new = geometric_contact(x, free, r, z, 1e-5, k, resistance(x, r)["A"])
            replays.append(dict(case=case, step=s["step"], result=new))
            for idx in selected:
                i, j = pairs[idx]
                v0 = x[j] - x[i]
                v1 = y[j] - y[i]
                d0 = np.linalg.norm(v0)
                d1 = np.linalg.norm(v1)
                events.append(
                    dict(
                        case=case,
                        step=s["step"],
                        time_s=s["step"] * 1e-5,
                        pair=[i, j],
                        actual_gap_m=d1 - 2 * r,
                        linear_gap_m=-old["frozen_normal_overlap_m"][idx],
                        initial_normal=v0 / d0,
                        updated_normal=v1 / d1,
                        normal_angle_rad=float(np.arccos(np.clip(np.dot(v0 / d0, v1 / d1), -1, 1))),
                        actual_distance_m=d1,
                        old_force_N=fr[idx],
                        new_gap_m=new["geometric_gap_m"][idx],
                        new_force_N=new["pair_force_N"][idx],
                        new_converged=new["converged"],
                    )
                )
    assert len(events) == 17
    write(
        out / "GEOMETRY_EVENTS.json",
        dict(
            events=events,
            replays=replays,
            geometry_mode="Same initial state, free force and fixed A/kc; independent single-step re-solves, not continued trajectories",
            old_event_count=17,
            new_false_reactions=sum(
                int(np.sum((v["result"]["geometric_gap_m"] >= 0) & (v["result"]["pair_force_N"] > 0)))
                for v in replays
            ),
            body_wall_s=perf_counter() - start,
        ),
    )
    print(
        "Two/four PASS; static",
        [
            (v["name"], v["budget_pass"], float(v["result"]["geometric_overlap_m"].max() * 1e9))
            for v in statics
        ],
        "old false events",
        len(events),
        flush=True,
    )


if __name__ == "__main__":
    main()
