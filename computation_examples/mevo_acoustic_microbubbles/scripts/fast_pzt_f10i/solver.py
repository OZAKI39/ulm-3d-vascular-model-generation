"""Same saved physical states, bounded tighter CPU solves and independent roots."""

from time import perf_counter
import sys
import numpy as np
import scipy
from mevo_mb.fast_sim.f10g_contact import resistance, normals
from mevo_mb.fast_sim.f10h_contact import geometric_contact, geometric_reference
from .common import ROOT, BASE, js, write, gz, identity


def metrics(a, b, r, kc):
    x = np.array(a["positions_m"])
    y = np.array(b["positions_m"])
    ia, ja, da, na = normals(x)
    _, _, db, nb = normals(y)
    ga = da - 2 * r
    gb = db - 2 * r
    fa = np.array(a["contact_force_N"])
    fb = np.array(b["contact_force_N"])
    pair = np.flatnonzero((ia == 7) & (ja == 24))[0]
    return dict(
        max_position_difference_m=float(np.linalg.norm(x - y, axis=1).max()),
        max_force_difference_N=float(np.linalg.norm(fa - fb, axis=1).max()),
        max_pair_force_difference_N=float(np.max(abs(kc * np.maximum(-ga, 0) - kc * np.maximum(-gb, 0)))),
        active_contact_set_equal=bool(np.array_equal(ga < 0, gb < 0)),
        edge_set_equal=bool(np.array_equal(ga <= 0.1 * r, gb <= 0.1 * r)),
        max_normal_difference=float(np.linalg.norm(na - nb, axis=1).max()),
        pair_7_24_gap_difference_m=float(gb[pair] - ga[pair]),
        pair_7_24_signed_edge_before_m=float(ga[pair] - 0.1 * r),
        pair_7_24_signed_edge_after_m=float(gb[pair] - 0.1 * r),
    )


def main():
    start = perf_counter()
    out = BASE / "analysis"
    protocol = js(ROOT / "configs/fast_pzt/f10i/protocol.json")
    r = protocol["radius_m"]
    zeta = 6 * np.pi * protocol["viscosity_Pa_s"] * r
    kc = protocol["kc_N_m"]
    states = gz(out / "CPU_SAVED_STATES.json.gz")
    rows = []
    for name, s in states.items():
        old = s["archived_result"]
        x = np.array(s["initial_positions_m"])
        free = np.array(old["free_m_s"])
        dt = old["dt_s"]
        t = perf_counter()
        op = resistance(x, r)
        assembly = perf_counter() - t
        A = op["A"]
        ref = geometric_reference(x, free, r, zeta, dt, kc, A)
        assert ref["reference_pass"]
        variants = []
        for tol in protocol["cpu_tolerances_m"]:
            result = geometric_contact(
                x, free, r, zeta, dt, kc, A, max_iterations=protocol["cpu_max_outer_iterations"], atol_m=tol
            )
            comparison = metrics(old, result, r, kc)
            reference_error = metrics(result, ref, r, kc)
            gap = np.array(result["geometric_gap_m"])
            force = np.array(result["pair_force_N"])
            assert np.all(force >= 0) and np.all(force[gap >= 0] == 0)
            variants.append(
                dict(tolerance_m=tol, result=result, vs_archived=comparison, vs_independent=reference_error)
            )
            if not result["converged"]:
                break
        rows.append(
            dict(
                name=name,
                mode=s["mode"],
                accepted_endpoint_s=s["accepted_endpoint_s"],
                archived_residual_m=old["mobility_equation_residual_m"],
                archived_outer_iterations=len(old["iterations"]),
                assembly_s=assembly,
                minimum_resistance_eigenvalue=float(np.linalg.eigvalsh(A)[0]),
                maximum_resistance_eigenvalue=float(np.linalg.eigvalsh(A)[-1]),
                reference=ref,
                variants=variants,
            )
        )
        print(
            name,
            [
                (v["tolerance_m"], v["result"]["converged"], v["vs_archived"]["max_position_difference_m"])
                for v in variants
            ],
            flush=True,
        )
    tight = [v for row in rows for v in row["variants"] if v["tolerance_m"] == 1e-15]
    passed = len(tight) == len(rows) and all(v["result"]["converged"] for v in tight)
    summary = dict(
        all_tight_converged=passed,
        saved_state_count=len(rows),
        max_tight_position_change_m=max(v["vs_archived"]["max_position_difference_m"] for v in tight),
        max_tight_contact_force_change_N=max(v["vs_archived"]["max_force_difference_N"] for v in tight),
        max_tight_pair_7_24_gap_change_m=max(
            abs(v["vs_archived"]["pair_7_24_gap_difference_m"]) for v in tight
        ),
        all_tight_active_sets_equal=all(v["vs_archived"]["active_contact_set_equal"] for v in tight),
        all_tight_edge_sets_equal=all(v["vs_archived"]["edge_set_equal"] for v in tight),
        max_tight_independent_position_error_m=max(
            v["vs_independent"]["max_position_difference_m"] for v in tight
        ),
        max_independent_reference_residual_m=max(row["reference"]["residual_m"] for row in rows),
    )
    write(
        out / "SOLVER_TOLERANCE.json",
        dict(
            source=identity(),
            summary=summary,
            rows=rows,
            environment=dict(
                python=sys.version, numpy=np.__version__, scipy=scipy.__version__, executable=sys.executable
            ),
            body_s=perf_counter() - start,
            new_trajectory=False,
            interpretation="Same input one-step re-solves; independent library root. Tiny local change is not a global trajectory error bound or proof of temporal order. No stiffness/model change.",
        ),
    )
    print(summary)


if __name__ == "__main__":
    main()
