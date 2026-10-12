"""Task A: same frozen inputs, independent resistance reference and old PBD."""

import argparse
from time import perf_counter
import numpy as np
from scipy import linalg, sparse
from scipy.sparse.linalg import spsolve
from mevo_mb.fast_sim.f10g_contact import resistance, jacobi, balance, pbd, overlaps, measure
from .common import ROOT, BASE, CASES, state, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="lubrication_local_v1")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    start = perf_counter()
    results = []
    compression = []
    write(out / "IDENTITY.json", identity())
    for case in CASES:
        data, t, path = state(case)
        x = np.array(t["initial_positions_m"])
        free = np.array(t["free_velocity_m_s"])
        flow = np.array(t["flow_velocity_m_s"])
        force = np.array(t["total_force_N"])
        r = data["radius_m"]
        dt = t["dt_s"]
        zeta = 6 * np.pi * 0.001 * r
        assert np.max(np.abs(free - flow - force / zeta)) < 1e-12
        op = resistance(x, r)
        A = op["A"]
        eigs = linalg.eigvalsh(A)
        ref = linalg.solve(A, free.ravel(), assume_a="pos").reshape(-1, 3)
        sparse_ref = spsolve(sparse.csr_matrix(A), free.ravel()).reshape(-1, 3)
        assert (
            np.max(np.abs(ref - sparse_ref)) < 1e-12
            and np.max(np.abs(A - A.T)) < 1e-13
            and eigs[0] > 0.9999999
        )
        i, j = op["pairs"].T
        n = op["normals"]
        scale = np.linalg.norm(free, axis=1).max()
        compression.extend(np.maximum(0, 0.5 * np.sum((force[i] - force[j]) * n, axis=1)))
        modes = {str(k): jacobi(op, free, k) for k in [2, 4, 8]}
        modes["reference"] = ref
        assert np.max(np.abs(modes["2"] - t["local_velocity_m_s"])) < 1e-12
        velocities = {}
        rows = []
        for mode, u in modes.items():
            balanced = balance(op, free, u) if mode != "reference" else u
            if mode == "2":
                assert np.max(np.abs(x + dt * balanced - t["predicted_positions_m"])) < 1e-15
            prediction = x + dt * balanced
            final, rounds = pbd(prediction, r)
            if mode == "2":
                for old, new in zip(t["rounds"], rounds):
                    assert np.max(np.abs(np.array(old["positions_m"]) - new["positions_m"])) < 1e-15
            residuals = {}
            for title, v in [("raw", u), ("balanced", balanced)]:
                res = (A @ v.ravel() - free.ravel()).reshape(-1, 3)
                pre = np.linalg.solve(op["blocks"], res[..., None])[..., 0]
                err = v - ref
                normalerr = np.sum((err[j] - err[i]) * n, axis=1)
                residuals[title] = dict(
                    residual_per_ID_m_s=res,
                    max_residual_m_s=np.linalg.norm(res, axis=1).max(),
                    max_residual_N=zeta * np.linalg.norm(res, axis=1).max(),
                    normalized_residual=np.linalg.norm(res, axis=1).max() / scale,
                    preconditioned_per_ID_m_s=pre,
                    normalized_preconditioned=np.linalg.norm(pre, axis=1).max() / scale,
                    max_velocity_error_m_s=np.linalg.norm(err, axis=1).max(),
                    normalized_velocity_error=np.linalg.norm(err, axis=1).max() / scale,
                    velocity_error_per_ID_m_s=err,
                    mean_velocity_m_s=v.mean(axis=0),
                    worst_normal_error_pair=op["pairs"][np.argmax(np.abs(normalerr))],
                    normal_error_per_pair_m_s=normalerr,
                )
            fn = (
                (lambda: linalg.solve(A, free.ravel(), assume_a="pos"))
                if mode == "reference"
                else (lambda k=int(mode): jacobi(op, free, k))
            )
            runtime = measure(fn)
            velocities[mode] = balanced
            row = dict(
                mode=mode,
                residuals=residuals,
                velocity_m_s=balanced,
                relative_normal_velocity_m_s=np.sum((balanced[j] - balanced[i]) * n, axis=1),
                prediction_positions_m=prediction,
                predicted_overlap_m=overlaps(prediction, r).max(),
                predicted_contact_pairs=int(np.sum(overlaps(prediction, r) > 0)),
                final_positions_m=final,
                pbd_rounds=rounds,
                final_overlap_m=overlaps(final, r).max(),
                max_single_pass_m=max(v["max_correction_m"] for v in rounds),
                contact_budget_pass=bool(
                    overlaps(final, r).max() <= 0.01 * r
                    and max(v["max_correction_m"] for v in rounds) <= 0.1 * r
                ),
                timing=runtime,
                pbd_timing=measure(lambda: pbd(prediction, r)),
            )
            rows.append(row)
        np.savez_compressed(
            out / (case + "_operator.npz"),
            A=A,
            free_m_s=free,
            x_m=x,
            flow_m_s=flow,
            force_N=force,
            reference_m_s=ref,
            **{"velocity_" + k: v for k, v in velocities.items()},
        )
        results.append(
            dict(
                case=case,
                start_s=t["start_s"],
                end_s=t["attempted_end_s"],
                trace=str(path.relative_to(ROOT)),
                resistance_units="dimensionless A; zeta*A in N s/m",
                zeta_N_s_m=zeta,
                symmetric=True,
                positive_definite=True,
                min_eigenvalue=eigs[0],
                max_eigenvalue=eigs[-1],
                condition_number=eigs[-1] / eigs[0],
                sparse_vs_cholesky_max_m_s=np.max(np.abs(ref - sparse_ref)),
                operator_initialization_timing=measure(lambda: resistance(x, r)),
                pairs=op["pairs"],
                normals=n,
                gaps_m=op["gaps_m"],
                coefficients=op["coefficients"],
                floor_flags=op["floored"],
                floor_pairs=int(op["floored"].sum()),
                groups=op["groups"],
                max_free_speed_m_s=scale,
                free_velocity_m_s=free,
                flow_velocity_m_s=flow,
                acoustic_free_velocity_m_s=force / zeta,
                pressure_abs_Pa=np.linalg.norm(t["pressure_complex_Pa"], axis=1),
                rows=rows,
            )
        )
    k = max(compression) / 8.9e-9
    write(
        out / "STIFFNESS_PREREGISTRATION.json",
        dict(
            rule="protocol.json: pooled maximum acoustic normal compression / 8.9nm",
            normal_compression_N=compression,
            max_compression_N=max(compression),
            target_overlap_m=8.9e-9,
            kc_N_m=k,
            multipliers=[0.5, 1, 2],
            registered_before_soft_results=True,
            scope="Frozen-state numerical penalty parameter, not calibrated shell contact",
        ),
    )
    write(
        out / "RESULTS.json",
        dict(
            cases=results,
            LUBRICATION_REFERENCE_PASS=True,
            LUBRICATION_CONVERGENCE_CHECKED=True,
            body_s=perf_counter() - start,
            training_ready=False,
        ),
    )
    print("Task A PASS; kc preregistered", k)
    for c in results:
        print(
            c["case"],
            [
                (
                    r["mode"],
                    r["residuals"]["balanced"]["normalized_preconditioned"],
                    r["predicted_overlap_m"] * 1e9,
                    r["final_overlap_m"] * 1e9,
                )
                for r in c["rows"]
            ],
        )


if __name__ == "__main__":
    main()
