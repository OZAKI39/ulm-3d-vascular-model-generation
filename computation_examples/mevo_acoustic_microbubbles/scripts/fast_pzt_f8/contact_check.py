"""Limited 2/4/actual32 contact study; no navigation promotion without admission."""

import argparse
import time
import numpy as np
import warp as wp
from scipy.optimize import nnls
from mevo_mb.fast_sim.f8_contact import prepare, solve
from mevo_mb.fast_sim import kernels as k
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import ROOT, BASE, identity, js, write


def pbd(free, p):
    begin = time.perf_counter()
    x = wp.array(free, dtype=wp.vec3d, device="cuda:0")
    a, b = wp.zeros_like(x), wp.zeros_like(x)
    diag = wp.zeros((len(free), 11), dtype=wp.float64, device="cuda:0")
    error = wp.zeros(1, dtype=int, device="cuda:0")
    wp.load_module(module=k, device="cuda:0")
    # Compile before measuring the warm kernel span, but retain full setup time.
    with wp.ScopedCapture(device="cuda:0") as cap:
        source = x
        for it in range(p.contact_iterations):
            out = a if it % 2 == 0 else b
            wp.launch(
                k.contact,
                dim=len(free),
                inputs=[source, out, len(free), p.radius_m, p.contact_relaxation, diag, error],
                device="cuda:0",
            )
            source = out
    setup = time.perf_counter() - begin
    start, end = wp.Event(enable_timing=True), wp.Event(enable_timing=True)
    begin = time.perf_counter()
    wp.record_event(start)
    wp.capture_launch(cap.graph)
    wp.record_event(end)
    wp.synchronize_event(end)
    return dict(
        positions_m=source.numpy(),
        diagnostics=diag.numpy(),
        error=int(error.numpy()[0]),
        setup_wall_s=setup,
        host_wall_s=time.perf_counter() - begin,
        gpu_launch_span_s=wp.get_event_elapsed_time(start, end) / 1000,
    )


def overlap(x, radius):
    d = np.linalg.norm(x[:, None] - x[None, :], axis=-1)
    np.fill_diagonal(d, np.inf)
    return max(0.0, float(2 * radius - d.min()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    args = ap.parse_args()
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    proto = js(ROOT / "configs/fast_pzt/f8/protocol.json")["contact_prototype"]
    write(out / "IDENTITY.json", identity())
    write(out / "REGISTERED_PROTOCOL.json", proto)
    p = F2Physics()
    radius = p.radius_m
    two = np.array([[-1.0, 0, 0], [1.0, 0, 0]]) * radius
    four = np.array([[0, 0, 0], [2.0, 0, 0], [0, 2.0, 0], [0.4, 0.3, 2.1]]) * radius
    fixed = np.load(ROOT / "results/fast_pzt/f7/fixed_v1/fixed_states.npz")
    actual = fixed["mature_positions_m"]
    cases = {
        "two": (two, two * 0.995),
        "four": (four, four + 0.005 * (four.mean(0) - four)),
        "actual32": (actual, actual + 1e-5 * fixed["mature_M2_velocity_m_s"]),
    }
    rows = {}
    arrays = {}
    begun = time.perf_counter()
    for name, (reference, free) in cases.items():
        prepared = prepare(reference, free, radius, p.viscosity_Pa_s)
        result = solve(prepared, proto["passes"])
        fine = solve(prepared, proto["diagnostic_refinement_passes"])
        baseline = pbd(free, p)
        dx = result["displacement_m"]
        lam = result["lambda_m"]
        impulse = (prepared["J"].T @ lam).reshape(-1, 3) / prepared["mu0"]
        np.testing.assert_allclose(dx.ravel(), prepared["B"] @ lam, rtol=1e-10, atol=1e-19)
        np.testing.assert_allclose(impulse.sum(0), 0, rtol=0, atol=1e-28)
        torque = np.cross(free - free.mean(0), impulse).sum(0)
        torque_scale = max(
            float(np.sum(np.linalg.norm(free - free.mean(0), axis=1) * np.linalg.norm(impulse, axis=1))),
            1e-100,
        )
        assert np.linalg.norm(torque) < 1e-12 * torque_scale
        assert np.all(lam >= 0) and np.isfinite(dx).all()
        reference_error = None
        if name != "actual32":
            A = prepared["J"] @ prepared["B"]
            L = np.linalg.cholesky(A)
            ref_lambda = nnls(L.T, -np.linalg.solve(L, prepared["gap"]))[0]
            ref_dx = prepared["B"] @ ref_lambda
            reference_error = float(np.max(np.abs(dx.ravel() - ref_dx)))
            np.testing.assert_allclose(dx.ravel(), ref_dx, rtol=1e-7, atol=1e-15)
        # Rotating the same physical configuration must rotate its response.
        theta = 0.37
        Q = np.array([[np.cos(theta), -np.sin(theta), 0], [np.sin(theta), np.cos(theta), 0], [0, 0, 1.0]])
        centered_ref = reference - reference.mean(0)
        centered_free = free - reference.mean(0)
        rotated = solve(
            prepare(centered_ref @ Q.T, centered_free @ Q.T, radius, p.viscosity_Pa_s), proto["passes"]
        )
        covariance = float(np.max(np.abs(rotated["displacement_m"] - dx @ Q.T)))
        assert covariance < 1e-15
        count = len(reference)
        # Reverse IDs tests order sensitivity of this deliberately bounded solver.
        reverse = solve(prepare(reference[::-1], free[::-1], radius, p.viscosity_Pa_s), proto["passes"])
        permutation = float(np.max(np.abs(reverse["displacement_m"][::-1] - dx)))
        total = prepared["setup_wall_s"] + result["host_wall_s"]
        original = baseline["setup_wall_s"] + baseline["host_wall_s"]
        rows[name] = dict(
            n=count,
            reference_kind="real saved F7 initial32 + actual free RPY 10us prediction"
            if name == "actual32"
            else "explicit synthetic mathematical contact test",
            overlap_before_m=overlap(free, radius),
            overlap_PBD_m=overlap(baseline["positions_m"], radius),
            overlap_candidate_m=overlap(free + dx, radius),
            projected_residual_m=result["trace"][-1, 0],
            refined_projected_residual_m=fine["trace"][-1, 0],
            max_increment_bound_m=np.sqrt(3) * result["trace"][:, 1].max(),
            mean_PBD_correction_m=(baseline["positions_m"] - free).mean(0),
            mean_candidate_correction_m=dx.mean(0),
            relative_correction_rms_m=np.sqrt(np.mean(np.sum((dx - dx.mean(0)) ** 2, axis=1))),
            correction_difference_rms_m=np.sqrt(
                np.mean(np.sum((free + dx - baseline["positions_m"]) ** 2, axis=1))
            ),
            nonnegative_multipliers=True,
            net_contact_impulse_N_s=impulse.sum(0),
            relative_net_torque=np.linalg.norm(torque) / torque_scale,
            CPU_NNLS_max_difference_m=reference_error,
            rotation_covariance_error_m=covariance,
            ID_permutation_difference_m=permutation,
            candidate_setup_wall_s=prepared["setup_wall_s"],
            candidate_solver_wall_s=result["host_wall_s"],
            candidate_GPU_launch_span_s=result["gpu_launch_span_s"],
            PBD_setup_wall_s=baseline["setup_wall_s"],
            PBD_solver_wall_s=baseline["host_wall_s"],
            PBD_GPU_launch_span_s=baseline["gpu_launch_span_s"],
            setup_inclusive_ratio=total / original,
            force_law="dx = M_RPY J^T impulse; central normal equal/opposite impulses; group mean is not forced to zero",
        )
        rows[name]["numerical_small_state_pass"] = bool(
            rows[name]["overlap_candidate_m"] <= proto["max_overlap_m"]
            and rows[name]["projected_residual_m"] <= proto["max_projected_residual_m"]
            and rows[name]["max_increment_bound_m"] <= proto["max_correction_increment_m"]
        )
        for key, value in dict(
            reference_m=reference,
            free_m=free,
            PBD_m=baseline["positions_m"],
            candidate_m=free + dx,
            trace=result["trace"],
            refined_trace=fine["trace"],
            lambda_m=lam,
            impulse_N_s=impulse,
        ).items():
            arrays[name + "_" + key] = value
    # Separated pair with identical translation must require no correction.
    separate = two * 1.5
    translated = separate + np.array([3e-9, -2e-9, 1e-9])
    zero = solve(prepare(separate, translated, radius, p.viscosity_Pa_s), proto["passes"])
    assert np.array_equal(zero["displacement_m"], np.zeros((2, 3)))
    admitted = (
        all(r["numerical_small_state_pass"] for r in rows.values())
        and rows["actual32"]["setup_inclusive_ratio"] <= proto["max_setup_inclusive_time_ratio_to_PBD"]
    )
    np.savez_compressed(out / "contact_states.npz", **arrays)
    write(
        out / "CHECKS.json",
        dict(
            status="PASS",
            rows=rows,
            zero_contact_translation_pass=True,
            production_admitted=admitted,
            classification="SMALL_STATE_CANDIDATE" if admitted else "NOT_ADMITTED",
            process_wall_s=time.perf_counter() - begun,
            reference="https://arxiv.org/abs/1909.06623v2",
            scope="Frozen free-space RPY with unilateral normal impulses. 16 ordered coordinate sweeps,64-sweep diagnostic only. CPU dense mobility/column assembly plus GPU contact updates, no full constraint matrix in iteration. No wall/lubrication or integrated navigation validation. CUDA spans include launch scheduling; not pure active kernel cost.",
        ),
    )
    print("Contact checks PASS; production admitted:", admitted, flush=True)


if __name__ == "__main__":
    main()
