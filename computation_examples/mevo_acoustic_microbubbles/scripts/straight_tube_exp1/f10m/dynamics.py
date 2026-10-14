"""Actual no-field OFF GPU trajectory in confirmed 32.5 mm lumen, no fabricated ON."""

import sys
import time
import numpy as np
import h5py
import warp as wp
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.metrics import trajectory_metrics
from mevo_mb.fast_sim.straight_tube_engine import cylinder_project, GPUZeroFlow
from mevo_mb.coupling.pdms_block import apparatus
from mevo_mb.fast_sim.pdms_engine import PDMSMainBoreSimulator, reject_port_neighborhood
from mevo_mb.fast_sim.flow_gpu import QueryBuffers
from .common import BASE, configs, write, read, sha, identity


def main():
    c = configs()
    b = c["protocol"]["bubbles"]
    lumen, _ = apparatus(c["geometry"])
    out = BASE / "dynamics"
    out.mkdir(exist_ok=False)
    opened = []

    def hook(event, args):
        if event == "open" and isinstance(args[0], str):
            f = args[0]
            if (
                f.endswith((".vtu", ".vtp"))
                or ("results/fast_pzt" in f and f.endswith((".h5", ".npz")))
                or ("/f10l/" in f and f.endswith((".h5", ".npz")))
            ):
                raise RuntimeError("Historical data access forbidden: " + f)
            if f.endswith((".h5", ".npz")):
                opened.append(f)

    sys.addaudithook(hook)
    p = F2Physics(radius_m=b["radius_m"], frequency_Hz=b["frequency_Hz"])
    path = BASE / "inputs/response.h5"
    assert sha(path) == b["response_sha256"]
    response = ResponseSection(path, p)
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    start = time.perf_counter()
    sim = PDMSMainBoreSimulator(
        initial,
        p,
        response,
        None,
        tube=lumen,
        top_port=c["geometry"]["top_port"],
        acoustic_enabled=False,
        graph_chunk=10,
    )
    # No basis object at all: analytic exact zero with original full-particle engine.
    result = sim.run(
        b["duration_s"],
        b["dt_s"],
        [dict(time_s=0.0, pressure_Pa=0.0)],
        sample_interval_s=b["sample_interval_s"],
        wall_budget_s=60,
    )
    summary = result.pop("summary")
    x = result["positions_m"]
    assert summary["status"] == "COMPLETE"
    np.testing.assert_array_equal(x, np.broadcast_to(initial, x.shape))
    assert lumen.contains(x, p.radius_m).all()
    assert (
        np.count_nonzero(result["flow"]["velocity_m_s"])
        == np.count_nonzero(result["flow"]["gradient_s_inv"])
        == 0
    )
    sim.fields()
    f1 = sim.force.numpy().copy()
    sim.pairs()
    f2 = sim.force.numpy() - f1
    assert np.count_nonzero(f1) == np.count_nonzero(f2) == np.count_nonzero(sim.pressure.numpy()) == 0
    metrics, events = trajectory_metrics(x, p.radius_m, result["time_s"])
    with h5py.File(out / "OFF_trajectory.h5", "w") as h:
        for name, v in result.items():
            if isinstance(v, np.ndarray):
                h.create_dataset(name, data=v, compression="gzip")
        h.create_dataset("particle_id", data=np.arange(b["count"]))
        h.create_dataset("radius_m", data=p.radius_m)
        h.create_dataset("primary_force_N", data=np.zeros_like(x), compression="gzip")
        h.create_dataset("secondary_force_N", data=np.zeros_like(x), compression="gzip")
        h.create_dataset("pressure_complex_Pa", data=np.zeros(x.shape[:-1], complex), compression="gzip")
        h.attrs["acoustic_basis"] = "NONE; analytic zero; apparatus field NOT_RUN"
        h.attrs["geometry"] = "32.5 mm centered PDMS lumen, ID 1.6 mm"
        h.attrs["execution_commit"] = identity()["execution_commit"]
    np.savez_compressed(out / "metrics.npz", **metrics)
    checks = {}

    def projection(q):
        xx = wp.array(q, dtype=wp.vec3d, device=sim.device)
        stats = wp.zeros((len(q), 5), dtype=wp.float64, device=sim.device)
        err = wp.zeros(1, dtype=int, device=sim.device)
        wp.launch(
            cylinder_project,
            dim=len(q),
            inputs=[
                xx,
                wp.vec3d(*lumen.center),
                wp.mat33d(*lumen.rotation.ravel()),
                lumen.radius_m,
                lumen.length_m / 2,
                p.radius_m,
                int(lumen.ends[0] == "closed"),
                int(lumen.ends[1] == "closed"),
                stats,
                err,
            ],
            device=sim.device,
        )
        return xx.numpy(), stats.numpy(), int(err.numpy()[0])

    xx, ss, err = projection([[0, lumen.radius_m - p.radius_m + 0.05 * p.radius_m, 0]])
    checks["small_radial_correction_within_original_budget"] = err == 0 and np.isclose(
        xx[0, 1], lumen.radius_m - p.radius_m, atol=1e-15
    )
    checks["wall_over_0p1R_rejected"] = bool(
        projection([[0, lumen.radius_m - p.radius_m + 0.2 * p.radius_m, 0]])[2] & 1024
    )
    checks["open_x_minus_endpoint_fail_stop"] = bool(projection([[-0.01625 - 1e-9, 0, 0]])[2] & 32768)
    xx, ss, err = projection([[0.01625 - p.radius_m + 0.05 * p.radius_m, 0, 0]])
    checks["confirmed_closed_x_plus_finite_radius"] = err == 0 and np.isclose(
        xx[0, 0], 0.01625 - p.radius_m, atol=1e-15
    )
    guard_x = wp.array([[0.01085, 0, 0]], dtype=wp.vec3d, device=sim.device)
    guard_error = wp.zeros(1, dtype=int, device=sim.device)
    wp.launch(
        reject_port_neighborhood,
        dim=1,
        inputs=[
            guard_x,
            wp.vec3d(*lumen.center),
            wp.mat33d(*lumen.rotation.ravel()),
            0.01085,
            0.00075 + p.radius_m,
            guard_error,
        ],
        device=sim.device,
    )
    checks["junction_guard_prevents_false_cylindrical_wall"] = bool(guard_error.numpy()[0] & 65536)
    checks["old_20mm_boundary_not_retained"] = projection([[0.015, 0, 0]])[2] == 0
    xx = wp.array(np.tile([0.015, 0.0001, 0], (4096, 1)), dtype=wp.vec3d, device=sim.device)
    buff = QueryBuffers(4096, sim.device)
    err = wp.zeros(1, dtype=int, device=sim.device)
    GPUZeroFlow().launch(
        xx,
        wp.ones(4096, dtype=int, device=sim.device),
        wp.zeros(4096, dtype=int, device=sim.device),
        buff,
        err,
    )
    checks["4096_GPU_flow_and_gradient_exact_zero"] = (
        np.count_nonzero(buff.velocity.numpy()) == np.count_nonzero(buff.gradient.numpy()) == 0
    )
    checks["OFF_bitwise_no_drift"] = True
    checks["all_finite_radius_inside"] = True
    checks["no_CFD_or_old_basis_loaded"] = True
    write(
        out / "SUMMARY.json",
        dict(
            summary,
            solver_geometry_version="F10M_PDMS_CLOSED_X_PLUS_TOP_OPEN_MAIN_BORE_M0",
            junction_scope="Top port included in geometry/materials; main-bore motion only; junction neighborhood rejects, no false reflection",
            physical_time_s=float(result["time_s"][-1]),
            maximum_displacement_m=0.0,
            whole_case_wall_s=time.perf_counter() - start,
            Rg_start_end_m=metrics["global_rg_m"][[0, -1], 0],
            components_start_end=metrics["cluster_count"][[0, -1], 0],
            checks=checks,
            identity=identity(),
        ),
    )
    write(
        out / "NO_CFD_ACCESS.json", dict(opened=sorted(set(opened)), CFD_access_count=0, basis_loaded=False)
    )
    write(
        out / "ON_NOT_RUN.json",
        dict(
            status="NOT_RUN",
            reasons=read(BASE / "analysis/ACOUSTIC_ADMISSION.json")["reasons"],
            training_ready=False,
        ),
    )
    write(out / "CLUSTER_EVENTS.json", events)
    assert all(checks.values())
    print(dict(status="PASS", OFF="COMPLETE", ON="NOT_RUN", checks=checks), flush=True)


if __name__ == "__main__":
    main()
