"""Small GPU wall-oracle comparison + actual short OFF; no apparatus acoustics."""

import os

os.environ["KWAVE_F10L_PYTHON_ONLY"] = "1"
import sys
import time
import traceback
import numpy as np
import h5py
import warp as wp
from mevo_mb.coupling.pdms_cavity import apparatus
from mevo_mb.fast_sim.pdms_cavity_engine import PDMSCavitySimulator, query_walls, cavity_project
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from mevo_mb.fast_sim.metrics import trajectory_metrics
from mevo_mb.acoustics.kwave_backend import KWaveCW, GPUTelemetry
from .common import BASE, configs, read, write, sha, identity


def source_input_check(out):
    """Four tiny independent numerical discs: CPU source weights only, no wave solve."""
    t = time.perf_counter()
    folder = out / "source_input"
    folder.mkdir()
    elements = []
    for name, pos, n in [
        ("Y-", [0, -0.00025, 0], [0, 1, 0]),
        ("Y+", [0, 0.00025, 0], [0, -1, 0]),
        ("X-", [-0.00025, 0, 0], [1, 0, 0]),
        ("X+", [0.00025, 0, 0], [-1, 0, 0]),
    ]:
        elements.append(
            dict(
                id=name,
                position_m=pos,
                normal=n,
                geometry="disc",
                size_m=[0.0002],
                amplitude_Pa=0.0,
                phase_rad=0.0,
                enabled=True,
            )
        )
    cfg = dict(
        grid=dict(shape=[112] * 3, spacing_m=[25e-6] * 3),
        medium=dict(sound_speed_m_s=1500, density_kg_m3=1000),
        pzt=elements,
        acoustics=dict(
            frequency_Hz=5e6,
            CFL=0.2,
            transit_multiples=1,
            ramp_periods=4,
            recording_periods=2,
            roi_half_width_m=[0.0003] * 3,
            bli_tolerance=0.01,
            upsampling_rate=4,
            pml_size=16,
            pml_inside=False,
            pml_alpha=2,
            backend="python",
            device="gpu",
            BonA=None,
        ),
    )
    write(folder / "CONFIG.json", cfg)
    adapter = KWaveCW(cfg, folder)
    occupied = np.argwhere(adapter.source_mask)
    clear = np.minimum(occupied.min(0), np.array([112] * 3) - 1 - occupied.max(0))
    u = np.array([100, 200j, -300, 0], complex)
    field_of_sources = adapter.weights.T @ u
    manual = sum(adapter.weights[j] * u[j] for j in range(4))
    gram = adapter.weights.astype(float) @ adapter.weights.astype(float).T
    assert min(clear) >= 1 and np.linalg.eigvalsh(gram).min() > 1e-6
    np.testing.assert_allclose(field_of_sources, manual, atol=1e-12)
    rec = dict(
        status="PASS",
        scope="CPU_SOURCE_WIRING_ONLY; 0.2mm numerical test discs, NO_PRESSURE_FIELD_SOLVE",
        source_count=4,
        independent_weight_rank=int(np.linalg.matrix_rank(gram)),
        edge_clearance_cells=clear,
        weights_sha256=sha(folder / "source_weights.npz"),
        wall_s=time.perf_counter() - t,
        experimental_prediction=False,
        actual_wave_solves=0,
        identity=identity(),
    )
    write(folder / "CHECK.json", rec)
    return rec


def main(out):
    started = time.perf_counter()
    c = configs()
    assert c == read(BASE / "CONFIG_SNAPSHOT.json")
    cavity, _ = apparatus(c["geometry"])
    b = c["protocol"]["bubbles"]
    physics = F2Physics(radius_m=b["radius_m"], frequency_Hz=b["frequency_Hz"])
    source = source_input_check(out)
    opened = []

    def hook(event, args):
        if event == "open" and isinstance(args[0], str):
            p = args[0]
            if p.endswith((".vtu", ".vtp")) or (
                "/results/" in p and p.endswith((".h5", ".npz")) and "/f10m_r3/" not in p
            ):
                raise RuntimeError("Unexpected historical scientific input: " + p)
            if p.endswith((".h5", ".npz")):
                opened.append(p)

    sys.addaudithook(hook)
    response_path = BASE / "inputs/response.h5"
    assert sha(response_path) == b["response_sha256"]
    response = ResponseSection(response_path, physics)
    write(out / "RESPONSE_SCOPE.json", response.metadata)
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    reference = np.load(BASE / "wall_reference.npz")
    wp.init()
    device = "cuda:0"
    checks = {}
    with GPUTelemetry() as telemetry:
        compile_start = time.perf_counter()
        wp.load_module(module=query_walls.module, device=device)
        wp.synchronize_device(device)
        compile_wall = time.perf_counter() - compile_start
        probes = reference["positions_m"]
        x = wp.array(probes, dtype=wp.vec3d, device=device)
        distances = wp.zeros(len(probes), dtype=wp.float64, device=device)
        normals = wp.zeros(len(probes), dtype=wp.vec3d, device=device)
        closest = wp.zeros(len(probes), dtype=wp.vec3d, device=device)
        t = time.perf_counter()
        wp.launch(
            query_walls,
            dim=len(probes),
            inputs=[x, cavity.R, cavity.a, cavity.L, cavity.H, cavity.ports[1], distances, normals, closest],
            device=device,
        )
        d, n, v = distances.numpy(), normals.numpy(), closest.numpy()
        query_wall = time.perf_counter() - t
        delta = float(np.max(abs(d - reference["signed_distance_m"])))
        checks["GPU_exposed_union_distance_matches_dense_CPU"] = (
            delta < c["protocol"]["wall_query_tolerance_m"]
        )
        checks["GPU_unit_normals_finite"] = np.isfinite(n).all() and np.allclose(
            np.linalg.norm(n, axis=1), 1, atol=1e-9
        )
        np.savez_compressed(
            out / "wall_comparison.npz",
            positions_m=probes,
            cpu_distance_m=reference["signed_distance_m"],
            gpu_distance_m=d,
            gpu_normal=n,
            gpu_closest_wall_m=v,
        )

        def project(q):
            xx = wp.array(np.atleast_2d(q), dtype=wp.vec3d, device=device)
            stats = wp.zeros((len(xx), 5), dtype=wp.float64, device=device)
            stats.assign(np.c_[np.full(len(xx), np.inf), np.zeros((len(xx), 4))])
            err = wp.zeros(1, dtype=int, device=device)
            wp.launch(
                cavity_project,
                dim=len(xx),
                inputs=[
                    xx,
                    wp.vec3d(*cavity.lumen.center),
                    wp.mat33d(*cavity.lumen.rotation.ravel()),
                    cavity.R,
                    cavity.a,
                    cavity.L,
                    cavity.H,
                    cavity.ports[1],
                    physics.radius_m,
                    cavity.ceiling,
                    stats,
                    err,
                ],
                device=device,
            )
            return xx.numpy(), stats.numpy(), int(err.numpy()[0])

        pos, st, err = project(reference["projection_input"])
        correction_delta = float(np.max(np.linalg.norm(pos - reference["projection_output"], axis=1)))
        checks["GPU_finite_radius_projection_ends_side_ports_junction"] = (
            err == 0 and correction_delta < 5e-10
        )
        checks["GPU_corrected_spheres_inside_actual_union"] = all(
            cavity.query(q, physics.radius_m)["valid"] for q in pos
        )
        for sign in [-1, 1]:
            q = np.array([sign * (cavity.L - physics.radius_m + 0.2 * physics.radius_m), 0, 0])
            pos, st, err = project(q)
            checks["X_" + str(sign) + "_over_budget_rejected_no_escape"] = (
                err & 1024
            ) != 0 and np.array_equal(pos[0], q)
        q = np.array([cavity.ports[0], 0, cavity.ceiling - physics.radius_m / 2])
        pos, st, err = project(q)
        checks["GPU_unknown_top_scope_stops_without_projection"] = (err & 131072) != 0 and np.array_equal(
            pos[0], q
        )
        path = np.array([[x, 0, z] for x in cavity.ports for z in np.linspace(0.0006, 0.0012, 13)])
        pos, st, err = project(path)
        checks["GPU_both_junction_paths_not_false_cylinder_walls"] = err == 0 and np.array_equal(pos, path)
        write(
            out / "WALL_CHECKS.json",
            dict(
                status="PASS" if all(checks.values()) else "FAIL",
                checks=checks,
                query_points=len(probes),
                distance_max_error_m=delta,
                projection_max_error_m=correction_delta,
                compile_s=compile_wall,
                query_wall_s=query_wall,
                scope="geometric finite-radius validity; no near-wall hydrodynamic calibration",
            ),
        )
        assert all(checks.values()), checks
        sim = PDMSCavitySimulator(
            initial, physics, response, None, cavity=cavity, acoustic_enabled=False, graph_chunk=10
        )
        for mode in ["three_controls", "nonzero_drive"]:
            try:
                event = dict(time_s=0.0, amplitude_Pa=[0] * 3 if mode == "three_controls" else [1] * 4)
                sim.control(event)
                checks["reject_" + mode] = False
            except ValueError:
                checks["reject_" + mode] = True
        trajectory = sim.run(
            b["duration_s"],
            b["dt_s"],
            [dict(time_s=0.0, pressure_Pa=0.0)],
            sample_interval_s=b["sample_interval_s"],
            wall_budget_s=30,
        )
        summary = trajectory.pop("summary")
        xx = trajectory["positions_m"]
        checks["actual_OFF_complete"] = summary["status"] == "COMPLETE"
        checks["actual_OFF_bitwise_static"] = np.array_equal(xx, np.broadcast_to(initial, xx.shape))
        checks["actual_OFF_zero_flow_gradient"] = (
            np.count_nonzero(trajectory["flow"]["velocity_m_s"])
            == np.count_nonzero(trajectory["flow"]["gradient_s_inv"])
            == 0
        )
        checks["actual_OFF_no_wall_correction"] = summary["quality"]["maximum_wall_correction_m"] == 0
        sim.fields()
        f1 = sim.force.numpy().copy()
        sim.pairs()
        f2 = sim.force.numpy() - f1
        checks["actual_OFF_zero_pressure_primary_secondary"] = (
            np.count_nonzero(f1) == np.count_nonzero(f2) == np.count_nonzero(sim.pressure.numpy()) == 0
        )
        metrics, events = trajectory_metrics(xx, physics.radius_m, trajectory["time_s"])
        with h5py.File(out / "OFF_trajectory.h5", "w") as h:
            for name, val in trajectory.items():
                if isinstance(val, np.ndarray):
                    h.create_dataset(name, data=val, compression="gzip")
            h.create_dataset("particle_id", data=np.arange(len(initial)))
            h.create_dataset("radius_m", data=physics.radius_m)
            h.create_dataset("primary_force_N", data=np.broadcast_to(f1, xx.shape), compression="gzip")
            h.create_dataset("secondary_force_N", data=np.broadcast_to(f2, xx.shape), compression="gzip")
            h.attrs["scope"] = "DEMO_INTERIOR_WETTING_ONLY; no free surface physics"
            h.attrs["execution_commit"] = identity()["execution_commit"]
            h.attrs["end_conditions"] = "closed,closed"
            h.attrs["acoustic_basis"] = "NONE; actual OFF only"
        np.savez_compressed(out / "metrics.npz", **metrics)
        write(out / "CLUSTER_EVENTS.json", events)
        write(out / "OFF_SUMMARY.json", summary)
    write(
        out / "GPU_TELEMETRY.json",
        dict(
            samples=telemetry.samples,
            errors=telemetry.errors,
            attribution="whole-device sampled peak, exclusive compute-app admission; includes compilation and small checks",
        ),
    )
    write(out / "NO_CFD_ACCESS.json", dict(CFD_access_count=0, opened_scientific_inputs=sorted(set(opened))))
    rec = dict(
        status="PASS" if all(checks.values()) else "FAIL",
        checks=checks,
        whole_validation_wall_s=time.perf_counter() - started,
        source_input_wall_s=source["wall_s"],
        sampled_device_peak_MiB=max(s["used_MiB"] for s in telemetry.samples),
        OFF_integrator_wall_s=summary.get("wall_s"),
        wall_distance_max_error_m=delta,
        projection_max_error_m=correction_delta,
        physical_duration_s=b["duration_s"],
        bubble_count=b["count"],
        identity=identity(),
        formal_ON="NOT_RUN",
        actual_wave_solves=0,
    )
    write(out / "VALIDATION.json", rec)
    assert all(checks.values()), checks
    write(
        BASE / "ACCEPTED_VALIDATION.json",
        dict(
            path=str(out.relative_to(BASE)),
            validation_sha256=sha(out / "VALIDATION.json"),
            identity=identity(),
        ),
    )
    print(rec, flush=True)


if __name__ == "__main__":
    name = os.environ.get("F10M_R3_ATTEMPT", "validation_v1")
    if not name.replace("_", "").isalnum():
        raise ValueError("Simple unique attempt label required")
    out = BASE / name
    out.mkdir(exist_ok=False)
    try:
        main(out)
    except Exception:
        write(out / "FAILURE.json", dict(traceback=traceback.format_exc(), identity=identity()))
        raise
