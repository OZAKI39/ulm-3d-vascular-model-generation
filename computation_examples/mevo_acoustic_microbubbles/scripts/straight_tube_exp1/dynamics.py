"""OFF / three-source short trajectories and actual same-ID observables."""

import sys
import time
import importlib.metadata
import numpy as np
from mevo_mb.coupling.straight_tube import apparatus
from mevo_mb.fast_sim.models import Basis, ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference
from mevo_mb.fast_sim.metrics import trajectory_metrics
from mevo_mb.fast_sim.straight_tube_engine import StraightTubeSimulator
from mevo_mb.fast_sim.consistent_field import query_numpy
from .common import ROOT, BASE, configs, read, write, sha, identity, validate_saved_acoustic_inputs


def run():
    validate_saved_acoustic_inputs()
    accesses = []

    def audit(event, args):
        if event == "open" and isinstance(args[0], str):
            path = args[0]
            if path.endswith((".vtu", ".vtp")) or (
                "results/fast_pzt" in path and path.endswith((".h5", ".npz"))
            ):
                raise RuntimeError("F10L must not load historical Y geometry, basis or CFD: " + path)
            if path.endswith((".h5", ".npz")):
                accesses.append(path)

    sys.addaudithook(audit)
    c = configs()
    b = c["bubbles"]
    a = c["acoustic"]
    tube, _ = apparatus(c["demo"]["geometry"], a["acoustics"]["frequency_Hz"])
    p = F2Physics(
        radius_m=b["radius_m"],
        frequency_Hz=a["acoustics"]["frequency_Hz"],
        density_kg_m3=a["medium"]["density_kg_m3"],
        viscosity_Pa_s=a["medium"]["viscosity_Pa_s"],
        radial_effective_viscosity_Pa_s=a["medium"]["radial_effective_viscosity_Pa_s"],
    )
    assert sha(ROOT / b["response_path"]) == b["response_sha256"]
    response = ResponseSection(ROOT / b["response_path"], p)
    basis = Basis(BASE / "acoustic/basis.h5", p.frequency_Hz)
    assert read(BASE / "acoustic/ACOUSTIC_GATE.json")["status"] == "PASS"
    initial = np.load(BASE / "initial_state.npz")["positions_m"]
    folder = BASE / "dynamics"
    folder.mkdir(exist_ok=False)
    write(
        folder / "ENVIRONMENT.json",
        {
            name: importlib.metadata.version(name)
            for name in ["warp-lang", "numpy", "scipy", "h5py", "k-wave-python", "cupy-cuda12x"]
        },
    )
    cases = {}
    for name, on in [("OFF", False), ("THREE_PZT_ON", True)]:
        started = time.perf_counter()
        out = folder / name
        out.mkdir()
        event = dict(a["demo_control"])
        if not on:
            event.update(amplitude_Pa=[0.0] * 3, enabled=[False] * 3)
        sim = StraightTubeSimulator(
            initial, p, response, basis, tube=tube, acoustic_enabled=on, graph_chunk=10
        )
        result = sim.run(
            b["duration_s"],
            b["dt_s"],
            [event],
            sample_interval_s=b["sample_interval_s"],
            wall_budget_s=b["wall_budget_s"],
        )
        summary = result.pop("summary")
        np.savez_compressed(
            out / "trajectory.npz",
            **{k: v for k, v in result.items() if isinstance(v, np.ndarray)},
            particle_id=np.arange(b["count"]),
            radius_m=p.radius_m,
        )
        write(
            out / "SUMMARY.json",
            dict(
                summary,
                identity=identity(),
                basis_sha256=sha(BASE / "acoustic/basis.h5"),
                response_sha256=sha(ROOT / b["response_path"]),
            ),
        )
        if summary["status"] != "COMPLETE":
            raise RuntimeError(name + " failed; retained trajectory prefix and summary")
        x = result["positions_m"]
        assert tube.contains(x, p.radius_m).all()
        assert np.all(result["flow"]["velocity_m_s"] == 0)
        assert np.all(result["flow"]["gradient_s_inv"] == 0)
        if not on:
            assert np.array_equal(x, np.broadcast_to(initial, x.shape))
        metrics, events = trajectory_metrics(x, p.radius_m, result["time_s"])
        np.savez_compressed(out / "metrics.npz", **metrics)
        write(out / "CLUSTER_EVENTS.json", events)
        pressures = []
        gradients = []
        primary = []
        secondary = []
        reference_errors = []
        # Re-evaluate unchanged GPU forces at each actual saved frame. These
        # are frame observables, not invented substep force histories.
        sim.control(event)
        for xx in x[:, 0]:
            sim.x.assign(xx)
            sim.fields()
            f1 = sim.force.numpy().copy()
            pressure = sim.pressure.numpy()
            pp = pressure[:, 0] + 1j * pressure[:, 1]
            gg = sim.consistent_gradients.numpy()
            gg = gg[..., 0] + 1j * gg[..., 1]
            sim.pairs()
            f2 = sim.force.numpy() - f1
            assert int(sim.error.numpy()[0]) == 0
            pressures.append(pp)
            gradients.append(gg if on else np.zeros_like(gg))
            primary.append(f1)
            secondary.append(f2)
        # One saved frame cross-checks CPU table, gradient and pair reduction.
        xx = x[0, 0]
        pp = np.array(pressures)[0]
        hh = response.query(pp)
        if on:
            cp, cg = query_numpy(sim.last_field.pressure_Pa, basis.origin, basis.spacing, xx)
            f1ref = -0.5 * np.real(hh[:, 0, None] * cg.conj())
            f2ref = pair_reference(xx, hh, pp, p, sim.chi_zero)
            for got, expected in [(pressures[0], cp), (primary[0], f1ref), (secondary[0], f2ref)]:
                err = float(
                    np.linalg.norm(np.asarray(got) - expected) / max(np.linalg.norm(expected), 1e-300)
                )
                reference_errors.append(err)
            assert max(reference_errors) < 1e-9
        np.savez_compressed(
            out / "frame_observables.npz",
            time_s=result["time_s"],
            pressure_complex_Pa=pressures,
            gradient_complex_Pa_m=gradients,
            primary_N=primary,
            secondary_N=secondary,
            background_velocity_m_s=np.zeros_like(x[:, 0]),
            wall_gap_m=tube.gaps(x, p.radius_m),
        )
        centroid = metrics["global_centroid_m"][:, 0]
        eig = metrics["global_eigenvalues_m2"][:, 0]
        aspect = np.sqrt(eig[:, 0] / eig[:, -1])
        cases[name] = dict(
            status=summary["status"],
            physical_time_s=float(result["time_s"][-1]),
            steps=summary["actual_substep_count"],
            recorded_frames=len(x),
            integration_wall_s=summary["wall_s"],
            whole_case_wall_s=time.perf_counter() - started,
            maximum_same_id_displacement_m=float(np.linalg.norm(x[:, 0] - initial, axis=-1).max()),
            centroid_displacement_m=(centroid[-1] - centroid[0]).tolist(),
            Rg_start_end_m=metrics["global_rg_m"][[0, -1], 0].tolist(),
            aspect_start_end=aspect[[0, -1]].tolist(),
            components_start_end=metrics["cluster_count"][[0, -1], 0].tolist(),
            pressure_peak_range_Pa=[float(np.abs(pressures).min()), float(np.abs(pressures).max())],
            maximum_primary_force_N=float(np.linalg.norm(primary, axis=-1).max()),
            maximum_secondary_force_N=float(np.linalg.norm(secondary, axis=-1).max()),
            quality=summary["quality"],
            CPU_GPU_reference_relative_errors=reference_errors,
            strictly_zero_flow=True,
            inside=True,
            training_ready=False,
        )
        write(out / "OBSERVABLE_SUMMARY.json", cases[name])
    # Limited local control directions at the actual initial cloud center.
    pos = initial.mean(0)[None, :]
    directions = []
    for i in range(3):
        weights = np.eye(3)[i] * a["demo_control"]["amplitude_Pa"][i]
        field, _, _ = basis.combine(weights, [0, 0, 0])
        pp, gg = query_numpy(field.pressure_Pa, basis.origin, basis.spacing, pos)
        vv = response.query(pp)[:, 0]
        force = -0.5 * np.real(vv[:, None] * gg.conj())
        directions.append(
            dict(
                source_id=basis.db.metadata["source_ids"][i],
                position_m=pos[0],
                pressure_complex_Pa=[float(pp[0].real), float(pp[0].imag)],
                grad_real_Pa_m=gg.real[0],
                grad_imag_Pa_m=gg.imag[0],
                primary_N=force[0],
            )
        )
    write(
        folder / "LOCAL_CONTROL_DIAGNOSTIC.json",
        dict(
            scope="Three individual source directions at one point only; no reachability or arbitrary shape claim",
            directions=directions,
        ),
    )
    write(
        folder / "NO_CFD_ACCESS.json",
        dict(
            status="PASS",
            guard="Python open audit rejects old Y assets and every VTK data file",
            data_files_opened=sorted(set(accesses)),
            CFD_access_count=0,
            background_mode="zero",
        ),
    )
    write(
        folder / "DYNAMICS_GATE.json",
        dict(
            status="PASS",
            cases=cases,
            identity=identity(),
            physical_calibration=False,
            nav_and_shape_control="NOT_TESTED",
            training_ready=False,
        ),
    )
    print(cases, flush=True)


if __name__ == "__main__":
    run()
