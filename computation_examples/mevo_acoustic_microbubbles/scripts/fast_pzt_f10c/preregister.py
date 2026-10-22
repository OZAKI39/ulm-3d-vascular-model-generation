"""Offline baseline audit and frozen numerical gates, before F10C experiments."""

import subprocess
import h5py
import numpy as np
from mevo_mb.fast_sim.morphodynamics import contrast, shape
from .common import ROOT, BASE, F10A, js, sha, write


def main():
    checked = {}
    for stage in ["f10a", "f10b"]:
        folder = ROOT / "results/fast_pzt" / stage
        manifest = js(folder / "ARTIFACT_MANIFEST.json")["files"]
        for rel, meta in manifest.items():
            assert sha(folder / rel) == meta["sha256"], (stage, rel)
        checked[stage] = len(manifest)
    for rel, meta in js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json").items():
        assert sha(ROOT / "results/fast_pzt/f5" / rel) == meta["sha256"], rel
    sensitivity = {}
    schemas = {}
    for key in ["OFF", "S_1000", "O_1000"]:
        paths = [F10A / phase / key / "trajectory.h5" for phase in ["extended_10_v1", "extended_5_v1"]]
        arrays = []
        for path in paths:
            with h5py.File(path) as h:
                arrays.append(h["positions_m"][:, 0])
                schemas[str(path.relative_to(ROOT))] = dict(
                    sha256=sha(path), keys=list(h.keys()), shape=h["positions_m"].shape
                )
        a, b = arrays
        rows = [contrast(x, y, np.arange(32)) for x, y in zip(a, b)]
        sensitivity[key] = {
            k: float(np.max([r[k] for r in rows]))
            for k in ["centroid_distance_m", "rg_relative", "tensor_relative", "aligned_internal_rmsd_m"]
        }
        sensitivity[key]["max_particle_error_m"] = float(np.linalg.norm(a - b, axis=-1).max())
        aspect = lambda x: np.sqrt(shape(x)[2][0] / shape(x)[2][-1])
        sensitivity[key]["max_aspect_relative"] = float(
            max(abs(aspect(y) / aspect(x) - 1) for x, y in zip(a, b))
        )
    protocol = dict(
        baseline_commit="81fe1206d623f04ea7e3d3a78021251890185c9e",
        duration_s=0.01,
        dt_s=1e-5,
        sample_s=5e-5,
        N=32,
        batch_sizes=[1, 4, 16],
        conditional_batch_size=64,
        batch_equality="bitwise all sampled state, forces, velocities, diagnostics and sorted events; zero tolerance",
        batch_controls="one common event schedule per family; asynchronous per-member controls unsupported",
        initial_variants="BENCHMARK_ONLY: translations along Y of 0, +0.2, -0.2, +0.4 micrometres, repeated as needed; no independent training samples",
        failure_policy="global device error aborts whole batch; all members rejected, individual diagnostic causes retained; no unaffected member marked PASS",
        modes={"Q0": [1, 1], "Q1": [2, 1], "Q2": [1, 2], "Q3": [2, 2]},
        multirate_group_gates=dict(max_centroid_distance_m=1e-8, max_paired_z_error_m=1e-8),
        multirate_morphology_gates=dict(
            max_particle_error_m=1.78e-7,
            max_aligned_rmsd_m=8.9e-8,
            max_rg_relative=0.01,
            max_tensor_relative=0.02,
            max_aspect_relative=0.02,
            largest_fraction_difference=0,
            edge_events_exact=True,
        ),
        rationale="Group 10 nm is above the observed 10/5 us centroid sensitivity but below one tenth of the paired Z signal. Morphology requires max ID error <=0.1R, aligned RMSD <=0.05R and 1-2% shape differences: deliberately stricter than unresolved F10A particle sensitivity. Exact edge-event agreement is required for morphology truth; a centroid pass alone cannot admit full morphology.",
        contact_budgets=dict(max_residual_over_R=0.01, max_single_pass_over_R=0.1),
        q4_condition="Only after q2 group AND morphology PASS; no full Cartesian matrix",
        engineering_speed_gate=2.0,
        training_ready=False,
        physics_scope="EXPLORATORY",
        F10A_dt_sensitivity=sensitivity,
        profiling="Warp event-instrumented per-launch activities are profiler spans, not CUPTI pure-kernel durations; unmeasured pure kernel time remains null. Unprofiled wall times are separate.",
    )
    target = ROOT / "configs/fast_pzt/f10c/protocol.json"
    assert not target.exists()
    write(target, protocol)
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            status="PASS",
            historical_files_verified=checked,
            input_files_verified=13,
            hdf5=schemas,
            worktrees=subprocess.check_output(["git", "worktree", "list"], text=True),
            f10b_gate=js(ROOT / "results/fast_pzt/f10b/analysis/FINAL_GATE.json"),
        ),
    )
    print(sensitivity)


if __name__ == "__main__":
    main()
