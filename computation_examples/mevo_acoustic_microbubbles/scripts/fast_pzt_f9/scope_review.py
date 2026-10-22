"""Small postprocessing-only mechanism, cost and RPY scope review."""

import h5py
import numpy as np
from .common import ROOT, BASE, js, sha, write


def main():
    result = js(BASE / "analysis/RESULTS.json")
    mesh_path = ROOT / "results/fast_pzt/f5/inputs/f4/prepared_geometry.npz"
    mesh = np.load(mesh_path)
    faces = mesh["faces"][mesh["labels"] == 1]
    mechanism = {}
    performance = {}
    for key, row in result["cases"].items():
        if "ledger" not in row:
            continue
        path = ROOT / row["path"] / "trajectory.h5"
        with h5py.File(path) as h:
            ev = h["motion_ledger/wall_contact_events"][:]
            faces_used = np.unique(ev[:, 19]).astype(int)
            triangles = {}
            for f in faces_used:
                vertices = mesh["vertices_m"][faces[f]]
                normal = np.cross(vertices[1] - vertices[0], vertices[2] - vertices[0])
                normal /= np.linalg.norm(normal)
                recorded = ev[ev[:, 19] == f, 16:19]
                np.testing.assert_allclose(
                    recorded, np.broadcast_to(normal, recorded.shape), atol=2e-12, rtol=0
                )
                triangles[str(f)] = dict(
                    vertices_m=vertices,
                    outward_normal=normal,
                    normal_source="Original vertex winding, independently checked against recorded repeated production query",
                )
            ids = {}
            for particle in np.unique(ev[:, 3]).astype(int):
                e = ev[ev[:, 3] == particle]
                ids[str(particle)] = dict(
                    first_step_start_s=e[0, 1],
                    last_step_start_s=e[-1, 1],
                    count=len(e),
                    signed_vector_m=e[:, 11:14].sum(0),
                    path_length_m=e[:, 14].sum(),
                )
        controls = js(path.parent / "SUMMARY.json")["controls"]
        mechanism[key] = dict(
            by_ID=ids,
            original_triangles=triangles,
            actual_control_event_times_s=[c["time_s"] for c in controls],
            interpretation="Face normals unchanged. Different histories include changed facet timing and particle arrangement, not new wall laws. Pair opposition does not establish an implementation error.",
        )
        summary = js(path.parent / "SUMMARY.json")
        inherited = (
            ROOT
            / "results/fast_pzt/f8"
            / ("y_half_v1" if key.endswith("_5") else "y_formal_v1")
            / key[:2]
            / "SUMMARY.json"
        )
        performance[key] = dict(
            integration_wall_s=summary["wall_s"],
            gpu_stream_elapsed_s=summary["gpu_stream_elapsed_s"],
            gpu_scope=summary["gpu_time_scope"],
            audit_readback_s=summary["audit"]["transfer_and_reduction_readback_s"],
            ordinary_sampled_download_s=summary["sampled_download_s"],
            output_and_postprocessing_s=summary["output_and_postprocessing_s"],
            hdf5_initial_write_s=summary["hdf5_wall_s"],
            process_body_wall_s=summary["process_body_wall_s"],
            trajectory_bytes=path.stat().st_size,
            original_F8_integration_s=js(inherited)["wall_s"] if inherited.exists() else None,
            original_F8_trajectory_bytes=inherited.with_name("trajectory.h5").stat().st_size
            if inherited.exists()
            else None,
        )
    pair = js(BASE / "wall_admission_v1/CHECKS.json")["performance"]
    off, on = pair["Y_wall_off"], pair["Y_wall_on"]
    write(
        BASE / "analysis/PERFORMANCE.json",
        dict(
            formal=performance,
            same_session_30ms=dict(
                audit_off_s=off["wall_s"],
                audit_on_s=on["wall_s"],
                ratio=on["wall_s"] / off["wall_s"],
                same_trajectory_bitwise=True,
            ),
            cost_scope="Single paired test, not a repeated performance distribution. Integration contains synchronization, normal saved-frame processing and audit readback. CUDA stream span includes host gaps, not exclusive active GPU compute. Historical F8 timings are not same-session benchmarks.",
            bottleneck="Additional diagnostic copies/kernels, ordered 32-particle reduction and repeated read-only wall query; sparse readback is a small recorded fraction. No profiler-derived kernel attribution claimed. Audit is explicitly off by default.",
            science_archive_bytes=js(BASE / "ARCHIVE_VERIFICATION.json")["remote_bytes"],
        ),
    )
    write(
        BASE / "analysis/WALL_MECHANISM.json",
        dict(
            cases=mechanism,
            geometry_sha256=sha(mesh_path),
            time_window_decision="No 2.5us trajectory: Y1 difference grows across multiple 90-110,150-160,190-210,250-260ms episodes. Exact component identity identifies wall-dominated sensitivity; two timesteps do not establish convergence order. Saved 0.5ms frames do not restore the full production caches, diagnostics, clocks and controller state. No validated complete-state restart exists for a common-state short window; no invented restart or new full run.",
            missing="No Y0 10us or Y2 5us exact ledger; their saved F8 trajectories supply positional sensitivity only.",
        ),
    )
    old = ROOT / "results/fast_pzt/f8/contact_v1/CHECKS.json"
    rpy = js(old)
    write(
        BASE / "analysis/RPY_SCOPE.json",
        dict(
            status="REVIEW_COMPLETE_NO_NEW_RUN",
            classification=rpy["classification"],
            sources={
                str(old.relative_to(ROOT)): sha(old),
                "src/mevo_mb/fast_sim/f8_contact.py": sha(ROOT / "src/mevo_mb/fast_sim/f8_contact.py"),
            },
            saved_actual32=rpy["rows"]["actual32"],
            reusable="Frozen-state residuals and signed displacement; mean RPY correction [1.231,-1.787,-1.245]nm versus central PBD near zero.",
            scope="F7 no-flow initial-state free-RPY prediction, not this Y trajectory. Cannot multiply one contact correction by a step count or infer Y wall causation.",
            visible_bottleneck="CPU dense 96x96 mobility and 96x496 columns; dim=1 GPU kernel performs 16 ordered passes over 496 constraints with 96-coordinate updates. Serial update dependency requires validated scheduling, not an untested dim increase.",
            F10_priority="First bound geometry/contact splitting sensitivity on these recorded near-wall states and implement complete-state restart. Only then reconsider a separate small parallel contact prototype; no full RPY navigation refactor is justified by F9.",
        ),
    )


if __name__ == "__main__":
    main()
