"""Independent data/identity checks and full MP4 decode before final gates."""

import argparse
import json
import re
import subprocess
from pathlib import Path
import h5py
import imageio_ffmpeg
import numpy as np
from mevo_mb.flow.fields import SimVascularFlowField
from .common import BASE, INPUT, ROOT, js, write, sha, F4_ARCHIVE


def main(label, checks):
    run = BASE / label
    visual = BASE / "visualization"
    check = js(BASE / checks / "RESULT.json")
    assert check["status"] == "PASS"
    source = js(BASE / "FLOW_SOURCE_IDENTITY.json")
    assert sha(INPUT / "result_121.vtu") == "131f7cb9b3e4ac9887bb0aa344d0eb2c0a2c8ac1d3b19fafb093ecd9a1045620"
    assert (
        sha(INPUT / "f4/labelled_lumen_m.vtp")
        == "e14b06f1164be3beab261ec052dbe7edad6c2394d80077787cb318cbece15f4b"
    )
    cfg = js(run / "main/CONFIG.json")
    a = np.load(run / "branch_A_actual_suffix.npz")
    with h5py.File(run / "main/trajectory.h5") as h, h5py.File(run / "branch_B/trajectory.h5") as b:
        x = h["positions_m"][:, 0]
        t = h["time_s"][:]
        np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
        np.testing.assert_array_equal(h["radii_m"][:], np.full(32, 1.78e-6))
        np.testing.assert_allclose(t, np.arange(1201) * 0.0005, atol=1e-15, rtol=0)
        np.testing.assert_array_equal(x[0], np.array(cfg["positions_m"]))
        np.testing.assert_array_equal(b["positions_m"][0, 0], x[600])
        np.testing.assert_array_equal(a["positions_m"], h["positions_m"][600:])
        assert h["active"][:].all() and b["active"][:].all()
        np.testing.assert_array_equal(h["flow/acoustic_enabled"][:], np.arange(1201) < 600)
        assert h["flow/flow_enabled"][:].all() and not b["flow/flow_enabled"][:].any()
        assert not b["flow/acoustic_enabled"][:].any()
        assert h["flow/valid"][:].all() and np.isfinite(x).all()
        cpu = SimVascularFlowField(INPUT / "result_121.vtu", spatial_method="tetra_p1")
        ids = np.arange(0, 1201, 60)
        ref = cpu.velocity(x[ids].reshape(-1, 3), 0).reshape(len(ids), 32, 3)
        recorded = h["flow/velocity_m_s"][ids, 0]
        error = float(np.linalg.norm(ref - recorded, axis=-1).max())
        assert error < 2e-12
        max_jump = float(np.linalg.norm(np.diff(x, axis=0), axis=-1).max())
        q = h["quality_samples/raw_F1_diagonal"][:]
        wall = h["wall_stats"][:]
        contact = dict(
            max_accepted_pair_overlap_m=float(q[..., 3].max()),
            max_pair_pass_correction_m=float(q[..., 6].max()),
            minimum_saved_wall_gap_m=float(h["wall_gap_m"][:].min()),
            max_wall_pass_correction_m=float(wall[..., 2].max()),
            pair_correction_ON_sum_m=float(q[600, ..., 4].sum()),
            pair_correction_OFF_sum_m=float(q[-1, ..., 4].sum() - q[600, ..., 4].sum()),
            wall_correction_ON_sum_m=float(wall[600, ..., 3].sum()),
            wall_correction_OFF_sum_m=float(wall[-1, ..., 3].sum() - wall[600, ..., 3].sum()),
            B_correction_over_free_path=None,
            B_ratio_reason="No free displacement; inherited raw ratio uses a 1e-300 denominator and is undefined physically. Use absolute contact relaxation, not the raw ratio.",
        )
        assert contact["max_accepted_pair_overlap_m"] < 0.01 * 1.78e-6
        assert contact["minimum_saved_wall_gap_m"] >= -0.01 * 1.78e-6
        assert contact["max_wall_pass_correction_m"] < 0.1 * 1.78e-6
    video = visual / "Video_F5_acoustic_release_flow_60fps.mp4"
    completed = subprocess.run(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-v",
            "info",
            "-i",
            str(video),
            "-map",
            "0:v:0",
            "-vf",
            "showinfo",
            "-f",
            "null",
            "-",
        ],
        capture_output=True,
        text=True,
        check=True,
    )
    log = completed.stderr
    (visual / "FULL_VIDEO_DECODE.log").write_text(log)
    frames = re.findall(r"n:\s*(\d+)\s+pts:\s*(\d+)\s+pts_time:([-+\d.e]+)", log)
    assert len(frames) == 1201, len(frames)
    ids = np.array([int(r[0]) for r in frames])
    pts = np.array([int(r[1]) for r in frames])
    np.testing.assert_array_equal(ids, np.arange(1201))
    assert (np.diff(pts) == np.diff(pts)[0]).all()
    unit = re.search(r"config in time_base:\s*(\d+)/(\d+),\s*frame_rate:\s*(\d+)/(\d+)", log)
    assert unit, "Decoder time base and frame rate must be recorded"
    num, den, frn, frd = map(int, unit.groups())
    assert frn / frd == 60
    np.testing.assert_allclose(pts * num / den, np.arange(1201) / 60, atol=1e-12, rtol=0)
    mapping = np.load(visual / "FRAME_MAP.npz")
    np.testing.assert_array_equal(mapping["source_frame"], np.arange(1201))
    np.testing.assert_array_equal(mapping["physical_time_s"], t)
    np.savez_compressed(visual / "DECODED_FRAME_TIMESTAMPS.npz", frame=ids, pts=pts, time_base=[num, den])
    for i in [1, 2, 3]:
        assert len(list(visual.glob(f"Fig_F5_{i}_*.png"))) >= 1
        assert len(list(visual.glob(f"Fig_F5_{i}_*.pdf"))) == 1
    vis = dict(
        status="PASS",
        frames=1201,
        fps=60,
        full_decode=True,
        constant_presentation_timestamps=True,
        source_frames_one_to_one=True,
        physical_sample_interval_s=0.0005,
        video_duration_s=1201 / 60,
        displayed_physical_duration_s=0.6,
        switch_frame=600,
        switch_video_time_s=10.0,
        position_interpolation=False,
        source_sha256=sha(run / "main/trajectory.h5"),
        video_sha256=sha(video),
        scope="Encoding, timestamps and source mapping; no claim of local physical validation",
    )
    write(visual / "VISUALIZATION_VALIDATION.json", vis)
    changed = subprocess.check_output(
        [
            "git",
            "diff",
            "--name-only",
            F4_ARCHIVE,
            "--",
            "src",
            "scripts/fast_pzt",
            "scripts/fast_pzt_f2",
            "scripts/fast_pzt_f3",
            "scripts/fast_pzt_f4",
            "configs/fast_pzt/f4",
        ],
        cwd=ROOT,
        text=True,
    ).splitlines()
    allowed = {"flow_gpu.py", "f5_engine.py"}
    assert all(Path(p).name in allowed for p in changed), changed
    numerical = dict(
        status="PASS",
        sampled_actual_trajectory_flow_CPU_error_m_s=error,
        CPU_check_points=21 * 32,
        maximum_saved_frame_particle_displacement_m=max_jump,
        contact=contact,
        archived_F1_F4_sources_unchanged=True,
        no_particle_exit_or_field_failure=True,
        same_switch_state_bitwise=True,
        branch_A_exact_main_suffix=True,
        main_steps=60000,
        branch_B_steps=30000,
        dt_s=1e-5,
        adaptivity=False,
    )
    write(BASE / "DELIVERY_VALIDATION.json", numerical)
    gate = dict(
        status="PASS_WITH_SCOPE",
        FLOW_IMPORT="PASS",
        GPU_FLOW_INTERPOLATION="PASS",
        FLOW_COUPLING="PASS",
        ATTRACTION_RELEASE="PASS",
        CLUSTER_TRANSPORT="PASS",
        VESSEL_CONFINEMENT="PASS_WITH_NUMERICAL_TOLERANCE",
        PERFORMANCE="PASS",
        VISUALIZATION="PASS",
        PHYSICS_SCOPE=["EXPLORATORY", "PRESCRIBED_STEADY_FLOW", "GEOMETRIC_WALL_ONLY"],
        checks=checks,
        science=label,
        numerical_validation=numerical,
        video=vis,
        limitations=[
            "Stored water-parameter CFD, not validated blood rheology or mesh-independent near-wall CFD",
            "One-way prescribed flow; truncated local pair mobility and central PBD; no wall hydrodynamic lubrication",
            "Distance-based cluster split, not physical bubble fragmentation; last split is right-censored",
            "No numerical convergence or biological validation follows from GPU/CPU agreement or a smooth video",
        ],
        source_identity_file="FLOW_SOURCE_IDENTITY.json",
        flow_identity_retained=bool(source),
        outcome_file="ON_OFF_FLOW_COMPARISON.json",
        performance_file="PERFORMANCE.json",
    )
    write(BASE / "FINAL_GATE.json", gate)
    print(json.dumps(dict(status=gate["status"], video_frames=1201, CPU_flow_error_m_s=error)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="science_v1")
    ap.add_argument("--checks", default="checks_v3")
    args = ap.parse_args()
    main(args.label, args.checks)
