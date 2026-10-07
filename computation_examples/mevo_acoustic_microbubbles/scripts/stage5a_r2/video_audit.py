"""Decode every MP4 frame and match displayed raw coordinates to source HDF5."""

import json
import subprocess
import h5py
import numpy as np
import imageio_ffmpeg
from scripts.stage5a_r2.seal import BASE
from scripts.stage5a_r2.validation_figures import local
from scripts.stage4l.dense_probe import sha, write


def main():
    root = BASE / "figures"
    video = root / "R7_OLD_NEW_ACTUAL_MOTION_60fps.mp4"
    source = json.loads((root / "R7_DATA_SOURCES.json").read_text())
    frames = dict(np.load(root / "R7_VIDEO_ACTUAL_FRAMES.npz"))
    clock = frames["requested_time_s"]
    raw = subprocess.check_output(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-v",
            "error",
            "-i",
            str(video),
            "-map",
            "0:v:0",
            "-f",
            "framemd5",
            "-",
        ],
        text=True,
    )
    (root / "R7_DECODED_FRAME_MD5.txt").write_text(raw)
    records = [s.split(",") for s in raw.splitlines() if s and not s.startswith("#")]
    timebase = next(s for s in raw.splitlines() if s.startswith("#tb 0:")).split(":", 1)[1].strip()
    pts = np.array([int(r[2]) for r in records])
    duration = np.array([int(r[3]) for r in records])
    checks = dict(
        decoded_frames=len(records),
        expected_frames=len(clock),
        timebase=timebase,
        constant_60_fps=bool(timebase == "1/60" and np.all(np.diff(pts) == 1) and np.all(duration == 1)),
        uniform_physical_clock=bool(np.allclose(np.diff(clock), np.diff(clock)[0], rtol=0, atol=1e-15)),
        physical_duration_s=float(clock[-1] - clock[0]),
        video_sha256=sha(video),
        render_source_snapshot_matches=sha(local(source["video"]["render_source_snapshot"]))
        == source["video"]["render_script_sha256"],
        video_hash_matches_provenance=sha(video) == source["video"]["sha256"],
        cases={},
    )
    for case, meta in source["cases"].items():
        for kind in ["reference", "candidate"]:
            key = f"{case}_{kind}"
            actual = frames[key + "_actual_time_s"]
            members = frames[key + "_particle_ids"]
            x = frames[key + "_positions_m"]
            with h5py.File(local(meta["report"][kind])) as h:
                ids = h["particle_id"][:]
                times = h["time_s"][:]
                index = [int(np.flatnonzero(ids == i)[0]) for i in members]
                exact = True
                for j, t in enumerate(actual):
                    k = int(np.searchsorted(times, t))
                    exact &= bool(times[k] == t and np.array_equal(x[j], h["positions_m"][k][index]))
            unique = len(np.unique(x.reshape(len(x), -1), axis=0))
            age = clock - actual
            checks["cases"][key] = dict(
                all_positions_exact_source_frames=exact,
                unique_actual_states=unique,
                actual_time_monotone=bool(np.all(np.diff(actual) >= 0)),
                max_sample_age_s=float(age.max()),
                future_frame_beyond_roundoff=bool(np.min(age) < -16 * np.finfo(float).eps * max(clock)),
                max_displacement_from_initial_m=float(np.linalg.norm(x - x[0], axis=2).max()),
            )
    checks["status"] = (
        "PASS"
        if checks["constant_60_fps"]
        and checks["render_source_snapshot_matches"]
        and checks["video_hash_matches_provenance"]
        and checks["uniform_physical_clock"]
        and len(records) == len(clock)
        and all(
            q["all_positions_exact_source_frames"]
            and q["unique_actual_states"] > 1
            and q["actual_time_monotone"]
            and not q["future_frame_beyond_roundoff"]
            for q in checks["cases"].values()
        )
        else "FAIL"
    )
    checks["scope"] = (
        "Every decoded video frame has a uniform timestamp. Every saved rendered coordinate is bitwise equal to an actual accepted source state; no coordinate interpolation. Raw frame archive and renderer source link pixels to these coordinates."
    )
    write(root / "R7_VIDEO_VALIDATION.json", checks)
    print(json.dumps(checks, indent=2))
    assert checks["status"] == "PASS"


if __name__ == "__main__":
    main()
