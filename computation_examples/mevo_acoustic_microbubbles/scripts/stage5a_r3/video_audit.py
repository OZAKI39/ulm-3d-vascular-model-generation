"""Decode every rendered frame using the installed imageio FFmpeg binary."""

import subprocess
import imageio.v2 as imageio
import imageio_ffmpeg
import numpy as np
from scripts.stage5a_r3.common import BASE, read, write, sha


def audit():
    folder = BASE / "videos"
    path = folder / "B2_actual_wall_stall_crossing_60fps.mp4"
    mapping = read(folder / "B2_actual_wall_stall_crossing_FRAME_MAP.json")
    decoded = subprocess.check_output(
        [
            imageio_ffmpeg.get_ffmpeg_exe(),
            "-v",
            "error",
            "-i",
            str(path),
            "-map",
            "0:v:0",
            "-f",
            "framemd5",
            "-",
        ],
        text=True,
    )
    (folder / "B2_actual_wall_stall_crossing_DECODED_FRAME_MD5.txt").write_text(decoded)
    records = [line.split(",") for line in decoded.splitlines() if line and not line.startswith("#")]
    timebase = (
        next(line for line in decoded.splitlines() if line.startswith("#tb 0:")).split(":", 1)[1].strip()
    )
    pts = np.array([int(r[2]) for r in records])
    duration = np.array([int(r[3]) for r in records])
    constant = timebase == "1/60" and bool(np.all(np.diff(pts) == 1) and np.all(duration == 1))
    reader = imageio.get_reader(path)
    metadata = reader.get_meta_data()
    for fraction in [0, 0.5, 1.0]:
        index = int(fraction * (len(records) - 1))
        imageio.imwrite(folder / f"B2_video_frame_{index:04d}.png", reader.get_data(index))
    reader.close()
    result = dict(
        status="PASS" if constant and len(records) == len(mapping["frames"]) else "FAIL",
        frames_decoded=len(records),
        frames_expected=len(mapping["frames"]),
        constant_60fps=constant,
        timebase=timebase,
        dimensions=metadata["size"],
        video_sha256=sha(path),
        frame_map_sha256=sha(folder / "B2_actual_wall_stall_crossing_FRAME_MAP.json"),
        rendered_source_sha256=sha(folder / "B2_RENDER_SOURCE.py")
        if (folder / "B2_RENDER_SOURCE.py").exists()
        else None,
        decoder=imageio_ffmpeg.get_ffmpeg_exe(),
        no_interpolation=mapping["trajectory_interpolation"] is False,
        note="Uniform decoded timestamps; repeated frames correspond to real stored-state holds, not dropped video frames",
    )
    write(folder / "B2_actual_wall_stall_crossing_VIDEO_AUDIT.json", result)
    print(result, flush=True)
    assert result["status"] == "PASS"


if __name__ == "__main__":
    audit()
