"""Read-only numerical/provenance/media checks for the view-05 animation."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'candidate_0/view05_flow_axes'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    manifest = json.loads((OUT/'data/render_manifest.json').read_text())
    protected = json.loads((OUT/'data/protected_files_sha256.json').read_text())
    assert all(sha(p) == digest for p, digest in protected.items())
    tracks = json.loads((OUT/'data/track_hashes.json').read_text())
    track_dir = ROOT/'analysis/realization_0000_contactfix/tracks'
    assert all(sha(track_dir/name) == digest for name, digest in tracks.items())
    with np.load(OUT/'data/animation_samples.npz') as samples:
        time = samples['frame_time_s']; positions = samples['position_m']
        speed = samples['mean_speed_m_s']; ids = samples['particle_id']
        displayed = samples['displayed']; prior = samples['prior_exclusion']
        assert np.allclose(np.diff(time), time[-1]/(len(time)-1), rtol=0, atol=1e-14)
        assert np.array_equal(displayed, np.isfinite(positions[:, :, 0]) & ~prior[None, :] & (speed >= .002))
        source = json.loads((ROOT/'candidate_0/static_t6s_white_frame/data/manifest.json').read_text())
        assert np.array_equal(ids[displayed[-1]], source['displayed_particle_ids'])
        mean_error = 0.
        position_error = 0.
        reference_rows = list(csv.DictReader((ROOT/'candidate_0/static_t6s_white_frame/data/display_selection.csv').open()))
        by_id = {int(r['particle_id']): r for r in reference_rows}
        for i, pid in enumerate(ids):
            row = by_id[int(pid)]
            if row['active_at_6s'] == 'True':
                mean_error = max(mean_error, abs(speed[-1, i]-float(row['mean_speed_mm_s'])/1000))
                position_error = max(position_error, float(np.linalg.norm(positions[-1, i]-[float(row[k]) for k in ['x_m', 'y_m', 'z_m']])))
        assert mean_error < 1e-11, mean_error
        assert position_error < 1e-12, position_error
    movie = OUT/'animations/microbubble_view05_physical_time.mp4'
    probe = json.loads(subprocess.check_output(['ffprobe', '-v', 'error', '-count_frames', '-select_streams', 'v:0',
        '-show_entries', 'stream=codec_name,width,height,r_frame_rate,nb_read_frames,duration', '-of', 'json', str(movie)]))['streams'][0]
    assert (probe['width'], probe['height']) == (1920, 1080)
    assert probe['r_frame_rate'] == '24/1'
    assert int(probe['nb_read_frames']) == manifest['frames']
    assert abs(float(probe['duration'])-manifest['video_duration_s']) < 1e-5
    result = subprocess.run(['ffmpeg', '-v', 'error', '-i', str(movie), '-f', 'null', '-'], capture_output=True, text=True)
    assert result.returncode == 0 and not result.stderr.strip(), result.stderr
    for name, index in [('start', 0), ('middle', manifest['frames']//2), ('end', manifest['frames']-1)]:
        subprocess.run(['ffmpeg', '-v', 'error', '-y', '-i', str(movie), '-vf', f'select=eq(n\\,{index})',
            '-frames:v', '1', str(OUT/'figures'/f'video_{name}.png')], check=True)
    report = dict(status='PASS', full_video_decode=True, probe=probe, movie_sha256=sha(movie),
        source_track_files_unchanged=len(tracks), protected_files_unchanged=len(protected),
        time_sampling_uniform=True, visibility_matches_filters_every_frame=True,
        terminal_particle_ids_match_previous_t6_snapshot=True,
        maximum_t6_speed_vs_previous_direct_window_average_error_m_s=mean_error,
        maximum_t6_position_vs_previous_snapshot_error_m=position_error,
        camera_direction_error=manifest['layout']['direction_error'],
        zoom_relative_to_reference=manifest['layout']['zoom_relative_to_reference'])
    (OUT/'data/QC.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
