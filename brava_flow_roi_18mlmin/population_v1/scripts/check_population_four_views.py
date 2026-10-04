"""Validate the four clips and independently check sampled particle speeds."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import subprocess
import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'candidate_0/views_0to1p2s'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def main():
    global OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=OUT)
    parser.add_argument('--compare-with', type=Path)
    args = parser.parse_args()
    OUT = args.output_dir.resolve()
    views = json.loads((OUT/'data/views.json').read_text())
    assert len(views) == 4
    baseline = None
    protections = {}
    checks = []
    for view in views:
        target = OUT/view['name']
        manifest = json.loads((target/'data/render_manifest.json').read_text())
        protections.update(json.loads((target/'data/protected_files_sha256.json').read_text()))
        with np.load(target/'data/animation_samples.npz') as file:
            current = {key: file[key] for key in file.files}
        if baseline is None:
            baseline = current
        else:
            assert all(np.array_equal(value, baseline[key], equal_nan=True) for key, value in current.items())
        t = current['frame_time_s']; shown = current['displayed']
        assert abs(t[0]) < 1e-14 and abs(t[-1]-1.2) < 1e-14
        assert np.allclose(np.diff(t), t[-1]/(len(t)-1), rtol=0, atol=1e-14)
        assert np.array_equal(shown, np.isfinite(current['position_m'][:, :, 0]) &
                              ~current['prior_exclusion'][None, :] & (current['mean_speed_m_s'] >= .002))
        layout = manifest['layout']
        annotations = json.loads((target/'data/coordinate_annotations.json').read_text())['frames'][0]['labels']
        visible_labels = [r for r in annotations if r['opacity'] >= .5]
        for i, label in enumerate(visible_labels):
            a = label['bbox_px']
            assert a[0] > 26 and a[2] < 1894 and a[1] > 204 and a[3] < 1054
            for other in visible_labels[i+1:] + layout['port_labels']:
                b = other['bbox_px']
                assert not (a[0] < b[2] and a[2] > b[0] and a[1] < b[3] and a[3] > b[1]), (label, other)
        if args.compare_with:
            old_target = args.compare_with/view['name']
            old_manifest = json.loads((old_target/'data/render_manifest.json').read_text())
            assert sha(old_target/old_manifest['movie_path']) == old_manifest['movie_sha256']
            with np.load(old_target/'data/animation_samples.npz') as old_samples:
                assert all(np.array_equal(current[k], old_samples[k], equal_nan=True) for k in current)
            for key in ['camera_position', 'camera_focal_point', 'parallel_scale', 'angle_labels']:
                assert layout[key] == old_manifest['layout'][key], key
            for key in ['background','bubble_color','marker_size_px','vessel_opacity','frames','fps','physical_time_s']:
                assert manifest[key] == old_manifest[key], key
            typography = layout['coordinate_typography']
            assert typography['axis_title_size_px'] == 40 and typography['tick_size_px'] == 30
        direction = np.array(layout['camera_position'])-layout['camera_focal_point']
        direction /= np.linalg.norm(direction)
        actual_azimuth = np.degrees(np.arctan2(direction[1], direction[0])) % 360
        actual_elevation = np.degrees(np.arcsin(direction[2]))
        az_error = abs((actual_azimuth-layout['world_azimuth_deg']+180) % 360-180)
        el_error = abs(actual_elevation-layout['elevation_deg'])
        assert az_error < 1e-10 and el_error < 1e-10
        assert layout['angle_labels'][0]['text'] == f"Azimuth offset: +{view['azimuth_offset_deg']:.1f}°"
        assert layout['angle_labels'][1]['text'] == f"Elevation: {view['elevation_deg']:.1f}°"
        assert all(b['bbox_px'][0] > 700 and b['bbox_px'][2] < 1920 and b['bbox_px'][3] < 200 for b in layout['angle_labels'])
        bbox = np.array(layout['vessel_bounds_px']); safe = np.array(layout['safe_vessel_bounds_px'])
        assert np.all(bbox[0] >= safe[0]-1e-7) and np.all(bbox[1] <= safe[1]+1e-7)
        for port in layout['port_labels']:
            box = port['bbox_px']
            assert box[0] >= 25 and box[2] <= 1895 and box[1] >= 207 and box[3] <= 1055
            assert 0 < port['cap_px'][0] < 1920 and 0 < port['cap_px'][1] < 1080
        movie = OUT/view['video']
        probe = json.loads(subprocess.check_output(['ffprobe', '-v', 'error', '-count_frames', '-select_streams', 'v:0',
            '-show_entries', 'stream=codec_name,width,height,r_frame_rate,nb_read_frames,duration', '-of', 'json', str(movie)]))['streams'][0]
        assert (probe['width'], probe['height']) == (1920, 1080) and probe['r_frame_rate'] == '24/1'
        assert int(probe['nb_read_frames']) == manifest['frames']
        assert abs(float(probe['duration'])-manifest['video_duration_s']) < 1e-5
        result = subprocess.run(['ffmpeg', '-v', 'error', '-i', str(movie), '-f', 'null', '-'], capture_output=True, text=True)
        assert result.returncode == 0 and not result.stderr.strip(), result.stderr
        decoded = target/'figures/video_end.png'
        subprocess.run(['ffmpeg', '-v', 'error', '-y', '-i', str(movie), '-vf', f"select=eq(n\\,{manifest['frames']-1})",
                        '-frames:v', '1', str(decoded)], check=True)
        a = np.array(Image.open(target/'figures/snapshot_end.png')).astype(float)
        b = np.array(Image.open(decoded)).astype(float)
        header_error = float(abs(a[:200, 700:]-b[:200, 700:]).mean())
        assert header_error < 2., header_error
        checks.append(dict(name=view['name'], status='PASS', full_video_decode=True, probe=probe,
                           movie_sha256=sha(movie), azimuth_error_deg=az_error, elevation_error_deg=el_error,
                           angle_label_encode_mean_abs_error_8bit=header_error,
                           all_views_same_particles_and_sampling=True,
                           visible_coordinate_label_count=len(visible_labels),
                           coordinate_labels_no_overlap_or_clipping=True,
                           coordinate_typography=layout.get('coordinate_typography'),
                           sampled_data_and_camera_unchanged_from_previous=bool(args.compare_with),
                           final_displayed_count=int(shown[-1].sum()), maximum_displayed_count=int(shown.sum(1).max())))
        print('CHECKED', view['name'], flush=True)

    assert all(sha(path) == digest for path, digest in protections.items())
    track_hashes = json.loads((OUT/views[0]['name']/'data/track_hashes.json').read_text())
    source = ROOT/'analysis/realization_0000_contactfix'
    check_frames = sorted({0, len(t)//4, len(t)//2, 3*len(t)//4, len(t)-1})
    max_speed_error = max_position_error = 0.
    comparisons = 0
    for i, pid in enumerate(baseline['particle_id']):
        path = source/'tracks'/f'mb_{pid:06d}.npz'
        assert sha(path) == track_hashes[path.name]
        with np.load(path) as file:
            ts, xs = file['physical_time_s'], file['x_m']
            segment_speed = np.linalg.norm(np.diff(xs, axis=0), axis=1)/np.diff(ts)
            for j in check_frames:
                if not ts[0] <= t[j] < ts[-1]:
                    continue
                weight = np.maximum(0., np.minimum(ts[1:], t[j]+.05)-np.maximum(ts[:-1], t[j]-.05))
                avg = np.dot(segment_speed, weight)/weight.sum()
                max_speed_error = max(max_speed_error, abs(avg-baseline['mean_speed_m_s'][j, i]))
                x = np.array([np.interp(t[j], ts, xs[:, k]) for k in range(3)])
                max_position_error = max(max_position_error, float(np.linalg.norm(x-baseline['position_m'][j, i])))
                comparisons += 1
    assert max_speed_error < 1e-10 and max_position_error < 1e-12
    old_layout = json.loads((ROOT/'candidate_0/view05_flow_axes/data/layout.json').read_text())
    new_layout = json.loads((OUT/views[0]['name']/'data/layout.json').read_text())
    current_camera_error = max(float(np.max(abs(np.array(old_layout[k])-new_layout[k])))
                               for k in ['camera_position', 'camera_focal_point', 'parallel_scale'])
    assert current_camera_error < 1e-10
    report = dict(status='PASS', video_count=4, views=checks,
        all_four_views_identical_sampled_positions_speeds_and_visibility=True,
        source_track_files_unchanged=len(track_hashes), protected_files_unchanged=len(protections),
        independent_window_checks=comparisons, physical_check_times_s=t[check_frames].tolist(),
        maximum_speed_vs_recorded_displacement_error_m_s=max_speed_error,
        maximum_position_vs_source_error_m=max_position_error,
        current_view_camera_matches_previous_maximum_error=current_camera_error,
        low_speed_filter='Rolling t +/- 0.05 s time-weighted 3-D speed >= 2 mm/s; use recorded lifetime overlap',
        simulation_started=False)
    dump(OUT/'data/QC.json', report)
    dump(OUT/'data/protected_files_sha256.json', protections)
    with (OUT/'data/view_summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['view','azimuth_offset_deg','world_azimuth_deg','elevation_deg','physical_start_s','physical_end_s','fps','frames','video_duration_s','final_displayed_count','video'])
        for view, check in zip(views, checks):
            writer.writerow([view['name'],view['azimuth_offset_deg'],view['world_azimuth_deg'],view['elevation_deg'],0,1.2,view['fps'],view['frames'],view['video_duration_s'],check['final_displayed_count'],view['video']])
    contact = Image.new('RGB', (1920, 1080), 'black')
    for i, view in enumerate(views):
        image = Image.open(OUT/view['name']/'figures/video_end.png').convert('RGB').resize((960, 540))
        contact.paste(image, (i % 2*960, i//2*540))
    contact.save(OUT/'CONTACT_SHEET.png')
    print('ALL_FOUR_PASS', json.dumps({k: report[k] for k in ['independent_window_checks','maximum_speed_vs_recorded_displacement_error_m_s','current_view_camera_matches_previous_maximum_error']}), flush=True)


if __name__ == '__main__':
    main()
