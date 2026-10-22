"""Validate the current presentation only; never launch physics calculations."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess

import numpy as np
from PIL import Image

root = Path(__file__).resolve().parents[1]
out = root / 'candidate_0'
backup = out / 'revisions/before_black_labels_maxfit_20260929'
record = json.loads((out / 'data/render_manifest.json').read_text())
prior = json.loads((backup / 'data/render_manifest.json').read_text())
assert not record['preview_only'] and record['candidate_id'] == 0
assert record['background'] == '#000000'
assert record['uniform_particle_color'] == '#3399ff'
for key in ['flow_sha256', 'source_identity_sha256', 'production_git_commit',
            'frames', 'fps', 'physical_time_end_s', 'marker_size_display_px',
            'nominal_dt_s', 'inlet_flow_m3_s']:
    assert record[key] == prior[key], key

old_direction = (np.array(prior['layout']['camera_focal_point'])
                 - np.array(prior['layout']['camera_position']))
old_direction /= np.linalg.norm(old_direction)
assert np.allclose(record['layout']['view_direction'], old_direction, atol=1e-14)
bounds = np.array(record['layout']['content_bounds_px'])
safe = np.array(record['layout']['safe_bounds_px'])
assert np.all(bounds[0] >= safe[0] - 1e-4)
assert np.all(bounds[1] <= safe[1] + 1e-4)
assert np.isclose(np.max((bounds[1]-bounds[0])/(safe[1]-safe[0])), 1., atol=1e-12)

catalog = list(csv.DictReader((root / 'analysis/realization_0000_contactfix/trajectory_catalog.csv').open()))
overlap = list(csv.DictReader((root / 'cfd_wall_audit/acceleration_overlap.csv').open()))
failures = {int(r['particle_id']) for r in catalog if r['status'] == 'SOLVER_FAILURE'}
risks = {int(r['particle_id']) for r in overlap if r['window_any_eta_gt_0p5'] == 'True'}
excluded = failures | risks
ids = {int(r['particle_id']) for r in catalog}
assert set(record['hidden_particle_ids']) == excluded
assert set(record['displayed_particle_ids']) == ids - excluded
selection = list(csv.DictReader((out / 'data/display_selection.csv').open()))
assert {int(r['particle_id']) for r in selection if r['displayed'] == 'False'} == excluded
assert len(selection) == len(ids)

protected = json.loads((backup / 'protected_files_sha256.json').read_text())
for name, expected in protected.items():
    with (root / name).open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == expected, name
for name, expected in record['display_filter']['sources_sha256'].items():
    with (root / name).open('rb') as handle:
        assert hashlib.file_digest(handle, 'sha256').hexdigest() == expected, name

samples = np.load(out / 'data/animation_samples.npz')
old_samples = np.load(backup / 'data/animation_samples.npz')
for key in samples.files:
    assert np.array_equal(samples[key], old_samples[key], equal_nan=True), key

geometry = json.loads((root.parent / 'inputs/geometry.json').read_text())
transform = np.array(geometry['pose_transforms']['0']['source_mm_to_candidate_mm'])
grid = record['coordinate_grid']
assert grid['type'] == 'SPATIAL_RIGHT_PANEL_3D_SIX_FACE_GRID'
assert grid['frame_color'] == '#FF5C5C' and grid['grid_opacity'] == .28
for i,(axis,count) in enumerate([('X',4),('Y',4),('Z',5)]):
    expected = np.linspace(grid['bounds_candidate_mm'][2*i], grid['bounds_candidate_mm'][2*i+1], count)
    assert np.array_equal(grid['tick_values_mm'][axis], expected)
    labels = [r['text'] for r in grid['annotations'] if r['axis']==axis]
    assert set(labels)=={f'{v:.2f}' for v in expected}|{f'{axis} (mm)'}
origin = np.array(record['layout']['camera_focal_point'])
right = np.cross(old_direction, [0.,0.,1.]);right /= np.linalg.norm(right)
up = np.cross(right, old_direction)
pixels_per_mm = 1080 / (2*record['layout']['parallel_scale'])
assert np.isclose(np.linalg.norm(right), 1) and np.isclose(np.linalg.norm(up), 1)
assert abs(np.dot(right, up)) < 1e-14
projection_error = 0.
for role, center in record['port_centers_m'].items():
    world = (np.array(center)*1000 + geometry['origin_source_mm']) @ transform[:3,:3].T + transform[:3,3]
    pixel = [960 + np.dot(world-origin, right)*pixels_per_mm,
             540 - np.dot(world-origin, up)*pixels_per_mm]
    projection_error = max(projection_error, float(np.max(abs(np.array(pixel)-record['ports_px'][role]))))
    assert projection_error < 1e-8
    annotation = record['port_annotations'][role]
    assert annotation['text'] == ('Inlet' if role == 'INLET' else 'Outlet '+str(int(role[-2:])))
    assert annotation['cap_px'] == record['ports_px'][role]
    bb = annotation['text_bbox_px']
    assert 20 <= bb[0] < bb[2] <= 1900 and 200 <= bb[1] < bb[3] <= 1056

video = out / 'animations/microbubble_physical_time.mp4'
probe = json.loads(subprocess.check_output([
    'ffprobe', '-v', 'error', '-select_streams', 'v:0', '-show_entries',
    'stream=width,height,pix_fmt,avg_frame_rate,nb_frames,duration', '-of', 'json', str(video)
]))['streams'][0]
assert (probe['width'], probe['height'], probe['avg_frame_rate']) == (1920, 1080, '24/1')
assert int(probe['nb_frames']) == 2880 and float(probe['duration']) == 120.
subprocess.run(['ffmpeg', '-v', 'error', '-xerror', '-threads', '2', '-i', str(video), '-f', 'null', '-'], check=True)
for name, time in [('start', 0), ('middle', 60), ('late', 115)]:
    target = out / 'figures' / f'qc_video_{name}.png'
    subprocess.run(['ffmpeg', '-v', 'error', '-y', '-ss', str(time), '-i', str(video),
                    '-frames:v', '1', str(target)], check=True)
    pixels = np.array(Image.open(target).convert('RGB'))
    assert int(pixels[0:15, 0:15].max()) == 0

qc = dict(status='PASS', candidate_id=0, video=probe,
          video_sha256=hashlib.file_digest(video.open('rb'), 'sha256').hexdigest(),
          full_video_decode_pass=True, projected_content_bounds_px=bounds.tolist(),
          maximum_safe_fit=True, viewing_direction_unchanged=True,
          source_particle_count=len(ids), displayed_particle_count=len(ids-excluded),
          hidden_particle_count=len(excluded), hidden_failure_count=len(failures),
          hidden_numerical_risk_count=len(risks), exclusion_ids_match_existing_audit=True,
          all_animation_sample_arrays_unchanged=True, protected_file_count=len(protected),
          candidate_15_unchanged=True, coordinate_grid='PHYSICAL_CANDIDATE_XYZ_MM_REFERENCE_RIGHT_PANEL',
          port_coordinate_projection_max_error_px=projection_error,
          corner_background_black=True)
(out / 'data/render_qc.json').write_text(json.dumps(qc, indent=2)+'\n')
print(json.dumps(qc, indent=2))
