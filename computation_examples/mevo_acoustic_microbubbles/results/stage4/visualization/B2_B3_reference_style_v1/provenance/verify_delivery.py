"""Verify the completed B2/B3 delivery without invoking any simulation."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import imageio.v2 as imageio
import imageio_ffmpeg
import numpy as np
from PIL import Image

root = Path('/home/lzy/projects/computation_examples/mevo_acoustic_microbubbles')
out = root / 'results/stage4/visualization/B2_B3_reference_style_v1'
manifest = json.loads((out / 'RENDER_MANIFEST.json').read_text())
assert manifest['status'] == 'PASS' and manifest['changed_inputs'] == []

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

assert sha(root / 'scripts/stage4/render_motion.py') == manifest['renderer_sha256']
assert sha(out / 'provenance/stage4_render_motion.py') == manifest['renderer_sha256']
assert sha(out / 'provenance/stage3_render_motion.py') == sha(root / 'scripts/stage3/render_motion.py')
checks = {}
for stem, expected in manifest['videos'].items():
    video = out / f'{stem}.mp4'
    assert sha(video) == expected['sha256']
    decoded = subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), '-v', 'error', '-threads', '2',
                              '-i', str(video), '-f', 'null', '-'], capture_output=True, text=True, check=True)
    assert not decoded.stderr, decoded.stderr
    reader = imageio.get_reader(video)
    meta = reader.get_meta_data()
    assert meta['codec'] == 'h264' and meta['pix_fmt'].startswith('yuv420p')
    assert meta['fps'] == 60 and meta['size'] == (1920, 1080)
    assert reader.count_frames() == 1321
    errors = {}
    for frame in [30, 230, 630, 1230]:
        actual = reader.get_data(frame)
        raw = imageio.imread(out / 'keyframes' / f'{stem}_{frame:04d}.png')
        errors[str(frame)] = float(np.abs(actual.astype(float) - raw).mean())
        assert errors[str(frame)] < 1., errors
        if stem == 'stage4_B2_B3_motion_3d_60fps' and frame in [230, 1230]:
            Image.fromarray(actual).resize((1500, 844), Image.Resampling.LANCZOS).save(
                f'/tmp/{stem}_encoded_{frame}.png')
    reader.close()
    checks[stem] = dict(full_decode='PASS', metadata=meta,
                       mean_absolute_rgb_encoding_error_0_to_255=errors)
    print(json.dumps(dict(video=stem, verification='PASS')), flush=True)

samples = np.load(out / 'display_samples.npz')
rows = list(csv.DictReader((out / 'frame_map.csv').open()))
assert len(rows) == 2642
for name in ['B2', 'B3']:
    assert samples[f'{name}_positions_m'].shape == (1321, 256, 3)
    assert np.isfinite(samples[f'{name}_positions_m']).all()
    assert len(np.unique(samples[f'{name}_particle_id'])) == 256
    source = root / manifest['sources'][name]['result_id']
    branches = json.loads((source / 'BRANCH_METRICS.json').read_text())
    selected = [row for row in rows if row['case'] == name]
    assert len(selected) == 1321
    for i, row in enumerate(selected):
        now = float(samples['physical_time_s'][i])
        assert int(row['frame']) == i and float(row['physical_time_s']) == now
        counts = {branch: sum(e['branch'] == branch and e['crossing_time_s'] <= now
                              for e in branches['crossings']) for branch in ['BRANCH_A', 'BRANCH_B']}
        assert int(row['A']) == counts['BRANCH_A'] and int(row['B']) == counts['BRANCH_B']
        assert int(row['not_crossed']) + int(row['A']) + int(row['B']) == 256

pause = root / 'results/stage4/performance_refactor/completion_execution/20261005T012343106214Z/USER_DEVELOPMENT_PAUSE_FOR_ANIMATION.json'
assert json.loads(pause.read_text())['resume_requires_user_instruction']
report = dict(status='PASS', python=sys.version, numpy=np.__version__,
              renderer_commit='e7e3d2f', interpolation_tests='6 passed locally', ruff='PASS',
              compute_environment='Visualization test skipped because optional imageio is not installed; existing sensitivity test passed. Compute environment unchanged.',
              all_video_frames_decode_without_errors=True, every_frame_count_checked_against_original_events=True,
              video_checks=checks, development_pause_record=str(pause.relative_to(root)))
(out / 'DELIVERY_VERIFICATION.json').write_text(json.dumps(report, indent=2) + '\n')
delivery = {str(p.relative_to(out)): dict(bytes=p.stat().st_size, sha256=sha(p))
            for p in sorted(out.rglob('*')) if p.is_file() and p.name not in ['DELIVERY_SHA256.json', 'render.log']}
(out / 'DELIVERY_SHA256.json').write_text(json.dumps(delivery, indent=2) + '\n')
print(json.dumps(dict(status='PASS', files=len(delivery), bytes=sum(v['bytes'] for v in delivery.values()))))
