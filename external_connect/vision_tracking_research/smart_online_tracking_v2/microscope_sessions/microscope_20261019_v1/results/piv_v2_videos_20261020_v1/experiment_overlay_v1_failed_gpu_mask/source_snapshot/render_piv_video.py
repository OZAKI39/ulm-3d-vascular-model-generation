"""Bounded-memory offline PIV on source pixels, with existing candidate overlays.

No physical scale or verified channel mask is assumed. Velocities are apparent
image motion in px/s; candidate detections and instrument readback are read-only.
"""
from pathlib import Path
from fractions import Fraction
from collections import deque
import argparse
import csv
import gzip
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import time

BASE = Path(__file__).resolve().parents[1]
RESEARCH = BASE.parent
LIVE = RESEARCH / 'realtime_collective_observation_v1'
F10X = RESEARCH / 'collective_motion_v4'
SESSION = BASE / 'microscope_sessions/microscope_20261019_v1'
DEFAULT_OUT = SESSION / 'results/piv_v2_videos_20261020_v1'
os.environ.setdefault('F10X_OUTPUT', str(DEFAULT_OUT / 'backend_import'))
sys.path[:0] = [str(LIVE / 'src'), str(F10X / 'src')]
import av
import cv2
import numpy as np
from openpiv import pyprocess
from backends import Backend
import piv_display as display
from render_calibrated_video import overlay, sha, dump
from video_text import draw_text, font_record

PROTOCOL = dict(
    algorithm='Existing F10X B4: OpenPIV GPU FFT cross correlation, Gaussian subpixel peaks',
    window_px=96, overlap_px=48, past_frame_gap=3, max_pair_interval_s=.25,
    preprocessing='local: positive Gaussian(sigma=4 px) background minus original grayscale',
    peak_to_peak_snr_min=1.2, min_texture_std=2., local_median_tolerance_px=2.5,
    maximum_displacement_px=16.,
    candidate_support='At least one saved candidate in BOTH interrogation windows',
    coordinate_system='Original image: x right, y down; px/s using actual source PTS',
    geometry='Camera crop [0,40,1280,984]; no verified lumen mask',
    reference='IMAGE_SPACE_UNREGISTERED; includes any stage/camera motion',
    time_assignment='Interval-averaged displacement / actual PTS interval; displayed at interval end',
    invalid_vectors='Rejected, never filled or extrapolated',
    display=dict(grid_px=72, causal_time_constant_s=.12, gain_s=8.,
                 maximum_arrow_px=40, minimum_arrow_px=3, reference_speed_px_s=2.,
                 smoothing='Current valid cells only; no missing-cell hold; raw data unchanged'),
    training_ready=False, control_permission=False)


def configure():
    cv2.setNumThreads(2)
    # Process-local visualization settings, without editing the live interface.
    display.GAIN_SECONDS = PROTOCOL['display']['gain_s']


def source_frame(path, index=0):
    with av.open(str(path)) as c:
        s = c.streams.video[0]
        s.codec_context.thread_count = 2
        for i, f in enumerate(c.decode(s)):
            if i == index:
                return cv2.cvtColor(f.to_ndarray(format='bgr24')[40:, :1280],
                                    cv2.COLOR_BGR2GRAY)
    raise ValueError('Missing frame')


def preflight(out):
    configure()
    out.mkdir(parents=True, exist_ok=True)
    dump(out / 'protocol.json', PROTOCOL)
    source = Path('/home/lzy/projects/external_connect/sample_videos/v2/experiment_overlay_v2.mp4')
    a = source_frame(source)
    ws, ov = PROTOCOL['window_px'], PROTOCOL['overlap_px']
    gpu = Backend('B4', ws, ov)
    checks = []
    for dx, dy in [(0., 0.), (3., -2.), (-2.5, 1.5)]:
        b = cv2.warpAffine(a, np.float32([[1, 0, dx], [0, 1, dy]]),
                          (a.shape[1], a.shape[0]), flags=cv2.INTER_LINEAR,
                          borderMode=cv2.BORDER_REFLECT_101)
        z = gpu.run(a, b, .075, mode='local', do_fit=False, actual_coordinates=True)
        valid = z['valid'] & (z['x'] > ws) & (z['x'] < a.shape[1]-ws)
        valid &= (z['y'] > ws) & (z['y'] < a.shape[0]-ws)
        error = np.hypot(z['u'][valid]-dx, z['v'][valid]-dy)
        row = dict(shift_xy_px=[dx, dy], valid=int(valid.sum()),
                   median_estimated_xy_px=[float(np.median(z[k][valid])) for k in ['u', 'v']],
                   median_error_px=float(np.median(error)), p95_error_px=float(np.percentile(error, 95)))
        checks.append(row)
        assert valid.sum() >= 30 and row['median_error_px'] < .15 and row['p95_error_px'] < .5, row
        if dx == dy == 0:
            assert row['p95_error_px'] < .001, row
    b = source_frame(source, 3)
    z = gpu.run(a, b, 1, mode='local', do_fit=False, actual_coordinates=True)
    cpu = Backend('B1', ws, ov).run(a, b, 1, mode='local', do_fit=False, actual_coordinates=True)
    assert np.allclose(z['x'], cpu['x']) and np.allclose(z['y'], cpu['y'])
    both = z['valid'] & cpu['valid']
    err = np.hypot(z['u'][both]-cpu['u'][both], z['v'][both]-cpu['v'][both])
    agreement = dict(common_valid=int(both.sum()), median_displacement_difference_px=float(np.median(err)),
                     p95_displacement_difference_px=float(np.percentile(err, 95)))
    assert agreement['common_valid'] >= 30 and agreement['p95_displacement_difference_px'] < .2, agreement
    empty = gpu.run(np.zeros_like(a), np.zeros_like(a), 1, mode='local',
                    do_fit=False, actual_coordinates=True)
    assert not empty['valid'].any()
    report = dict(status='PASS', translations=checks, cpu_gpu=agreement,
                  blank_valid_windows=int(empty['valid'].sum()), source_sha256=sha(source),
                  interpretation='Numerical/sign checks, not physical ground truth for these recordings')
    dump(out / 'preflight.json', report)
    print(json.dumps(report), flush=True)


def support(centers, shape, x, y):
    points = np.zeros(shape, np.float32)
    q = np.rint(np.asarray(centers, float).reshape(-1, 2)-[0, 40]).astype(int)
    if len(q):
        ok = (q[:, 0] >= 0) & (q[:, 0] < shape[1]) & (q[:, 1] >= 0) & (q[:, 1] < shape[0])
        q = q[ok]
        np.add.at(points, (q[:, 1], q[:, 0]), 1)
    integral = cv2.integral(points, sdepth=cv2.CV_64F)
    ws = PROTOCOL['window_px']
    x0 = np.rint(x-ws/2).astype(int)
    y0 = np.rint(y-ws/2).astype(int)
    return (integral[y0+ws, x0+ws]-integral[y0, x0+ws]
            -integral[y0+ws, x0]+integral[y0, x0]) >= 1


def render_piv(image, candidates, vectors):
    objects = [{'xy_local': [x, y-40]} for x, y in candidates]
    result = overlay(image, objects, [0, 40, 1280, 984])
    camera = result[40:, :1280]
    arrows = display.piv_arrows(vectors, sample=False)
    display.draw_cv_arrows(camera, arrows)
    # Essential quantitative scale only; source timestamp and PZT are retained.
    draw_text(result, 'PIV', 740, 23, 18, (230, 235, 235))
    display.draw_cv_arrows(result[:30, :1280],
                          [display.arrow_shape([785, 16], [2., 0.], gain=8.)])
    draw_text(result, '2 px/s | cap 5', 812, 23, 18, (230, 235, 235))
    return result, len(arrows)


def run(number, root):
    configure()
    if json.loads((root / 'preflight.json').read_text())['status'] != 'PASS':
        raise ValueError('Preflight must pass first')
    if json.loads((root / 'protocol.json').read_text()) != PROTOCOL:
        raise ValueError('Frozen protocol mismatch')
    previous = SESSION / f'results/opencv_v2_videos_20261020_v1/experiment_overlay_v{number}_final'
    old = json.loads((previous / 'manifest.json').read_text())
    source = Path(old['source'])
    detections = previous / 'detections.jsonl.gz'
    assert old['status'] == 'COMPLETED' and sha(source) == old['source_sha256']
    assert sha(detections) == old['detections_sha256']
    out = root / f'experiment_overlay_v{number}'
    out.mkdir(exist_ok=False)
    started = time.monotonic()
    meta = dict(status='RUNNING', source=str(source), source_sha256=sha(source),
                candidate_video=old['output'], candidate_video_sha256=old['output_sha256'],
                detections_sha256=sha(detections), protocol=PROTOCOL, font=font_record(),
                parent_git_head=old['git_head'],
                git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=BASE, text=True).strip(),
                environment=dict(opencv=cv2.__version__, **{k: importlib.metadata.version(k)
                                 for k in ['av', 'openpiv', 'cupy-cuda12x']}))
    paths = [Path(__file__), BASE/'scripts/render_calibrated_video.py',
             F10X/'src/backends.py', LIVE/'src/piv_display.py', LIVE/'src/video_text.py']
    paths += list((F10X/'third_party/openpiv-python-gpu/openpiv_gpu').glob('*.py'))
    meta['source_hashes'] = {str(p): sha(p) for p in paths}
    snapshot = out / 'source_snapshot'
    snapshot.mkdir()
    for p in paths:
        shutil.copy2(p, snapshot/p.name)
    dump(out/'manifest.json', meta)
    ws, ov = PROTOCOL['window_px'], PROTOCOL['overlap_px']
    x, y = pyprocess.get_coordinates((984, 1280), ws, ov, center_on_field=False)
    n = old['decoded_input_frames']
    raw = np.full((n, *x.shape, 2), np.nan, np.float32)
    snr = np.full((n, *x.shape), np.nan, np.float32)
    valid = np.zeros((n, *x.shape), bool)
    selected = np.zeros_like(valid)
    times = np.zeros(n)
    dt_all = np.full(n, np.nan)
    partners = np.full(n, -1, np.int32)
    histories = deque(maxlen=3)
    filter_ = display.PIVDisplayFilter()
    backend = Backend('B4', ws, ov)
    stats = []
    expected = []
    durations = {}
    unchanged = True
    movie = out / f'experiment_overlay_v{number}_opencv_piv.mp4'
    partial = movie.with_suffix('.partial.mp4')
    with av.open(str(source)) as src, gzip.open(detections, 'rt') as labels, \
         gzip.open(out/'display_vectors.jsonl.gz', 'wt') as log, \
         av.open(str(partial), 'w', options={'movflags': '+faststart'}) as dst:
        stream = src.streams.video[0]
        stream.codec_context.thread_count = 2
        tb = stream.time_base
        encoder = dst.add_stream('libx264', rate=Fraction(30, 1))
        encoder.width, encoder.height, encoder.pix_fmt = 1964, 1024, 'yuv420p'
        encoder.time_base = encoder.codec_context.time_base = tb
        encoder.options = dict(crf='18', preset='fast', threads='2', bf='0', g='60')

        def mux(packet):
            pts = Fraction(packet.pts)*packet.time_base
            packet.duration = round(durations[pts]/packet.time_base)
            dst.mux(packet)

        origin = None
        for index, frame in enumerate(src.decode(stream)):
            assert index < n
            record = json.loads(next(labels))
            assert record['frame'] == index and record['source_pts'] == frame.pts
            absolute = Fraction(frame.pts)*frame.time_base
            if origin is None:
                origin = absolute
            pts = absolute-origin
            assert pts == Fraction(record['output_pts'])*Fraction(record['output_time_base'])
            assert not expected or pts > expected[-1]
            assert frame.duration
            expected.append(pts)
            durations[pts] = Fraction(frame.duration)*frame.time_base
            times[index] = t = float(pts)
            image = frame.to_ndarray(format='bgr24')
            gray = cv2.cvtColor(image[40:, :1280], cv2.COLOR_BGR2GRAY)
            current_support = support(record['centers_xy'], gray.shape, x, y)
            vectors = np.empty((0, 4))
            compute_ms = 0.
            dt = None
            reason = 'WARMUP'
            if len(histories) == 3:
                partner, prior_t, a, prior_support = histories[0]
                dt = t-prior_t
                if 0 < dt <= .25:
                    select = current_support & prior_support
                    selected[index] = select
                    dt_all[index], partners[index] = dt, partner
                    reason = 'NO_SUPPORTED_WINDOWS'
                    if select.any():
                        z = backend.run(a, gray, dt, mode='local', selection=select,
                                        do_fit=False, actual_coordinates=True)
                        assert np.allclose(x, z['x']) and np.allclose(y, z['y'])
                        raw[index, ..., 0], raw[index, ..., 1] = z['u'], z['v']
                        snr[index] = z['snr']
                        ok = z['valid'] & select & (np.hypot(z['u'], z['v']) <= 16.)
                        valid[index] = ok
                        vectors = np.c_[x[ok], y[ok], z['u'][ok]/dt, z['v'][ok]/dt]
                        compute_ms = z['total_ms']
                        reason = 'OK' if len(vectors) >= 8 else 'LOW_SUPPORT'
                else:
                    reason = 'PTS_GAP'
                    histories.clear()
            histories.append((index, t, gray, current_support))
            state = dict(valid=len(vectors) >= 8, vectors=vectors, arrival_monotonic=t,
                         geometry=None, session_id=f'OFFLINE_V{number}',
                         motion_reference='IMAGE_SPACE_UNREGISTERED')
            smoothed = filter_.update(state, ttl_ms=250)
            result, arrows = render_piv(image, record['centers_xy'], smoothed)
            unchanged &= np.array_equal(result[:, 1280:], image[:, 1280:])
            speed = np.linalg.norm(vectors[:, 2:], axis=1)
            row = dict(frame=index, source_pts=frame.pts, time_s=t, dt_s=dt,
                       partner_frame=int(partners[index]), reason=reason, candidates=record['count'],
                       supported_windows=int(selected[index].sum()), valid_windows=len(vectors),
                       displayed_arrows=arrows, median_speed_px_s=float(np.median(speed)) if len(speed) else None,
                       p95_speed_px_s=float(np.percentile(speed, 95)) if len(speed) else None,
                       piv_ms=compute_ms)
            stats.append(row)
            log.write(json.dumps(dict(**row, display_vectors_xyuv_local=smoothed.tolist()),
                                 allow_nan=False, separators=(',', ':'))+'\n')
            output = av.VideoFrame.from_ndarray(result, format='bgr24')
            output.pts, output.time_base = round(pts/tb), tb
            for packet in encoder.encode(output):
                mux(packet)
            if index in [0, n//2, n-1]:
                cv2.imwrite(str(out/f'preview_{index:05d}.jpg'), result)
            if index % 300 == 0:
                progress = dict(frame=index+1, total=n, valid_windows=len(vectors), arrows=arrows,
                                wall_s=time.monotonic()-started)
                dump(out/'progress.json', progress)
                print(json.dumps(progress), flush=True)
        assert len(expected) == n and next(labels, None) is None
        for packet in encoder.encode():
            mux(packet)
    np.savez_compressed(out/'piv_fields.npz', x_image_px=x, y_image_px=y+40,
                        displacement_px=raw, valid=valid, supported=selected, snr=snr,
                        time_s=times, dt_s=dt_all, partner_frame=partners)
    with (out/'frame_statistics.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=stats[0].keys())
        writer.writeheader()
        writer.writerows(stats)
    count = 0
    error = 0.
    with av.open(str(partial)) as c:
        s = c.streams.video[0]
        s.codec_context.thread_count = 2
        for f in c.decode(s):
            assert count < n
            error = max(error, float(abs(Fraction(f.pts)*f.time_base-expected[count])))
            count += 1
        duration = float(s.duration*s.time_base)
    validation = dict(decoded_frames=count, expected_frames=n, maximum_pts_error_s=error,
                      duration_s=duration, expected_duration_s=float(expected[-1]+durations[expected[-1]]),
                      pzt_pixels_unchanged_before_lossy_encode=bool(unchanged))
    dump(out/'validation.json', validation)
    assert count == n and error == 0 and unchanged
    assert abs(duration-validation['expected_duration_s']) < 1e-7
    partial.replace(movie)
    finite = [r for r in stats if r['valid_windows']]
    meta.update(status='COMPLETED', output=str(movie), validation=validation,
                wall_s=time.monotonic()-started,
                frames_with_arrows=sum(r['displayed_arrows'] > 0 for r in stats),
                frames_with_valid_piv=len(finite),
                median_of_frame_median_speed_px_s=float(np.median([r['median_speed_px_s'] for r in finite])) if finite else None,
                piv_ms_median=float(np.median([r['piv_ms'] for r in finite])) if finite else None,
                hashes={p.name: sha(p) for p in [movie, out/'piv_fields.npz',
                         out/'display_vectors.jsonl.gz', out/'frame_statistics.csv']})
    dump(out/'manifest.json', meta)
    print(json.dumps({k: meta[k] for k in ['status', 'output', 'wall_s', 'frames_with_arrows',
                                         'frames_with_valid_piv', 'median_of_frame_median_speed_px_s']}), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-root', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--preflight', action='store_true')
    parser.add_argument('--video', type=int, choices=[1, 2])
    args = parser.parse_args()
    if args.preflight:
        preflight(args.output_root)
    if args.video:
        run(args.video, args.output_root)
