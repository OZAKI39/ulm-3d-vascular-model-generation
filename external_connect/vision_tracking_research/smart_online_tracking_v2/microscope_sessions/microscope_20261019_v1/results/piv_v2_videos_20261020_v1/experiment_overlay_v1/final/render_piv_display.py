"""Final display pass from saved PIV measurements; no repeated PIV or inference."""
from fractions import Fraction
from pathlib import Path
import argparse
import gzip
import json
import shutil
import time
import av
import cv2
import numpy as np
import render_piv_video as p


def run(number):
    started = time.monotonic()
    parent = p.DEFAULT_OUT / f'experiment_overlay_v{number}'
    meta = json.loads((parent/'manifest.json').read_text())
    assert meta['status'] == 'COMPLETED'
    raw_path = parent/'piv_fields.npz'
    assert p.sha(raw_path) == meta['hashes']['piv_fields.npz']
    assert p.sha(meta['source']) == meta['source_sha256']
    archive = np.load(raw_path)
    raw, valid = archive['displacement_px'], archive['valid']
    x, y = archive['x_image_px'], archive['y_image_px']-40
    dt, times = archive['dt_s'], archive['time_s']
    out = parent/'final'
    out.mkdir(exist_ok=False)
    p.configure()
    p.display.GAIN_SECONDS = 24.
    p.display.MAX_LENGTH_PX = 48.
    filter_ = p.display.PIVDisplayFilter()
    movie = out/f'experiment_overlay_v{number}_opencv_piv.mp4'
    partial = movie.with_suffix('.partial.mp4')
    metadata = dict(status='RUNNING', source=meta['source'], source_sha256=meta['source_sha256'],
                    measurement_manifest_sha256=p.sha(parent/'manifest.json'),
                    raw_fields_sha256=p.sha(raw_path), script_sha256=p.sha(__file__),
                    gain_s=24., minimum_arrow_length_px=6., maximum_arrow_length_px=48.,
                    reference_speed_px_s=1., saturated_above_px_s=2.,
                    raw_data_changed=False, font=p.font_record(),
                    display_filter='Existing causal 72 px cell mean / 120 ms EMA; no missing-cell hold')
    p.dump(out/'manifest.json', metadata)
    shutil.copy2(__file__, out/'render_piv_display.py')
    n = len(times)
    expected, durations = [], {}
    unchanged = True
    arrow_counts = []
    labels = p.SESSION/f'results/opencv_v2_videos_20261020_v1/experiment_overlay_v{number}_final/detections.jsonl.gz'
    assert p.sha(labels) == meta['detections_sha256']
    with av.open(meta['source']) as src, gzip.open(labels, 'rt') as records, \
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
            stamp = Fraction(packet.pts)*packet.time_base
            packet.duration = round(durations[stamp]/packet.time_base)
            dst.mux(packet)

        origin = None
        for i, frame in enumerate(src.decode(stream)):
            assert i < n
            record = json.loads(next(records))
            absolute = Fraction(frame.pts)*frame.time_base
            if origin is None:
                origin = absolute
            stamp = absolute-origin
            assert record['frame'] == i and record['source_pts'] == frame.pts
            assert float(stamp) == times[i]
            expected.append(stamp)
            durations[stamp] = Fraction(frame.duration)*frame.time_base
            ok = valid[i]
            vectors = np.c_[x[ok], y[ok], raw[i, ..., 0][ok]/dt[i], raw[i, ..., 1][ok]/dt[i]]
            state = dict(valid=len(vectors) >= 8, vectors=vectors, arrival_monotonic=times[i],
                         geometry=None, session_id=f'PIV_FINAL_V{number}',
                         motion_reference='IMAGE_SPACE_UNREGISTERED')
            smooth = filter_.update(state, ttl_ms=250)
            displayed = smooth[np.linalg.norm(smooth[:, 2:], axis=1)*24 >= 6]
            image = frame.to_ndarray(format='bgr24')
            result = p.overlay(image, [dict(xy_local=[a, b-40]) for a, b in record['centers_xy']],
                               [0, 40, 1280, 984])
            arrows = p.display.piv_arrows(displayed, sample=False)
            p.display.draw_cv_arrows(result[40:, :1280], arrows)
            p.draw_text(result, 'PIV', 740, 23, 18, (230, 235, 235))
            p.display.draw_cv_arrows(result[:30, :1280],
                                    [p.display.arrow_shape([785, 16], [1., 0.], gain=24.)])
            p.draw_text(result, '1 px/s | cap 2', 820, 23, 18, (230, 235, 235))
            unchanged &= np.array_equal(result[:, 1280:], image[:, 1280:])
            arrow_counts.append(len(arrows))
            log.write(json.dumps(dict(frame=i, source_pts=frame.pts, dt_s=float(dt[i]) if np.isfinite(dt[i]) else None,
                                      time_s=times[i], displayed_xyuv_local=displayed.tolist()),
                                 allow_nan=False, separators=(',', ':'))+'\n')
            output = av.VideoFrame.from_ndarray(result, format='bgr24')
            output.pts, output.time_base = round(stamp/tb), tb
            for packet in encoder.encode(output):
                mux(packet)
            if i % 600 == 0:
                print(json.dumps(dict(frame=i+1, total=n, arrows=len(arrows),
                                      wall_s=time.monotonic()-started)), flush=True)
        assert len(expected) == n and next(records, None) is None
        for packet in encoder.encode():
            mux(packet)
    error, count = 0., 0
    with av.open(str(partial)) as c:
        s = c.streams.video[0]
        s.codec_context.thread_count = 2
        for frame in c.decode(s):
            assert count < n
            error = max(error, float(abs(Fraction(frame.pts)*frame.time_base-expected[count])))
            if count in [n//4, n//2, n-1]:
                cv2.imwrite(str(out/f'preview_{count:05d}.png'), frame.to_ndarray(format='bgr24'))
            count += 1
        duration = float(s.duration*s.time_base)
    assert count == n and error == 0 and unchanged
    assert abs(duration-float(expected[-1]+durations[expected[-1]])) < 1e-7
    assert p.sha(raw_path) == metadata['raw_fields_sha256']
    partial.replace(movie)
    metadata.update(status='COMPLETED', output=str(movie), output_sha256=p.sha(movie),
                    frames=count, duration_s=duration, maximum_pts_error_s=error,
                    pzt_pixels_unchanged_before_lossy_encode=bool(unchanged),
                    raw_fields_unchanged=True, wall_s=time.monotonic()-started,
                    frames_with_arrows=sum(a > 0 for a in arrow_counts),
                    arrow_count_median=float(np.median(arrow_counts)),
                    display_vectors_sha256=p.sha(out/'display_vectors.jsonl.gz'))
    p.dump(out/'manifest.json', metadata)
    print(json.dumps(metadata), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--video', type=int, choices=[1, 2], required=True)
    run(parser.parse_args().video)
