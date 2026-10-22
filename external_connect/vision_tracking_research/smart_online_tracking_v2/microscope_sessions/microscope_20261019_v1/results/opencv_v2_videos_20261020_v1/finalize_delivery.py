"""Verify deliverable identities and extract decoded previews without buffering a video."""
from pathlib import Path
from fractions import Fraction
from itertools import zip_longest
import gzip
import hashlib
import json
import av
import cv2

ROOT = Path(__file__).resolve().parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    results = []
    preview = ROOT / 'visual_inspection'
    preview.mkdir(exist_ok=True)
    for number in (1, 2):
        folder = ROOT / f'experiment_overlay_v{number}_final'
        manifest = json.loads((folder / 'manifest.json').read_text())
        assert manifest['status'] == 'COMPLETED'
        for file_key, hash_key in [('source', 'source_sha256'), ('output', 'output_sha256'),
                                   ('profile_path', 'profile_sha256')]:
            assert sha(manifest[file_key]) == manifest[hash_key], file_key
        assert sha(folder / 'detections.jsonl.gz') == manifest['detections_sha256']
        assert manifest['decoded_input_frames'] == manifest['decoded_output_frames']
        assert manifest['maximum_pts_error_s'] == 0
        assert manifest['pzt_pixels_unchanged_before_lossy_encode']
        assert manifest['expected_duration_s'] == manifest['output_duration_s']
        count = manifest['decoded_output_frames']
        selected = {0, count // 2, count - 1}
        if number == 2:
            selected.add(2751)  # The previously rejected circle at the camera edge.
        with gzip.open(folder / 'detections.jsonl.gz', 'rt') as f:
            rows = [json.loads(line) for i, line in enumerate(f) if i in selected]
        previews = []
        with av.open(manifest['output']) as c:
            stream = c.streams.video[0]
            stream.codec_context.thread_count = 2
            for row in rows:
                target = Fraction(row['output_pts']) * Fraction(row['output_time_base'])
                c.seek(int(target / stream.time_base), stream=stream, backward=True)
                for frame in c.decode(stream):
                    stamp = Fraction(frame.pts) * frame.time_base
                    if stamp < target:
                        continue
                    assert stamp == target
                    name = preview / f'v{number}_decoded_{row["frame"]:04d}.jpg'
                    assert cv2.imwrite(str(name), frame.to_ndarray(format='bgr24'),
                                       [cv2.IMWRITE_JPEG_QUALITY, 95])
                    previews.append(str(name.relative_to(ROOT)))
                    break
                else:
                    raise RuntimeError('Preview frame missing')
        results.append(dict(video=manifest['output'], sha256=manifest['output_sha256'],
                            frames=count, duration_s=manifest['output_duration_s'],
                            source_unchanged=True, maximum_pts_error_s=0,
                            pzt_overlay_preserved=True, previews=previews,
                            candidate_counts=manifest['candidate_counts'],
                            processing_wall_s=manifest['wall_s']))
    old = ROOT / 'experiment_overlay_v2/detections.jsonl.gz'
    new = ROOT / 'experiment_overlay_v2_final/detections.jsonl.gz'
    compared = 0
    with gzip.open(old, 'rt') as a, gzip.open(new, 'rt') as b:
        for previous, final in zip_longest(a, b):
            assert previous is not None and final is not None
            assert json.loads(previous) == json.loads(final)
            compared += 1
    delivery = dict(status='VERIFIED', outputs=results, source_timing_preserved=True,
                    repair_changed_detections=False, unchanged_detection_frames=compared,
                    profile_use='Explicitly requested offline candidate; realtime default unchanged',
                    channel_mask='NOT_DRAWN', training_ready=False,
                    scope='Unreviewed candidates across the camera area, not verified lumen bubble counts')
    (ROOT / 'delivery.json').write_text(json.dumps(delivery, indent=2))
    print(json.dumps(delivery, indent=2))


if __name__ == '__main__':
    main()
