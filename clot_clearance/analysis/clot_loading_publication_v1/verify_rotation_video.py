"""Decode the delivered movie and verify stream, camera path and stills."""
from pathlib import Path
import hashlib
import json
import sys
import argparse

import imageio_ffmpeg
import numpy as np
from PIL import Image

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'results/clot_loading_publication_v1/z_rotation'


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def main():
    movie = OUT/'Fig_clot_loading_Z_rotation.mp4'
    provenance = json.loads((OUT/'PROVENANCE.json').read_text())
    clean_revision = provenance.get('grid_face_count') == 6
    if clean_revision:
        assert all(face['grid_segments'] >= 9 for face in provenance['grid_faces'])
        assert provenance['layout']['page_clipping'] == 0 and not provenance['layout']['visible_text_overlaps']
        assert provenance['layout']['maximum_per_frame_opacity_change'] < .14
        assert provenance['layout']['maximum_per_frame_visible_label_motion_px'] < 35
        assert set(provenance['layout']['text_inventory']) == {
            'X (mm)', 'Y (mm)', 'Z (mm)', '-6.8', '-6.0', '-5.2', '-0.5', '0', '0.5',
            '0.6', '1', '1.4', '2 µN', 'Clot nodes', 'Fixed anchors', 'Surviving bonds',
            'Load force', 'Fixed load center'}
    reader = imageio_ffmpeg.read_frames(str(movie), pix_fmt='rgb24')
    meta = next(reader)
    assert meta['size'] == (1920, 1920)
    assert meta['fps'] == 30.
    assert abs(meta['duration']-12.) < .001
    assert meta['codec'] == 'h264' and 'yuv420p' in meta['pix_fmt']
    camera = json.loads((OUT/'CAMERA_PATH.json').read_text())
    if clean_revision:
        previous_camera = json.loads((OUT.parent/'z_rotation/CAMERA_PATH.json').read_text())
        assert camera == previous_camera
    eye = np.array([a['eye_mm'] for a in camera]); focal = np.array([a['focal_mm'] for a in camera])
    assert len(camera) == 360 and all(a['physical_step'] == 540 for a in camera)
    np.testing.assert_allclose(focal, np.tile(focal[0], (360, 1)), rtol=0, atol=1e-14)
    np.testing.assert_allclose(eye[:, 2], eye[0, 2], rtol=0, atol=1e-14)
    np.testing.assert_allclose(np.linalg.norm(eye-focal, axis=1), 20., rtol=0, atol=1e-13)
    elevations = np.rad2deg(np.arcsin((eye-focal)[:, 2]/20.))
    np.testing.assert_allclose(elevations, -20., rtol=0, atol=1e-12)
    np.testing.assert_allclose(np.unwrap(np.arctan2(eye[:, 1]-focal[:, 1], eye[:, 0]-focal[:, 0])),
                               np.deg2rad(np.arange(360))+np.arctan2(*(eye[0]-focal[0])[[1, 0]]),
                               rtol=0, atol=1e-13)
    sample_errors = []; hashes = []; consecutive_changes = []; previous = None
    first = None; last = None
    for i, raw in enumerate(reader):
        a = np.frombuffer(raw, np.uint8).reshape(1920, 1920, 3)
        if clean_revision:
            # These bands held the old title/angle, formula/parameters, and
            # footer. Confirm their removal in every encoded video frame.
            assert a[:170].min() > 246
            assert a[1590:1730].min() > 246
            assert a[1820:].min() > 246
        roi = np.ascontiguousarray(a[450:1150, 450:1450])
        hashes.append(hashlib.sha256(roi.tobytes()).hexdigest())
        if previous is not None:
            consecutive_changes.append(float(np.abs(roi.astype(np.int16)-previous.astype(np.int16)).mean()))
        previous = roi
        if i == 0:
            first = roi.copy()
        if i % 45 == 0 or (clean_revision and (OUT/f'view_{i:03d}.png').exists()):
            reference = np.asarray(Image.open(OUT/f'view_{i:03d}.png').convert('RGB'))
            error = a.astype(float)-reference.astype(float)
            mae = float(abs(error).mean()); psnr = float(10*np.log10(255**2/np.mean(error**2)))
            assert mae < 3 and psnr > 30, (i, mae, psnr)
            sample_errors.append(dict(frame=i, mean_absolute_encoding_error=mae, psnr_dB=psnr))
        last = roi
    count = i+1
    assert count == 360 and len(set(hashes)) == 360
    assert min(consecutive_changes) > .01
    seam = float(np.abs(last.astype(np.int16)-first.astype(np.int16)).mean())
    assert seam < 1.5*max(consecutive_changes)
    assert all(sha(Path(path)) == expected for path, expected in provenance['source_sha256'].items())
    manifest = json.loads((OUT/'DELIVERY_SHA256.json').read_text())
    assert all(sha(OUT/name) == record['sha256'] for name, record in manifest.items())
    report = dict(status='PASS', video=str(movie), bytes=movie.stat().st_size, sha256=sha(movie),
                  stream_metadata=meta, decoded_frames=count, unique_scene_frames=len(set(hashes)),
                  constant_camera_elevation_and_radius=True, constant_physical_step=540,
                  requested_side_view_elevation_deg=float(elevations[0]),
                  six_grid_faces_verified=clean_revision,
                  removed_text_bands_blank_in_all_frames=clean_revision,
                  constant_azimuth_increment_deg=1, source_files_unchanged=True,
                  consecutive_scene_change_range=[min(consecutive_changes), max(consecutive_changes)],
                  loop_seam_scene_difference=seam, selected_encoding_checks=sample_errors)
    (OUT/'VIDEO_VERIFICATION.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(status='PASS', frames=count, fps=meta['fps'], duration_s=meta['duration'],
                         resolution=meta['size'], bytes=movie.stat().st_size,
                         minimum_psnr_dB=min(a['psnr_dB'] for a in sample_errors)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--directory', type=Path, default=OUT)
    OUT = parser.parse_args().directory
    main()
