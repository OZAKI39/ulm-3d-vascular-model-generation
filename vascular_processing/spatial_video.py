"""Export the unchanged human dual-panel scene as a reproducible rotation MP4."""
from datetime import datetime, timezone
import time

from vtkmodules.vtkRenderingCore import vtkCamera

from utils.rodent_vasculature import interactive as ui
from .compact_static_views import camera_state
from .topbrain_qc import sha256, write_json


VIDEO_FPS = 25
VIDEO_DURATION_SECONDS = 36.


def export_spatial_video(viewer, *, fps=VIDEO_FPS, duration_seconds=VIDEO_DURATION_SECONDS):
    """Consume an independent off-screen viewer, retaining interactive framing.

    Fixed simulation time produces the requested duration even if rendering
    takes longer than real time. Geometry, selection, fonts, overlays and both
    camera focal points remain those of SpatialROIViewer. No GUI events or
    additional wheel steps are injected, and the GUI is constructed afterwards.
    """
    manifest = None
    run = None
    try:
        import imageio.v2 as imageio
        import imageio_ffmpeg

        if fps <= 0 or duration_seconds <= 0:
            raise ValueError('Video frame rate and duration must be positive')
        frame_count = round(fps * duration_seconds)
        if frame_count < 2:
            raise ValueError('A rotation video requires at least two frames')
        stamp = datetime.now(timezone.utc)
        parent = viewer.run.resolve() / 'videos' / viewer.initial_roi / viewer.view
        run = parent / stamp.strftime('%Y%m%dT%H%M%S_%fZ')
        run.mkdir(parents=True, exist_ok=False)
        path = run / 'human_spatial_rotation.mp4'
        partial = run / 'human_spatial_rotation.partial.mp4'
        speed = ui.HORIZONTAL_ROTATION_DEGREES_PER_STEP * 1000. / ui.HORIZONTAL_ROTATION_INTERVAL_MS
        width, height = map(int, viewer.plotter.window_size)
        manifest = dict(status='EXPORTING', created_utc=stamp.isoformat(), path=str(path),
                        source_provenance=viewer.static_provenance, roi=viewer.initial_roi,
                        format='MP4', codec='H.264', pixel_format='yuv420p', fps=fps,
                        duration_seconds=frame_count / fps, expected_frames=frame_count,
                        window_size=[width, height], rotation_axis='Y',
                        rotation_degrees_per_second=speed, rotation_degrees=speed * frame_count / fps,
                        framing='unchanged initial interactive cameras; no extra zoom',
                        display_units='mm', geometry_units='um', frames_written=0,
                        renderer=type(viewer).__name__, camera_checkpoints=[])
        write_json(run / 'manifest.json', manifest)
        print(f'SPATIAL_VIDEO_EXPORT: {path} ({frame_count / fps:g} s, {fps} fps)', flush=True)
        viewer.plotter.show(interactive=False, auto_close=False)
        for _ in range(3):
            viewer.plotter.render_window.Render()
        baseline = []
        for renderer in viewer.plotter.renderers:
            camera = vtkCamera()
            camera.DeepCopy(renderer.camera)
            baseline.append(camera)
        # Do not allow imageio's default 16-pixel rounding to change 1800x900.
        # Padding is only needed for an odd custom YAML resolution (H.264 4:2:0).
        output_params = ['-crf', '18', '-preset', 'fast', '-movflags', '+faststart']
        if width % 2 or height % 2:
            output_params += ['-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2']
        start = time.monotonic()
        with imageio.get_writer(str(partial), format='FFMPEG', mode='I', fps=fps,
                                codec='libx264', pixelformat='yuv420p', quality=None,
                                macro_block_size=1, ffmpeg_log_level='error',
                                output_params=output_params) as writer:
            for index in range(frame_count):
                angle = speed * index / fps
                for renderer, camera in zip(viewer.plotter.renderers, baseline):
                    renderer.camera.DeepCopy(camera)
                    renderer.camera.Azimuth(angle)
                    renderer.ResetCameraClippingRange()
                # A second render settles the projected triad/scale after the
                # camera change; this is the same annotation pipeline as PNG.
                for _ in range(2):
                    viewer.plotter.render_window.Render()
                frame = viewer.plotter.screenshot(return_img=True)
                writer.append_data(frame[:, :, :3])
                if index in (0, frame_count // 4, frame_count // 2, 3 * frame_count // 4, frame_count - 1):
                    manifest['camera_checkpoints'].append(dict(frame=index, time_seconds=index / fps,
                        azimuth_deg=angle, cameras=[camera_state(r.camera) for r in viewer.plotter.renderers]))
                if (index + 1) % (fps * 5) == 0 or index == frame_count - 1:
                    manifest['frames_written'] = index + 1
                    write_json(run / 'manifest.json', manifest)
                    print(f'SPATIAL_VIDEO_PROGRESS: {index + 1}/{frame_count} frames', flush=True)
        # Decode the finished stream before publishing the latest pointer.
        actual_frames, actual_seconds = imageio_ffmpeg.count_frames_and_secs(partial)
        if actual_frames != frame_count or abs(actual_seconds - frame_count / fps) > 1. / fps:
            raise RuntimeError(f'Incomplete video: {actual_frames} frames, {actual_seconds} seconds')
        with imageio.get_reader(str(partial), format='FFMPEG') as reader:
            metadata = reader.get_meta_data()
        expected_size = (width + width % 2, height + height % 2)
        if tuple(metadata['size']) != expected_size or abs(metadata['fps'] - fps) > .001:
            raise RuntimeError(f'Unexpected encoded size/frame rate: {metadata}')
        partial.replace(path)
        manifest.update(status='COMPLETE', frames_written=actual_frames,
                        encoded_duration_seconds=actual_seconds, encoded_size=list(metadata['size']),
                        bytes=path.stat().st_size, sha256=sha256(path),
                        encoding_wall_seconds=time.monotonic() - start)
        write_json(run / 'manifest.json', manifest)
        result = dict(status='COMPLETE', path=str(path), manifest=str(run / 'manifest.json'),
                      fps=fps, duration_seconds=actual_seconds, frames=actual_frames,
                      size=list(metadata['size']), bytes=path.stat().st_size)
        write_json(parent / 'latest.json', result)
        print(f'SPATIAL_VIDEO_COMPLETE: {path}', flush=True)
        return result
    except Exception as exc:
        if manifest is not None:
            manifest.update(status='FAILED', error=f'{type(exc).__name__}: {exc}')
            write_json(run / 'manifest.json', manifest)
        raise
    finally:
        viewer.close()
