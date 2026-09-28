"""Shared multi-angle exports for compact ROIs and the saved manufacturing STL."""
from datetime import datetime, timezone

from vtkmodules.vtkRenderingCore import vtkCamera

from .compact_display_adapter import CompactBraVaViewer
from .topbrain_qc import sha256, write_json


STATIC_AZIMUTHS_DEG = (0, 60, 120, 180, 240, 300)


def camera_state(camera):
    return dict(position=camera.GetPosition(), focal_point=camera.GetFocalPoint(),
                view_up=camera.GetViewUp(), view_angle_deg=camera.GetViewAngle(),
                parallel_projection=bool(camera.GetParallelProjection()),
                parallel_scale=camera.GetParallelScale(), clipping_range=camera.GetClippingRange())


def export_compact_static_views(directory, *, initial_roi=None, view='strict', output_dir=None):
    """Export compact SWC/surface data with the existing viewer."""
    viewer = CompactBraVaViewer(directory, initial_roi=initial_roi, view=view,
                               show=False, output_dir=output_dir)
    return _export_static_views(viewer)


def export_manufacturing_static_views(path, *, initial_roi=None, output_dir=None):
    """Export the actual saved STL with the identical camera/annotation workflow."""
    from .manufacturing_stl_display import ManufacturingSTLViewer

    viewer = ManufacturingSTLViewer(path, initial_roi=initial_roi,
                                    show=False, output_dir=output_dir)
    return _export_static_views(viewer)


def _export_static_views(viewer):
    """Save full dual-panel PNGs and close the supplied off-screen viewer.

    The disposable off-screen viewer has no rotation timer. The interactive
    viewer is constructed separately afterwards, so neither its camera nor its
    callbacks/selection state can be changed by this export.
    """
    try:
        view = viewer.view
        stamp = datetime.now(timezone.utc)
        parent = viewer.run.resolve()/'static_views'/viewer.initial_roi/view
        run = parent/stamp.strftime('%Y%m%dT%H%M%S_%fZ')
        run.mkdir(parents=True, exist_ok=False)
        manifest = dict(status='EXPORTING', created_utc=stamp.isoformat(),
                        output_directory=str(run), data_source=str(viewer.dataset_dir),
                        roi=viewer.initial_roi, view=view, format='PNG',
                        window_size=list(viewer.plotter.window_size),
                        angle_reference='Azimuth about Y, relative to the existing initial camera in each panel',
                        renderer=f'{type(viewer).__name__}; existing dual-panel style unchanged',
                        mouse_wheel_backward_steps_per_panel=1,
                        display_units='mm', geometry_units='um', images=[])
        if hasattr(viewer, 'stl_provenance'):
            manifest['stl_provenance'] = viewer.stl_provenance
        if hasattr(viewer, 'static_provenance'):
            manifest['source_provenance'] = viewer.static_provenance
        write_json(run/'manifest.json', manifest)
        print(f'STATIC_VIEWS_EXPORT: {run}', flush=True)
        try:
            # Initialize exactly the same overlays as run_window, without
            # starting an interactor or creating a visible window.
            for _ in range(3):
                viewer.plotter.render_window.Render()
            viewer.plotter.show(interactive=False, auto_close=False)
            # Use the actual UI wheel handler: one zoom-out notch in each
            # static panel, once before copying the cameras for all angles.
            # The separately constructed interactive viewer is unaffected.
            interactor = viewer.plotter.iren.interactor
            for renderer in viewer.plotter.renderers:
                x, y = renderer.GetOrigin()
                width, height = renderer.GetSize()
                interactor.SetEventPosition(x + width // 2, y + height // 2)
                interactor.MouseWheelBackwardEvent()
            baseline = []
            for renderer in viewer.plotter.renderers:
                camera = vtkCamera()
                camera.DeepCopy(renderer.camera)
                baseline.append(camera)
            for angle in STATIC_AZIMUTHS_DEG:
                for renderer, camera in zip(viewer.plotter.renderers, baseline):
                    renderer.camera.DeepCopy(camera)
                    if angle:
                        renderer.camera.Azimuth(angle)
                        # Keep the same framing/zoom but prevent clipping after rotation.
                        renderer.ResetCameraClippingRange()
                # Projected axes, triad letters and scalar-bar text settle after rendering.
                for _ in range(3):
                    viewer.plotter.render_window.Render()
                prefix = (viewer.static_filename_prefix if hasattr(viewer, 'static_filename_prefix')
                          else f'BG001_{viewer.side}_{viewer.initial_roi}_{view}')
                path = run/f'{prefix}_azimuth_{angle:03d}.png'
                viewer.plotter.screenshot(str(path))
                manifest['images'].append(dict(path=str(path), azimuth_deg=angle,
                    sha256=sha256(path), bytes=path.stat().st_size,
                    cameras=[camera_state(r.camera) for r in viewer.plotter.renderers]))
            manifest['status'] = 'COMPLETE'
            write_json(run/'manifest.json', manifest)
            result = dict(status='COMPLETE', directory=str(run), manifest=str(run/'manifest.json'),
                          count=len(manifest['images']), azimuths_deg=list(STATIC_AZIMUTHS_DEG))
            # Publish the latest pointer only once all six images are saved.
            write_json(parent/'latest.json', result)
            return result
        except Exception as exc:
            manifest.update(status='FAILED', error=f'{type(exc).__name__}: {exc}')
            write_json(run/'manifest.json', manifest)
            raise
    finally:
        viewer.plotter.renderers[1].RemoveObserver(viewer.scale_alignment_observer)
        for controller in (*viewer.controllers, *viewer.orientation_labels):
            if controller is not None:
                controller.dispose()
        viewer.plotter.close()
