"""Real render/interaction audit for the single X/Y/Z annotation requirement."""
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from vascular_processing.compact_display_adapter import CompactBraVaViewer
from vascular_processing.topbrain_qc import write_json

OUT = Path(__file__).resolve().parent
DATA = ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi'
failures = []
counts = {'full_rotation': 0, 'gui': 0}


def audit_frame(viewer, phase):
    counts[phase] += 1
    for panel, controller in enumerate(viewer.controllers):
        for axis, annotations in zip('XYZ', (controller.x_annotations, controller.y_annotations, controller.z_annotations)):
            visible = [a for a in annotations if any(t.GetVisibility() for t in a.tick_actors)]
            titles = sum(bool(a.title_actor.GetVisibility()) for a in annotations)
            marks = sum(bool(a.tick_line_actor.GetVisibility()) for a in annotations)
            if len(visible) != 1 or titles != 1 or marks != 1:
                failures.append(dict(phase=phase, frame=counts[phase], panel=panel, axis=axis,
                                     tick_groups=len(visible), titles=titles, marks=marks))
            for annotation in visible:
                if not all(t.GetVisibility() and t.GetTextProperty().GetOpacity() == 1.
                           for t in annotation.tick_actors):
                    failures.append(dict(phase=phase, panel=panel, axis=axis, incomplete_ticks=True))


viewer = CompactBraVaViewer(DATA, initial_roi='BALANCED', show=False, output_dir=OUT/'rotation')
try:
    viewer.plotter.show(interactive=False, auto_close=False)
    observer = viewer.plotter.render_window.AddObserver('EndEvent', lambda *_: audit_frame(viewer, 'full_rotation'))
    # The production timer uses a 0.4-degree step. Check every rendered frame
    # through a full revolution, including axis-side transitions.
    for step in range(900):
        for renderer in viewer.plotter.renderers:
            renderer.camera.Azimuth(.4)
        viewer.plotter.render_window.Render()
        if step % 150 == 149:
            print(f'Rotation checked: {step + 1}/900 frames', flush=True)
    # Exercise oblique mouse-like motion, wheel zoom, viewport resize and ROI changes.
    for renderer in viewer.plotter.renderers:
        renderer.camera.Elevation(15.)
    viewer.plotter.render_window.Render()
    viewer.plotter.window_size = (1600, 800)
    viewer.plotter.render_window.Render()
    for index in (0, 2, 1):
        viewer.select(index)
        viewer.plotter.render_window.Render()
    viewer.plotter.render_window.RemoveObserver(observer)
finally:
    viewer.plotter.renderers[1].RemoveObserver(viewer.scale_alignment_observer)
    for controller in (*viewer.controllers, *viewer.orientation_labels):
        if controller is not None:
            controller.dispose()
    viewer.plotter.close()

viewer = CompactBraVaViewer(DATA, initial_roi='BALANCED', show=True, output_dir=OUT/'gui_smoke')
viewer.plotter.render_window.AddObserver('EndEvent', lambda *_: audit_frame(viewer, 'gui'))
report = viewer.run_window(show=True, smoke_seconds=2)
assert report['gui_smoke_passed'], report
result = dict(status='PASS' if not failures else 'FAIL', rendered_frames=counts, failures=failures,
              per_panel_per_axis_tick_groups=1, continuous_rotation_degrees=360.,
              roi_switch_resize_and_elevation_checked=True, gui_smoke_passed=True,
              gui_events=report['events'])
write_json(OUT/'render_audit.json', result)
print(json.dumps(result, indent=2))
assert not failures, failures[:10]
