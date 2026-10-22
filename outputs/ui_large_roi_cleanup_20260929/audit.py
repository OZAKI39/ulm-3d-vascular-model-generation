"""Inspect actual STL actors and every rendered GUI frame after overlay cleanup."""
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from vascular_processing.manufacturing_stl_display import ManufacturingSTLViewer
from vascular_processing.compact_static_views import camera_state
from vascular_processing.topbrain_qc import write_json, sha256

OUT = Path(__file__).resolve().parent
STL = ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/manufacturing_roi/BG001_RMCA_print_candidate.stl'
before = json.loads((OUT/'before/state.json').read_text())
viewer = ManufacturingSTLViewer(STL, show=True, output_dir=OUT/'gui_smoke')
# Match the captured baseline's post-render state; VTK updates the surface
# mapper's upstream normals pipeline lazily on its first render.
for _ in range(3):
    viewer.plotter.render_window.Render()
assert json.loads(json.dumps([camera_state(r.camera) for r in viewer.plotter.renderers])) == before['cameras']
surface = viewer.stl_actor.mapper.dataset
for name, array in [('point', surface.points), ('face', surface.faces), ('scalar', surface.point_data['diameter_um'])]:
    assert hashlib.sha256(np.asarray(array).tobytes()).hexdigest() == before[f'surface_{name}_sha256']
assert sha256(STL) == before['stl_provenance']['stl_sha256']

def scene_properties(renderer):
    items = []
    for actor in renderer.actors.values():
        if not actor.IsA('vtkActor') or not actor.GetMapper():
            continue
        mesh = actor.GetMapper().GetInput()
        if hasattr(mesh, 'GetNumberOfPoints'):
            items.append(dict(color=list(actor.GetProperty().GetColor()), opacity=actor.GetProperty().GetOpacity(),
                points=mesh.GetNumberOfPoints(), cells=mesh.GetNumberOfCells(),
                scalar_visibility=bool(actor.GetMapper().GetScalarVisibility())))
    return items

blue = np.array([81, 229, 255]) / 255.
for panel in (0, 1):
    expected = [{k:v for k,v in a.items() if k != 'name'} for a in before['geometry'][str(panel)]
                if not (panel == 1 and np.allclose(a['color'], blue) and not a['scalar_visibility'])]
    actual = scene_properties(viewer.plotter.renderers[panel])
    # Actor dictionary insertion order is not a rendering property.
    assert sorted(actual, key=lambda a: json.dumps(a,sort_keys=True)) == sorted(expected, key=lambda a: json.dumps(a,sort_keys=True)), (panel, actual, expected)

frames = 0
failures = []
def audit_frame(*_):
    global frames
    frames += 1
    left, right = viewer.plotter.renderers[0], viewer.plotter.renderers[1]
    texts = {name:actor.GetInput() for name,actor in left.actors.items() if name.startswith('full_scene_legend_text_')}
    if list(texts.values()) != ['Centerline orientation', 'Selected ROI vessels', 'ROI bounding box']:
        failures.append(dict(frame=frames, left_legend=texts))
    symbol = left.actors['full_scene_legend_symbol_1']
    if symbol.GetProperty().GetColor() != (1.,0.,0.) or symbol.GetMapper().GetInput().GetNumberOfLines() != 1:
        failures.append(dict(frame=frames, selected_vessel_symbol_invalid=True))
    if 'sampling_roi_centerline' in right.actors or any(a.GetInput() == 'ROI centerline'
            for n,a in right.actors.items() if n.startswith('sampling_roi_legend_text_')):
        failures.append(dict(frame=frames, centerline_or_legend_present=True))
    if any(np.allclose(a['color'],blue) and not a['scalar_visibility'] for a in scene_properties(right)):
        failures.append(dict(frame=frames, blue_overlay_present=True))
    if viewer.stl_actor.mapper.dataset.n_cells != 504984:
        failures.append(dict(frame=frames, wrong_surface=True))
    for panel,controller in enumerate(viewer.controllers):
        for axis,group in zip('XYZ',(controller.x_annotations,controller.y_annotations,controller.z_annotations)):
            if sum(any(t.GetVisibility() for t in a.tick_actors) for a in group) != 1:
                failures.append(dict(frame=frames,panel=panel,axis=axis,duplicate_ticks=True))

viewer.plotter.render_window.AddObserver('EndEvent', audit_frame)
# Exercise legend reflow before the actual timer/key smoke test as well.
viewer.plotter.window_size = (1600,800)
viewer.plotter.render_window.Render()
viewer.plotter.window_size = (1800,900)
viewer.plotter.render_window.Render()
report = viewer.run_window(show=True, smoke_seconds=2)
result = dict(status='PASS' if not failures and report['gui_smoke_passed'] else 'FAIL',
    rendered_frames=frames, failures=failures, gui_smoke_passed=report['gui_smoke_passed'],
    cameras_unchanged=True, stl_points_faces_and_colour_scalars_unchanged=True,
    other_scene_actor_properties_unchanged=True, source_stl_unchanged=sha256(STL)==before['stl_provenance']['stl_sha256'],
    no_interpolated_blue_overlay=True, left_legend_matches_red_vessels_and_box=True,
    resize_roi_refresh_and_keyboard_checked=True)
write_json(OUT/'gui_audit.json', result)
print(json.dumps(result,indent=2))
assert not failures and report['gui_smoke_passed']
