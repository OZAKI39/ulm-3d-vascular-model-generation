"""Single native TetGen operation. Compatible with installed SV Python 3.5."""
import hashlib
import json
import os
import sys
import time
import shutil

import numpy as np
import sv
import vtk
from vtk.util.numpy_support import vtk_to_numpy

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
case = os.environ.get('BRAVA_PIPE_CASE', 'R1_P0')
source = os.path.join(root, 'straight_pipe', case)
out = os.path.join(root, os.environ['BRAVA_PIPE_OUTPUT_REL'])
os.makedirs(out, exist_ok=False)
config = json.load(open(os.path.join(source, 'INPUT.json')))
assert case == 'R1_P0' or os.environ.get('BRAVA_PIPE_CALIBRATION_AUTHORIZED') == '1'
policy = config['policy']
surface_remesh = os.environ.get('BRAVA_PIPE_SURFACE_REMESH', '1' if policy['layers'] else '0') == '1'
mesher = sv.meshing.TetGen()
mesher.load_model(os.path.join(source, 'model.vtp'))
set_walls = os.environ.get('BRAVA_PIPE_SET_WALLS', '1') == '1'
if set_walls:
    mesher.set_walls([1])
assert sorted(mesher.get_model_face_ids()) == [1, 2, 3]
options = sv.meshing.TetGenOptions(global_edge_size=config['h_bulk_m'],
                                surface_mesh_flag=surface_remesh, volume_mesh_flag=True)
options.use_mmg = os.environ.get('BRAVA_PIPE_USE_MMG', '1') == '1'
assert bool(options.surface_mesh_flag) == surface_remesh
assert bool(options.volume_mesh_flag)
if policy['layers']:
    assert surface_remesh, 'Installed native layer path requires surface remeshing'
    mesher.set_boundary_layer_options(number_of_layers=policy['layers'],
                                     edge_size_fraction=policy['edge_size_fraction'],
                                     layer_decreasing_ratio=policy['layer_decreasing_ratio'],
                                     constant_thickness=policy['constant_thickness'])
shutil.copyfile(__file__, os.path.join(out, 'executed_adapter.py'))
api = dict(python=sys.version, vtk=vtk.vtkVersion.GetVTKVersion(), sv_home=os.environ.get('SV_HOME'),
           requested_policy=policy, boundary_layer_options_doc=mesher.set_boundary_layer_options.__doc__,
           set_walls_called=set_walls, set_walls_doc=mesher.set_walls.__doc__,
           surface_mesh_flag_doc=sv.meshing.TetGenOptions.surface_mesh_flag.__doc__,
           volume_mesh_flag_doc=sv.meshing.TetGenOptions.volume_mesh_flag.__doc__,
           no_bisect_doc=sv.meshing.TetGenOptions.no_bisect.__doc__,
           actual_options={str(k): str(v) for k, v in options.get_values().items()},
           face_ids_before=mesher.get_model_face_ids(), model_sha256=hashlib.sha256(open(os.path.join(source, 'model.vtp'), 'rb').read()).hexdigest(),
           adapter_sha256=hashlib.sha256(open(__file__, 'rb').read()).hexdigest())
with open(os.path.join(out, 'INSTALLED_API_AND_INPUT.json'), 'w') as stream:
    json.dump(api, stream, indent=2)
print('NATIVE_PIPE_INPUT=' + json.dumps(api), flush=True)
start = time.time()


def write_surface(path, polydata):
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(polydata)
    assert writer.Write() == 1


try:
    mesher.generate_mesh(options)
except Exception as error:
    failure = dict(status='NATIVE_GENERATION_FAILED', error=str(error), elapsed_s=time.time() - start)
    try:
        surf = mesher.get_surface()
        write_surface(os.path.join(out, 'failed_native_surface.vtp'), surf)
        failure['surface_recovered'] = True
    except Exception as surface_error:
        failure['surface_recovered'] = False
        failure['surface_error'] = str(surface_error)
    with open(os.path.join(out, 'GENERATION.json'), 'w') as stream:
        json.dump(failure, stream, indent=2)
    raise
mesh = mesher.get_mesh()
mesher.write_mesh(os.path.join(out, 'volume.vtu'))
write_surface(os.path.join(out, 'surface.vtp'), mesher.get_surface())
for face, name in [(1, 'WALL'), (2, 'INLET'), (3, 'OUTLET_01')]:
    write_surface(os.path.join(out, name + '.vtp'), mesher.get_face_polydata(face))
types, counts = np.unique(vtk_to_numpy(mesh.GetCellTypesArray()), return_counts=True)
record = dict(status='GENERATED_PENDING_QA', case=case, source_model=source,
              surface_remesh=surface_remesh, elapsed_s=time.time() - start,
              points=mesh.GetNumberOfPoints(), cells=mesh.GetNumberOfCells(),
              cell_types={str(int(k)): int(v) for k, v in zip(types, counts)},
              face_ids_after=mesher.get_model_face_ids(), input_identity=api)
with open(os.path.join(out, 'GENERATION.json'), 'w') as stream:
    json.dump(record, stream, indent=2)
print('NATIVE_PIPE_RESULT=' + json.dumps(record), flush=True)
