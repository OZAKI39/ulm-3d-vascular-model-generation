# Executed only by the official SimVascular embedded Python (3.5 compatible).
import json
import os
import sv
import vtk

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
with open(os.path.join(ROOT, 'reports/mesh_and_flow/official_smoke.json')) as stream:
    assert json.load(stream)['status'] == 'PASS'
path = os.path.join(ROOT, 'outputs/mesh_and_flow/model/source.vtp')
reader = vtk.vtkXMLPolyDataReader()
reader.SetFileName(path)
reader.Update()
source = reader.GetOutput()
assert source.GetNumberOfCells() > 0
model = sv.modeling.PolyData()
model.set_surface(surface=source)
ids = sorted(model.get_face_ids())
assert ids == [1, 2, 3, 4, 5], ids
surface = model.get_polydata()
mesher = sv.meshing.TetGen()
mesher.set_model(model)
mesher.set_walls([1])
assert sorted(mesher.get_model_face_ids()) == ids
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName(os.path.join(ROOT, 'outputs/mesh_and_flow/model/imported.vtp'))
writer.SetInputData(surface)
assert writer.Write() == 1
with open(os.path.join(ROOT, 'inputs/fem_reference/port_contract.json')) as stream:
    roles = json.load(stream)['facet_names']
mapping = {'status': 'PASS', 'method': 'Original ModelFaceID labels preserved by official PolyData API',
           'faces': [{'sv_face_id': i, 'original_entity_id': i, 'role': roles[str(i)]} for i in ids]}
with open(os.path.join(ROOT, 'configs/face_map.json'), 'w') as stream:
    json.dump(mapping, stream, indent=2)
result = {'status': 'PASS', 'points': surface.GetNumberOfPoints(), 'triangles': surface.GetNumberOfCells(),
          'model_class': str(type(model)), 'mesher_class': str(type(mesher)),
          'face_ids': ids, 'mesh_generated': False, 'units': 'm'}
with open(os.path.join(ROOT, 'reports/mesh_and_flow/geometry_import.json'), 'w') as stream:
    json.dump(result, stream, indent=2)
print('SURFACE_MODEL_IMPORT=' + json.dumps(result))
