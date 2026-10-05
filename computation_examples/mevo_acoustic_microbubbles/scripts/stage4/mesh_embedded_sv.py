# Official SimVascular 2023.05 embedded Python 3.5, using the existing TetGen workflow.
import json
import os
import time
import sv
import vtk

with open(os.environ['MEVO_STAGE4_MESH_REQUEST']) as stream:
    request = json.load(stream)
out = request['output']
os.makedirs(out, exist_ok=True)
assert not os.path.exists(os.path.join(out, 'volume.vtu'))
mesher = sv.meshing.TetGen()
mesher.load_model(request['surface'])
mesher.set_walls([1])
assert sorted(mesher.get_model_face_ids()) == [1, 2, 3, 4]
options = sv.meshing.TetGenOptions(global_edge_size=request['edge_size_m'],
                                 surface_mesh_flag=True, volume_mesh_flag=True)
options.use_mmg = False
started = time.time()
mesher.generate_mesh(options)
mesh = mesher.get_mesh()
assert mesh.GetNumberOfCells() > 0
mesher.write_mesh(os.path.join(out, 'volume.vtu'))
for face_id in mesher.get_model_face_ids():
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(os.path.join(out, 'face_' + str(face_id) + '.vtp'))
    writer.SetInputData(mesher.get_face_polydata(face_id))
    assert writer.Write() == 1
report = dict(generator='SimVascular 2023.05 sv.meshing.TetGen', points=mesh.GetNumberOfPoints(),
              cells=mesh.GetNumberOfCells(), elapsed_s=time.time()-started,
              global_edge_size_m=request['edge_size_m'], status='GENERATED_PENDING_AUDIT')
with open(os.path.join(out, 'GENERATION.json'), 'w') as stream:
    json.dump(report, stream, indent=2)
print(json.dumps(report))
