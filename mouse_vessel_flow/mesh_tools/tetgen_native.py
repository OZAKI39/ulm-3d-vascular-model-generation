"""Original SimVascular TetGen body, compatible with embedded Python 3.5.

Executed by an explicit installed native launcher; not an old-workspace adapter.
"""
from pathlib import Path
import json, os, time

def validate_root(prepared_root):
    root = Path(prepared_root).resolve()
    project = Path(__file__).resolve().parents[2]
    protected = [project / n for n in ('inputs','cases','archive')]
    protected += [project / 'results' / n for n in ('roi_only','best_feasible','network_h0','mean_2mmps','stage_q')]
    if any(root == p or p in root.parents for p in protected):
        raise ValueError('Cannot mesh into retained scientific directories')
    return root

def generate(prepared_root, factor=1.0):
    root = validate_root(prepared_root)
    import sv
    import vtk
    ROOT = str(root)
    with open(os.path.join(ROOT, 'reports/mesh_and_flow/geometry_import.json')) as stream:
        assert json.load(stream)['status'] == 'PASS'
    with open(os.path.join(ROOT, 'configs/mesh_policy.json')) as stream:
        policy = json.load(stream)
    assert factor in (1.0, 0.8)
    if factor == 0.8:
        primary = os.path.join(ROOT, 'outputs/mesh_and_flow/mesh_generation/primary')
        assert os.path.isdir(primary), 'The primary generation must be attempted first'
        assert not os.path.exists(os.path.join(primary, 'volume.vtu')), 'No fallback after a generated primary volume'
    attempt = 'primary' if factor == 1.0 else 'fallback'
    output = os.path.join(ROOT, 'outputs/mesh_and_flow/mesh_generation', attempt)
    os.makedirs(output, exist_ok=True)
    assert not os.path.exists(os.path.join(output, 'volume.vtu'))
    mesher = sv.meshing.TetGen()
    mesher.load_model(os.path.join(ROOT, 'outputs/mesh_and_flow/model/imported.vtp'))
    mesher.set_walls([1])
    assert sorted(mesher.get_model_face_ids()) == [1, 2, 3, 4, 5]
    options = sv.meshing.TetGenOptions(global_edge_size=policy['h_ref_m'] * factor,
                                      surface_mesh_flag=True, volume_mesh_flag=True)
    options.use_mmg = False
    start = time.time()
    mesher.generate_mesh(options)
    mesh = mesher.get_mesh()
    assert mesh.GetNumberOfCells() > 0
    mesher.write_mesh(os.path.join(output, 'volume.vtu'))
    def write_surface(data, filename):
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(os.path.join(output, filename))
        writer.SetInputData(data)
        assert writer.Write() == 1
    write_surface(mesher.get_surface(), 'surface.vtp')
    for face_id in mesher.get_model_face_ids():
        write_surface(mesher.get_face_polydata(face_id), 'face_' + str(face_id) + '.vtp')
    record = {'status': 'PASS', 'attempt': attempt, 'edge_factor': factor,
              'global_edge_size_m': policy['h_ref_m'] * factor,
              'surface_meshing': True, 'volume_meshing': True,
              'points': mesh.GetNumberOfPoints(), 'cells': mesh.GetNumberOfCells(),
              'elapsed_s': time.time() - start, 'generator': 'sv.meshing.TetGen'}
    with open(os.path.join(output, 'generation.json'), 'w') as stream:
        json.dump(record, stream, indent=2)
    print('TETRA_MESH_GENERATION=' + json.dumps(record))
    return record

if __name__ == "__main__":
    generate(os.environ["MOUSE_MESH_ROOT"], float(os.environ.get("MESH_EDGE_FACTOR", "1.0")))
