from stage03_helpers import *
from fem3d.audit import sha256

def test_selected_source_and_all_cell_quality_are_verified():
    c=config();m=preflight()['mesh']
    assert m['source_sha256']==c['mesh']['volume_mesh_sha256']==sha256(ROOT/c['mesh']['source']/'mesh/volume_mesh.npz')
    check=read('outputs/stage03/mesh/qc.json')
    assert check['status']=='PASS'
    assert m['cells']==len(np.load(ROOT/c['mesh']['source']/'mesh/volume_mesh.npz')['tetra'])
    assert m['inlets']==1 and m['outlets']==3 and m['unlabelled_exterior_facets']==0
    assert set(m['boundary_counts'])==set('12345') and m['cell_tag']==100
    assert m['geometry_modified'] is False
