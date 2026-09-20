import numpy as np
import pyvista as pv
from fem_sync_support import ROOT, read, digest


def test_actual_face_ids_and_complete_surface_partition():
    b=read('frozen_reference/boundary_manifest.json')
    faces={f['role']:f['sv_face_id'] for f in read('configs/face_map.json')['faces']}
    assert set(b['boundaries'])=={'WALL','INLET','OUTLET_01','OUTLET_02','OUTLET_03'}
    assert b['complete_disjoint_surface_partition']
    assert sum(x['facets'] for x in b['boundaries'].values())==read('frozen_reference/mesh_manifest.json')['surface_facets']
    for name,f in b['boundaries'].items():
        assert f['sv_face_id']==faces[name] and digest(ROOT/f['path'])==f['sha256']


def test_outward_normals_verified_against_actual_tetrahedra():
    b=read('frozen_reference/boundary_manifest.json')['boundaries']
    g=pv.read(ROOT/read('frozen_reference/mesh_manifest.json')['path'])
    tetra=g.cells.reshape(-1,5)[:,1:]
    for f in b.values():
        s=pv.read(ROOT/f['path']);xyz=s.points[s.faces.reshape(-1,4)[:,1:]]
        n=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
        owner=g.points[tetra[s['GlobalElementID']-1]].mean(1)
        assert np.all(np.sum(n*(xyz.mean(1)-owner),axis=1)>0)
        assert f['outward_facets']==s.n_cells and not f['normal_array_stored']
        assert np.isclose(np.linalg.norm(n,axis=1).sum()/2,f['area_m2'])
    assert b['INLET']['signed_flux_m3_s']<0
    assert all(b[n]['signed_flux_m3_s']>0 for n in b if n.startswith('OUTLET'))
