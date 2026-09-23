"""Saved streamline provenance, independent cap geometry and P1 velocity checks."""
from pathlib import Path
import sys,json,hashlib
import numpy as np
import pyvista as pv
import vtk
from PIL import Image

ROOT=Path(__file__).resolve().parents[2]
CASE=ROOT/'flow_cases/mean-2p0-mmps';OUT=CASE/'streamlines'


def read(name):return json.loads((OUT/name).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def on_surface(point,surface):
    triangles=np.asarray(surface.points)[surface.faces.reshape(-1,4)[:,1:]]
    a=triangles[:,1]-triangles[:,0];b=triangles[:,2]-triangles[:,0];p=point-triangles[:,0]
    aa=np.einsum('ij,ij->i',a,a);bb=np.einsum('ij,ij->i',b,b);ab=np.einsum('ij,ij->i',a,b)
    ap=np.einsum('ij,ij->i',a,p);bp=np.einsum('ij,ij->i',b,p);den=aa*bb-ab*ab
    u=(ap*bb-bp*ab)/den;v=(bp*aa-ap*ab)/den
    normal=np.cross(a,b);distance=np.abs(np.einsum('ij,ij->i',p,normal))/np.linalg.norm(normal,axis=1)
    return np.any((distance<1e-12)&(u>=-1e-8)&(v>=-1e-8)&(u+v<=1+1e-8))


def test_correct_new_field_and_original_artifacts_unchanged():
    c=read('COMPUTE_VALIDATION.json');assert c['all_pass'] and c['inlet_mean_mm_s']==2.
    assert c['source_field_sha256']==sha(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    for file,expected in read('SOURCE_LOCK.json').items():assert sha(CASE/file)==expected
    for file,expected in c['outputs_sha256'].items():assert sha(OUT/file)==expected


def test_selected_curves_are_saved_forward_integrations_connecting_actual_caps():
    c=read('COMPUTE_VALIDATION.json');catalog=read('data/catalog.json')
    assert c['selected_count']==len(catalog)==96
    assert c['selected_outlet_counts']=={'OUTLET_01':16,'OUTLET_02':56,'OUTLET_03':24}
    paths=np.load(OUT/'data/selected_paths.npz');pool=np.load(OUT/'data/candidates.npz')
    caps={role:pv.read(CASE/'SV_MESH/mesh-surfaces'/f'{role}.vtp') for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']}
    for row in catalog:
        path=paths[f"line_{row['line_id']:03d}"]
        assert np.array_equal(path,pool[f"seed_{row['seed_id']:04d}"])
        assert np.isfinite(path).all() and np.all(np.diff(path[:,0])>0)
        assert on_surface(path[0,1:],caps['INLET'])
        assert on_surface(path[-1,1:],caps[row['outlet']])
        assert not any(on_surface(path[-1,1:],cap) for name,cap in caps.items() if name!=row['outlet'])


def test_saved_velocity_matches_independent_barycentric_interpolation():
    paths=np.load(OUT/'data/selected_paths.npz');points=[];velocities=[]
    for row in read('data/catalog.json'):
        key=f"{row['line_id']:03d}";path=paths['line_'+key]
        # Interior samples avoid implementation-specific boundary ownership tolerances.
        indices=np.unique(np.linspace(1,len(path)-2,20).astype(int))
        points.extend(path[indices,1:]);velocities.extend(paths['velocity_'+key][indices])
    field=pv.read(CASE/'frozen_flow/steady_flow_mean_2p0_mmps.vtu')
    xyz=np.asarray(field.points,dtype=float);tet=field.cells.reshape(-1,5)[:,1:]
    locator=vtk.vtkStaticCellLocator();locator.SetDataSet(field);locator.BuildLocator()
    sampled=[]
    for position in points:
        # VTK enumerates boxes only. Its permissive point-probe containment can
        # choose a neighboring tetra with negative barycentric weights.
        ids=vtk.vtkIdList();locator.FindCellsWithinBounds(np.column_stack((position-1e-13,position+1e-13)).ravel(),ids)
        candidates=np.array(sorted(ids.GetId(i) for i in range(ids.GetNumberOfIds())))
        nodes=tet[candidates];vertices=xyz[nodes]
        matrix=np.swapaxes(vertices[:,1:]-vertices[:,:1],1,2)
        tail=np.linalg.solve(matrix,(position-vertices[:,0])[...,None])[...,0]
        weights=np.column_stack((1-tail.sum(1),tail));valid=np.flatnonzero(weights.min(1)>=-1e-10)
        assert len(valid)>0
        k=valid[0];sampled.append(weights[k]@field['Velocity'][nodes[k]])
    assert np.max(np.linalg.norm(np.array(sampled)-velocities,axis=1))<1e-12


def test_refined_paths_match_geometry_and_final_endpoint():
    coarse=np.load(OUT/'data/selected_paths.npz');fine=np.load(OUT/'data/refined_paths.npz')
    sys.path.insert(0,str(ROOT/'scripts/flow_2mmps'))
    from compute_streamlines import resampled_path
    for row in read('data/catalog.json'):
        key=f"line_{row['line_id']:03d}";a,b=coarse[key],fine[key]
        assert np.linalg.norm(a[-1,1:]-b[-1,1:])<.05e-6
        assert np.max(np.linalg.norm(resampled_path(a)-resampled_path(b),axis=1))<.05e-6


def test_streamline_media_has_full_decode_and_unclipped_uniform_camera():
    media=read('MEDIA_VALIDATION.json');assert media['all_pass']
    video=media['video'];assert video['decoded_frames']==432 and video['distinct_frames']>=410
    assert sha(CASE/video['file'])==video['sha256']
    trace=read('CAMERA_TRACE.json');assert len(trace)==432
    assert np.allclose(np.diff([r['azimuth_deg'] for r in trace]),360/432)
    assert max(abs(x) for r in trace for x in r['projected_bounds'])<.94
    for r in media['figures']:
        assert sha(CASE/r['file'])==r['sha256']
        with Image.open(CASE/r['file']) as im:im.verify()
