"""The accepted four-port core and its upstream printing provenance, in SI units."""
import numpy as np
from .paths import ROOT, sha256, read_json, fresh_directory, write_json
from .transforms import source_mm_to_cfd_m

def describe(root=ROOT):
    data = read_json(root/'inputs/geometry/geometry.json')
    return dict(model='BraVa BG001 RMCA', source_geometry=data['geometry'],
                source_units='mm', computation_units='m', ports=data['ports'],
                candidate_0_and_15='Rigid poses of the same physical solution',
                source_sha256=sha256(root/'inputs/geometry/source_core_mm.stl'))

def verify_geometry(root=ROOT):
    import trimesh
    record = read_json(root/'inputs/geometry/geometry.json')
    mesh = trimesh.load_mesh(root/'inputs/geometry/source_core_mm.stl', process=True)
    assert sha256(root/'inputs/geometry/source_core_mm.stl') == record['sha256']
    assert mesh.is_watertight and mesh.is_winding_consistent and len(mesh.split()) == 1
    with np.load(root/'inputs/geometry/source_surface.npz') as data:
        assert np.array_equal(source_mm_to_cfd_m(mesh.vertices, root), data['points_m'])
        assert np.array_equal(mesh.faces, data['triangles'])
        assert set(np.unique(data['facet_tags'])) == {1,2,3,4,5}
    assert {k:v['face_id'] for k,v in record['ports'].items()} == {'I1':4,'O1':3,'O2':5,'O3':2}
    from .transforms import cfd_m_to_candidate_mm, rotate_vectors
    poses={}
    for candidate in (0,15):
        original=read_json(root/f'inputs/geometry/candidate_{candidate}_print_transform.json')
        matrix=np.asarray(record['pose_transforms'][str(candidate)]['source_mm_to_candidate_mm'])
        np.testing.assert_array_equal(matrix,original['transform_4x4'])
        np.testing.assert_allclose(matrix[:3,:3]@matrix[:3,:3].T,np.eye(3),rtol=0,atol=1e-14)
        np.testing.assert_allclose(np.linalg.det(matrix[:3,:3]),1.,rtol=0,atol=1e-14)
        displayed=cfd_m_to_candidate_mm(source_mm_to_cfd_m(mesh.vertices,root),candidate,root)
        expected=mesh.vertices@matrix[:3,:3].T+matrix[:3,3]
        np.testing.assert_allclose(displayed,expected,rtol=0,atol=1e-12)
        restored=np.c_[displayed,np.ones(len(displayed))]@np.asarray(original['casting_restore_transform']).T
        np.testing.assert_allclose(restored[:,:3],mesh.vertices,rtol=0,atol=1e-12)
        vectors=np.array([[.001,.002,-.003],[-.02,.005,.01]])
        np.testing.assert_allclose(np.linalg.norm(rotate_vectors(vectors,candidate,root),axis=1),np.linalg.norm(vectors,axis=1),rtol=1e-14,atol=0)
        poses[str(candidate)]=dict(status='PASS',proper_rigid_rotation=True,
            maximum_position_roundtrip_error_mm=float(np.max(np.abs(restored[:,:3]-mesh.vertices))),
            speed_preserved=True,independent_CFD_case=False)
    return dict(status='PASS', closed=True, source_points=len(mesh.vertices),source_triangles=len(mesh.faces),
                mm_to_m_exact=True, ports=['INLET','OUTLET_01','OUTLET_02','OUTLET_03'],poses=poses)

def prepare_surface(output, root=ROOT):
    """Original four-port tagging recipe, writing a new output, never the inputs.

    The accepted source is already a closed printing core. This function does
    not rerun the upstream printing reconstruction or invent a new lumen.
    """
    import trimesh
    import pyvista as pv
    output=fresh_directory(output,root)
    source=root/'inputs/geometry/source_core_mm.stl'
    original=read_json(root/'inputs/geometry/geometry.json')
    assert sha256(source)==original['sha256']
    attachment=read_json(root/'inputs/geometry/port_attachment_summary.json')
    mesh=trimesh.load_mesh(source,process=True)
    assert mesh.is_watertight and mesh.is_winding_consistent and len(mesh.split())==1
    origin=mesh.bounds.mean(axis=0);points=(mesh.vertices-origin)*1e-3
    triangles=np.asarray(mesh.faces);tags=np.ones(len(triangles),np.int32);ports={}
    role_ids={'I1':4,'O1':3,'O2':5,'O3':2}
    for name,tag in role_ids.items():
        value=attachment['per_port'][name];center=np.asarray(value['outside_end_mm']);wall=np.asarray(value['wall_target_mm'])
        normal=(center-wall)/np.linalg.norm(center-wall)
        plane=np.max(np.abs((mesh.vertices[triangles]-center)@normal),axis=1)
        near=np.max(np.linalg.norm(mesh.vertices[triangles]-center,axis=2),axis=1)
        mask=(plane<2e-5)&(near<1.2*value['radius_mm'])
        assert mask.sum()>3 and (tags[mask]==1).all()
        tags[mask]=tag;area=mesh.area_faces[mask].sum()*1e-6
        assert np.min(mesh.face_normals[mask]@normal)>.999
        ports[name]=dict(face_id=tag,role={'I1':'INLET','O1':'OUTLET_01','O2':'OUTLET_02','O3':'OUTLET_03'}[name],
            center_m=((center-origin)*1e-3).tolist(),center_source_mm=center.tolist(),outward_normal=normal.tolist(),
            radius_m=value['radius_mm']*1e-3,area_m2=float(area),triangles=int(mask.sum()))
    surface=pv.PolyData(points,np.c_[np.full(len(triangles),3),triangles].ravel());surface.cell_data['ModelFaceID']=tags
    surface.save(output/'labelled_core_m.vtp')
    np.savez_compressed(output/'source_surface.npz',points_m=points,triangles=triangles,facet_tags=tags)
    record=dict(source=str(source),source_sha256=sha256(source),units={'source':'mm','simulation':'m'},
        origin_source_mm=origin.tolist(),ports=ports,pose_transforms=original['pose_transforms'],
        same_physical_solution_for_candidate_0_and_15=True,CFD_executed=False)
    write_json(output/'GEOMETRY_IDENTITY.json',record)
    return record

