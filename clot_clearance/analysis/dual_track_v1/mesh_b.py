from analysis.dual_track_v1.common import *
import os,sys,subprocess,time,pyvista as pv
from scipy.spatial import cKDTree
def compact(points, triangles):
    used, inverse = np.unique(triangles, return_inverse=True)
    face = pv.PolyData(points[used], np.column_stack([np.full(len(triangles),3), inverse.reshape(-1,3)]))
    face.point_data["GlobalNodeID"] = (used+1).astype(np.int32)
    return face


def audit(raw, output, workflow):
    sys.path.insert(0, str(workflow))
    from vascular_validation.mesh_diagnostics import volume_validity, tetra_sicn, quality_summary
    grid = pv.read(raw / "volume.vtu")
    if set(grid.cells_dict) != {pv.CellType.TETRA}:
        raise ValueError("Expected linear tetrahedra from the existing mesher")
    points = np.array(grid.points, dtype=float)
    tetra = np.asarray(grid.cells_dict[pv.CellType.TETRA], np.int64)
    validity = volume_validity(points, tetra)
    quality = tetra_sicn(points, tetra)
    stats = quality_summary(quality)
    if not stats["q_min"] > 0:
        raise ValueError("Nonpositive signed tetrahedron quality")
    all_faces = np.concatenate([tetra[:, indices] for indices in ([0,1,2],[0,1,3],[0,2,3],[1,2,3])])
    owners = np.tile(np.arange(len(tetra)),4)
    _, first, counts = np.unique(np.sort(all_faces,axis=1),axis=0,return_index=True,return_counts=True)
    assert np.max(counts) <= 2
    boundary = all_faces[first[counts == 1]].copy()
    adjacent = owners[first[counts == 1]]
    xyz = points[boundary]
    direction = xyz.mean(axis=1)-points[tetra[adjacent]].mean(axis=1)
    wrong = np.einsum("ij,ij->i",np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),direction) < 0
    boundary[wrong] = boundary[wrong][:,[0,2,1]]
    lookup = {tuple(sorted(face)):i for i,face in enumerate(boundary)}
    tags = np.zeros(len(boundary),np.int32)
    tree = cKDTree(points)
    for face_id in (1,2,3,4,5):
        face = pv.read(raw / ("face_"+str(face_id)+".vtp"))
        distances, index = tree.query(face.points)
        assert distances.max() < 1e-10
        for triangle in index[face.faces.reshape(-1,4)[:,1:]]:
            j = lookup[tuple(sorted(triangle))]
            assert tags[j] == 0
            tags[j] = face_id
    assert np.all(tags > 0)
    # SV's vtkPoints export is float32. Restore only the unchanged PLC nodes to
    # original double coordinates after proving a bijective, exact float32 match.
    plc=pv.read(B/'clot_present_lumen.vtp');used=np.unique(boundary);distance,index=cKDTree(plc.points).query(points[used])
    assert len(used)==plc.n_points and len(np.unique(index))==plc.n_points
    assert np.array_equal(points[used].astype(np.float32),plc.points[index].astype(np.float32))
    restore=dict(native_coordinate_dtype=str(grid.points.dtype),boundary_nodes=len(used),max_roundoff_restore_m=float(distance.max()),bijection=True,exact_float32_roundtrip=True)
    points[used]=plc.points[index];grid.points=points
    validity=volume_validity(points,tetra);quality=tetra_sicn(points,tetra);stats=quality_summary(quality)
    assert stats['q_min']>0
    (output / "mesh-surfaces").mkdir(parents=True)
    grid.point_data["GlobalNodeID"] = np.arange(1,len(points)+1,dtype=np.int32)
    grid.cell_data["GlobalElementID"] = np.arange(1,len(tetra)+1,dtype=np.int32)
    grid.cell_data["minSICN"] = quality
    grid.save(output / "mesh-complete.mesh.vtu")
    surface = compact(points,boundary)
    surface.cell_data["ModelFaceID"] = tags
    surface.cell_data["GlobalElementID"] = (adjacent+1).astype(np.int32)
    surface.cell_data["BranchID"] = np.where(surface.cell_centers().points[:,0] <= 0,0,
                                            np.where(surface.cell_centers().points[:,1] >= 0,1,2))
    surface.save(output / "mesh-complete.exterior.vtp")
    assert surface.n_open_edges == 0
    face_info = {}
    for name,face_id in [("WALL",1),("INLET",2),("OUTLET_A",3),("OUTLET_B",4),("CLOT",5)]:
        selected = tags == face_id
        face = compact(points,boundary[selected])
        face.cell_data["GlobalElementID"] = (adjacent[selected]+1).astype(np.int32)
        face.save(output / "mesh-surfaces" / (name+".vtp"))
        face_info[name] = dict(triangles=face.n_cells,area_m2=float(face.area))
    np.savez_compressed(output / "mesh_arrays.npz", points_m=points,tetra=tetra,
                        boundary_triangles=boundary,facet_tags=tags,adjacent_tetra=adjacent)
    return dict(precision_restore=restore,validity=validity,quality=stats,faces=face_info,points=len(points),tetrahedra=len(tetra),
                diagnostic_source=str(workflow / "vascular_validation/mesh_diagnostics.py"),
                diagnostic_source_sha256=sha(workflow / "vascular_validation/mesh_diagnostics.py"))


def main(level):
 design=json.loads((B/'DESIGN_FREEZE.json').read_text());m=design['mesh'];out=B/level;out.mkdir(exist_ok=False)
 if level=='refined':assert json.loads((B/'baseline/CFD_RESULT.json').read_text())['status']=='PASS_WITH_SCOPE'
 source=B/'clot_present_lumen_int32.vtp';representation=json.loads((B/'GEOMETRY_REPRESENTATION_ADDENDUM.json').read_text());assert sha(source)==representation['new_sha256'];assert representation['old_sha256']==design['geometry_sha256']
 req=dict(surface=str(source),output=str(out/'raw'),global_edge_m=m['global_edge_m'],local_edge_m=m[level+'_local_edge_m'],sphere_radius_m=m['local_sphere_radius_m'],sphere_center_m=m['local_sphere_center_m']);save(out/'MESH_REQUEST.json',req)
 env=dict(os.environ,DUAL_TRACK_MESH_REQUEST=str(out/'MESH_REQUEST.json'),PYTHONNOUSERSITE='1',QT_QPA_PLATFORM='offscreen',XDG_CONFIG_HOME=str(out/'app_config'),XDG_CACHE_HOME=str(out/'app_cache'));env.pop('PYTHONPATH',None)
 t=time.monotonic();cmd=['/workspace/y_clot/bin/y-clot-simvascular','-python','--',str(ROOT/'analysis/dual_track_v1/mesh_embedded.py')]
 with (out/'mesh.log').open('w') as log:r=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=m['budgets']['per_mesh_wall_s'])
 save(out/'MESH_EXECUTION.json',dict(returncode=r.returncode,wall_s=time.monotonic()-t,command=cmd));assert r.returncode==0
 report=audit(out/'raw',out/'SV_MESH',Path('/workspace/y_clot/deps/y_cfd_vendor'))
 geom=json.loads((B/'GEOMETRY_AUDIT.json').read_text());grid=pv.read(out/'SV_MESH/mesh-complete.mesh.vtu');volume=float(grid.volume);relative=abs(volume/geom['fluid_volume_m3']-1)
 report.update(volume_m3=volume,relative_volume_PLC=relative,PLC_sha256=design['geometry_sha256'],local_target_m=req['local_edge_m'],mesh_sha256=sha(out/'SV_MESH/mesh-complete.mesh.vtu'))
 # Every generated boundary point remains on the frozen PLC, without remeshing.
 surface=pv.read(out/'SV_MESH/mesh-complete.exterior.vtp');plc=pv.read(B/'clot_present_lumen.vtp');dist,_=cKDTree(plc.points).query(surface.points);report['max_boundary_node_distance_from_PLC_m']=float(dist.max())
 checks=dict(positive_quality=report['quality']['q_min']>=m['quality']['min_SICN'],volume=relative<=m['quality']['relative_volume_to_PLC'],points=report['points']<=m['budgets']['max_points'],elements=report['tetrahedra']<=m['budgets']['max_tetrahedra'],closed=surface.n_open_edges==0,unchanged_boundary=dist.max()<1e-10)
 if level=='refined':
  coarse=np.load(B/'baseline/SV_MESH/mesh_arrays.npz');fine=np.load(out/'SV_MESH/mesh_arrays.npz');local=[]
  for data in [coarse,fine]:
   tet=data['tetra'];pp=data['points_m'];xyz=pp[tet];inside=np.linalg.norm(xyz.mean(1)-m['local_sphere_center_m'],axis=1)<m['local_sphere_radius_m'];vv=np.linalg.det(xyz[:,1:]-xyz[:,:1])/6;local.append(dict(local_tetrahedra=int(inside.sum()),local_median_cell_volume_m3=float(np.median(vv[inside]))))
  checks['actual_local_refinement']=local[1]['local_tetrahedra']>local[0]['local_tetrahedra'] and local[1]['local_median_cell_volume_m3']<local[0]['local_median_cell_volume_m3'];report['local_refinement_measured']=local
 report.update(checks=checks,status='PASS' if all(checks.values()) else 'FAIL',wall_s=time.monotonic()-t);save(out/'MESH_QUALITY_AUDIT.json',report);save(B/'MESH_QUALITY_AUDIT.json',dict(level=level,**report));assert all(checks.values()),checks
if __name__=='__main__':main(sys.argv[1])
