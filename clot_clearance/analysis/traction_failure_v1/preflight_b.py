from analysis.traction_failure_v1.common import *
from scipy.spatial import cKDTree
from analysis.dual_track_v1.measure_b import Measurements

def main():
 t=time.perf_counter();cpu=time.process_time();out=B/'GEOMETRY_PREFLIGHT.json';assert not out.exists()
 z=dict(np.load(V3/'V3_MATERIAL_POINTS_AND_BONDS.npz'));n=int(z['clot_count']);X=z['X'][:n];s=dict(np.load(PREV/'refined/CLOT_SURFACE_TRACTIONS.npz'))
 rr=json.loads((PREV/'refined/CFD_RESULT.json').read_text());a=np.load(PREV/'refined/SV_MESH/mesh_arrays.npz');P=a['points_m'];T=a['tetra'];faces=a['boundary_triangles'];tag=a['facet_tags'];sel=tag==5;f=faces[sel];adj=a['adjacent_tetra'][sel]
 mm=Measurements(PREV/'refined/SV_MESH/mesh_arrays.npz',2.6354060534101725e-8,.005);u,p=mm.read(ROOT/rr['snapshot']);xyz=P[f];cent=xyz.mean(1);av=.5*np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);area=np.linalg.norm(av,axis=1);nf=av/area[:,None];nc=-nf
 assert np.all(np.einsum('ij,ij->i',nf,cent-P[T[adj]].mean(1))>0)
 grad=np.linalg.solve(P[T[adj,1:]]-P[T[adj,0]][:,None],u[T[adj,1:]]-u[T[adj,0]][:,None]).transpose(0,2,1)
 visc=.001*np.einsum('nij,nj->ni',grad+grad.transpose(0,2,1),nc);pressure=-p[f].mean(1)[:,None]*nc;vn=np.einsum('ij,ij->i',visc,nc)[:,None]*nc;shear=visc-vn
 compare={'triangles_xyz_m':xyz,'triangle_mesh_node_ids':f+1,'adjacent_mesh_tetra_ids':adj+1,'triangle_area_m2':area,'n_clot':nc,'n_fluid':nf,'pressure_Pa':p[f].mean(1),'velocity_gradient_s_inv':grad,'pressure_traction_Pa':pressure,'viscous_normal_traction_Pa':vn,'tangential_traction_Pa':shear,'total_traction_Pa':pressure+visc}
 errors={k:float(np.max(abs(s[k]-v))) for k,v in compare.items()};assert max(errors.values())<1e-12
 assert sha(ROOT/rr['snapshot'])==rr['snapshot_sha256'];assert sha(PREV/'refined/SV_MESH/mesh-complete.mesh.vtu')==rr['mesh_sha256']
 body=V3/'cases/C1_GAUSSIAN_A200/body.data';text=body.read_text();atoms=text.split('Atoms')[1].split('Bonds')[0].splitlines();rows=[]
 for line in atoms:
  v=line.split()
  if len(v)==6:
   try:rows.append([int(v[0]),*map(float,v[3:])])
   except ValueError:pass
 rows=np.array(sorted(rows));assert len(rows)==len(z['X']);assert np.array_equal(rows[:,0],np.arange(1,len(rows)+1));assert np.max(abs(rows[:,1:]-z['X']))<1e-17
 near,tri,d=closest_surface(X,xyz);face_dist,face_node=cKDTree(X).query(cent);quad_dist,_=cKDTree(X).query(s['quadrature_coordinates_m'].reshape(-1,3))
 geom=np.load(PREV/'geometry_arrays.npz');expected=geom['points_m'][geom['triangles'][geom['tags']==5]]
 # Coordinate multiset equality independently excludes vessel wall / attachment faces.
 def keys(x):return {tuple(sorted(tuple(y) for y in tri)) for tri in x}
 assert keys(expected)==keys(xyz)
 inputpaths=[body,V3/'V3_MATERIAL_POINTS_AND_BONDS.npz',V3/'V3_GEOMETRY_CONFIG.json',PREV/'refined/CLOT_SURFACE_TRACTIONS.npz',PREV/'geometry_arrays.npz',PREV/'refined/SV_MESH/mesh_arrays.npz',ROOT/rr['snapshot'],PREV/'refined/CFD_CLOT_TRACTION_AUDIT.json',PREV/'refined/CFD_MASS_BALANCE.json',PREV/'MESH_SENSITIVITY.json']
 save(B/'INPUT_SHA.json',{str(p.relative_to(ROOT)):sha(p) for p in inputpaths})
 np.savez_compressed(B/'GEOMETRY_PREFLIGHT.npz',material_ids=rows[:n,0].astype(np.int64),X=X,anchors=z['X'][n:],closest_surface_m=near,closest_surface_triangle=tri,distance_to_surface_m=d,face_nearest_material_distance_m=face_dist,face_nearest_material_index=face_node,original_surface=z['surface'][:n],interface_material_indices=np.unique(z['bonds'][z['bond_types']==2][z['bonds'][z['bond_types']==2]<n]))
 save(B/'NATIVE_SURFACE_IDENTITY.json',dict(snapshot=rr['snapshot'],checks=errors,all_original_faces_exact=True,original_field_roundtrip=mm.identity,no_attachment_faces=True,normal_verified_by_adjacent_tet=True,pressure_reference='Original zero-reference outlet traction; relative gauge pressure',original_surface_sha256=sha(PREV/'refined/CLOT_SURFACE_TRACTIONS.npz')))
 save(out,dict(material_nodes=n,anchors=int(z['anchors'].sum()),original_surface_nodes=int(z['surface'].sum()),material_distance_m=stats(d),face_nearest_material_distance_m=stats(face_dist),quadrature_nearest_material_distance_m=stats(quad_dist),within_one_dx=int((d<=DX).sum()),within_sqrt3_half_dx=int((d<=3**.5/2*DX).sum()),body_coordinate_error_m=float(abs(rows[:,1:]-z['X']).max()),runtime=cost(t,cpu)))
 print(out.read_text())
if __name__=='__main__':main()
