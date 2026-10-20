"""Compare existing signed affine map with a separate nonnegative local kernel.
Neither clipping signed weights nor claiming all force-conserving maps are local.
No full V3 dynamics is run.
"""
from common import *
from acoustic_streaming_shear import *
from scipy.sparse import load_npz as sparse_load,csr_matrix,save_npz
from scipy.spatial import cKDTree

def main():
 z=load_npz(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz');geom=load_npz(MAP/'GEOMETRY_PREFLIGHT.npz');X=geom['X'];n=z['n_clot'];cent=z['surface_coordinates_m'];area=z['triangle_area_m2'];W=sparse_load(INV/'AFFINE_SCALAR_MAPPING.npz')
 # Local positive alternative is independently constructed from geometry,
 # NOT a zero-clipped signed affine map. It preserves partition of unity only.
 dist,ids=cKDTree(X).query(cent,k=24);w=np.exp(-.5*(dist/HORIZON)**2)/(dist+DX/4)**2;w/=w.sum(1)[:,None]
 Wp=csr_matrix((w.ravel(),(ids.ravel(),np.repeat(np.arange(len(cent)),24))),shape=W.shape)
 flow=StokesMode();tau=float(flow.amplitude(HOLD));yc=.5*(X[:,1].min()+X[:,1].max());sine=np.sin(WAVENUMBER*(cent[:,1]-yc));direction=np.array([1.,0,0])-n*n[:,0,None]
 ts=tau*sine[:,None]*direction;ff=area[:,None]*ts;v=np.random.default_rng(43).normal(size=X.shape)*.001;records={};forces={};moment=np.cross(cent,ff).sum(0);total=float(np.linalg.norm(ff,axis=1).sum());mscale=float(np.sum(np.linalg.norm(np.cross(cent-X.mean(0),ff),axis=1)))
 for name,M in [('signed_affine',W),('nonnegative_local_kernel',Wp)]:
  fn=M@ff;forces[name]=fn;vn=M.T@v;moment_node=np.cross(X,fn).sum(0);reconstructed=M.T@X
  neighborhood=cKDTree(X).query(X,k=7)[1][:,1:];arean=M@area;density=fn/np.maximum(abs(arean[:,None]),DX**2/100);rough=float(np.mean(np.linalg.norm(density[:,None]-density[neighborhood],axis=2)))
  record=dict(force_abs_error_N=float(np.linalg.norm(fn.sum(0)-ff.sum(0))),force_error_over_total_abs_force=float(np.linalg.norm(fn.sum(0)-ff.sum(0))/total),moment_abs_error_Nm=float(np.linalg.norm(moment_node-moment)),moment_error_over_total_abs_moment=float(np.linalg.norm(moment_node-moment)/mscale),power_abs_error_W=float(abs(np.sum(fn*v)-np.sum(ff*vn))),min_weight=float(M.data.min()),max_weight=float(M.data.max()),negative_weight_fraction=float(np.mean(M.data<0)),max_column_abs_sum=float(abs(M).sum(0).max()),max_centroid_reproduction_error_m=float(np.linalg.norm(reconstructed-cent,axis=1).max()),sum_node_force_norm_over_sum_face_norm=float(np.linalg.norm(fn,axis=1).sum()/total),max_nodal_force_N=float(np.linalg.norm(fn,axis=1).max()),nearest_neighbor_force_per_area_variation_Pa=rough)
  assert record['force_error_over_total_abs_force']<1e-10 and record['power_abs_error_W']<1e-20
  if name=='signed_affine':assert record['moment_error_over_total_abs_moment']<1e-10
  records[name]=record
 f=forces['signed_affine'];p=forces['nonnegative_local_kernel'];good=(np.linalg.norm(f,axis=1)>1e-18)&(np.linalg.norm(p,axis=1)>1e-18);dot=np.sum(f*p,axis=1);angles=np.degrees(np.arccos(np.clip(dot[good]/(np.linalg.norm(f[good],axis=1)*np.linalg.norm(p[good],axis=1)),-1,1)))
 signed_absolute_support=abs(W)@np.linalg.norm(ff,axis=1);positive_budget=Wp@np.linalg.norm(ff,axis=1);active=positive_budget>1e-18;ratio=np.linalg.norm(f[active],axis=1)/positive_budget[active]
 records['local_comparison']=dict(common_active_nodes=int(good.sum()),nodes_direction_difference_over_90deg=int(np.sum(angles>90)),direction_difference_deg_quantiles=np.quantile(angles,[0,.5,.95,1]),signed_force_over_positive_local_force_budget_quantiles=np.quantile(ratio,[0,.5,.95,1]),no_negative_weight_clipping=True,warning='Direction differences are aggregate actual nodal forces for this field, not inferred from weight signs. A positive map sacrifices exact reference moment; signed map may amplify local loads despite exact global balance.')
 records['decision']=dict(status='V3_MAPPING_NOT_ADMITTED_FOR_DYNAMIC_STREAMING',selected_small_specimen='Actual-node top triangles use nonnegative barycentric weights with exact force/moment and adjoint virtual power',V3='Signed map retained for diagnostics, positive kernel is comparison only. Need boundary-conforming material/support refinement or a validated regularized traction transfer before claiming local V3 accuracy.',scope='Actual CFD exposed triangles; imposed local tangent-plane Stokes stress, not computed curved V3 streaming')
 records['source_sha256']={str(p):sha(p) for p in [CFD/'refined/CLOT_SURFACE_TRACTIONS.npz',MAP/'GEOMETRY_PREFLIGHT.npz',INV/'AFFINE_SCALAR_MAPPING.npz']}
 save_npz(O/'V3_NONNEGATIVE_COMPARISON_MAPPING.npz',Wp)
 np.savez_compressed(O/'V3_SURFACE_MAPPING_DATA.npz',material_atom_id=geom['material_ids'],X=X,triangles_xyz_m=z['triangles_xyz_m'],triangle_mesh_node_ids=z['triangle_mesh_node_ids'],centroid_m=cent,area_m2=area,n_clot=n,shear_traction_Pa=ts,normal_traction_Pa=np.zeros_like(ts),face_force_N=ff,signed_node_force_N=f,nonnegative_node_force_N=p,positive_local_force_budget_N=positive_budget)
 dump(O/'SURFACE_MAPPING_AUDIT.json',dict(v3=records,small=json.loads((O/'SMALL_SURFACE_MAPPING_AUDIT.json').read_text())))
 print('V3_MAPPING_AUDIT_DONE',len(cent),records['local_comparison'])
if __name__=='__main__':main()
