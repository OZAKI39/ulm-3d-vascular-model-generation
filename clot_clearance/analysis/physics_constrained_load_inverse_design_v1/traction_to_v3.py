"""Scalar affine-consistent weights, adapted from accepted local surface support.

Each column sums to one and reproduces its face centroid. This guarantees
force and moment for constant face tractions; the transpose defines power.
Signed extrapolation weights are required where nodes lie below the surface.
"""
from common import *
from scipy.sparse import csr_matrix,save_npz
from scipy.spatial import cKDTree
class FrozenV3Provider:
 """Callable SI force port for original V3 material IDs, independent of x(t)."""
 def __init__(self,params):
  from scipy.sparse import load_npz as load_sparse
  from traction_library import Couette
  self.flow=Couette(params['drive_speed_m_s'],params['hold_until_s']);self.W=load_sparse(O/'AFFINE_SCALAR_MAPPING.npz');self.surface=load_npz(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz');n=self.surface['n_clot']
  self.unit_force=self.W@(self.surface['triangle_area_m2'][:,None]*(np.array([1.,0,0])-n*n[:,0,None]));self.ids=load_npz(MAP/'GEOMETRY_PREFLIGHT.npz')['material_ids']
 def force(self,time_s,original_atom_ids=None):
  force=float(self.flow.shear(time_s))*self.unit_force
  if original_atom_ids is None:return force
  ids=np.asarray(original_atom_ids,dtype=int);out=np.zeros((len(ids),3));material=(ids>=1)&(ids<=len(self.ids));out[material]=force[ids[material]-1];return out
 def replace_load(self,existing_external_state,time_s):
  """Keep the accepted V3 wall/contact/background/numeric ports unchanged.

  This is an offline-checked insertion interface; no new full-V3 run is
  admitted while the local patch fluid geometry and fracture response fail.
  """
  z={k:v.copy() for k,v in existing_external_state.items()};new=self.force(time_s);z['external_force']=z['external_force']-z['load']+new;z['load']=new;return z
def high_fidelity_gate(summaries):
 localized=any(z['status']=='COMPLETED' and z['last_row']['new_softened']>0 and z['last_row']['boundary_damage']>z['last_row']['internal_damage'] for z in summaries)
 z=dict(status='NOT_RUN_EVIDENCE_INSUFFICIENT',new_full_V3_runs=0,small_localized_damage=localized,curved_V3_fluid_solution_validated=False,existing_solver=str(R/'analysis/v3_chunk_load_direction_control_v1/run_protocol.py'),mechanical_port='FrozenV3Provider.force / replace_load; original force units and material IDs',reason='No demonstrated localized cohesive damage; independent plane-patch flow is not a global curved no-slip solution',next_required='A validated spatially localized near-wall flow source on actual geometry, then repeat the small-mechanism admission before full V3',fidelity_separation='No small-specimen response is entered as a V3 observation')
 dump(O/'FULL_V3_ADMISSION.json',z);return z
def build():
 old=load_npz(MAP/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz');geom=load_npz(MAP/'GEOMETRY_PREFLIGHT.npz');z=load_npz(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz')
 W=csr_matrix((old['data'],old['indices'],old['indptr']),shape=tuple(old['shape'])).tocsc();X=geom['X'];cent=z['surface_coordinates_m'];conditions=[];errors=[];rr=[];cc=[];vv=[];expanded=[];tree=cKDTree(X)
 for f in range(W.shape[1]):
  sl=slice(W.indptr[f],W.indptr[f+1]);ids=W.indices[sl];w=W.data[sl]
  C=np.vstack([np.ones(len(ids)),((X[ids]-cent[f])/HORIZON).T]);target=np.array([1.,0,0,0])
  if np.linalg.matrix_rank(C)<4 or np.linalg.cond(C@C.T)>1e7:
   # Surface-only supports can be coplanar: augment only those with nearby
   # original material nodes, never anchors or invented quadrature particles.
   distance,near=tree.query(cent[f],k=24);newids=np.union1d(ids,near);neww=np.zeros(len(newids));neww[np.searchsorted(newids,ids)]=w;ids=newids;w=neww
   C=np.vstack([np.ones(len(ids)),((X[ids]-cent[f])/HORIZON).T]);expanded.append(dict(face_id=f+1,node_count=len(ids),max_distance_m=float(np.linalg.norm(X[ids]-cent[f],axis=1).max())))
  D=np.maximum(w,.05/len(ids));Q=(C*D)@C.T;conditions.append(np.linalg.cond(Q));dw=D*(C.T@np.linalg.solve(Q,target-C@w));v=w+dw;errors.append(np.max(abs(C@v-target)));rr.extend(ids);cc.extend([f]*len(ids));vv.extend(v)
 W=csr_matrix((vv,(rr,cc)),shape=W.shape);save_npz(O/'AFFINE_SCALAR_MAPPING.npz',W);dump(O/'MAPPING_SUPPORT_AUGMENTATION.json',dict(reason='Initial old surface-only support had singular affine constraints; minimal nearby original-material augmentation, prior mapping untouched',faces=expanded))
 lib=load_npz(O/'PHYSICAL_TRACTION_LIBRARY.npz');tf=lib['total_traction_increment_Pa'][np.argmax(lib['lower_wall_shear_Pa'])];ff=tf*z['triangle_area_m2'][:,None];fn=W@ff;x0=z['reference_point_m']
 vel=np.random.default_rng(1016).normal(size=X.shape)*.01;vf=W.T@vel
 Ff=ff.sum(0);Fn=fn.sum(0);Mf=np.cross(cent-x0,ff).sum(0);Mn=np.cross(X-x0,fn).sum(0);Pf=float(np.sum(ff*vf));Pn=float(np.sum(fn*vel))
 errF=float(np.linalg.norm(Ff-Fn)/max(np.linalg.norm(Ff),1e-30));errM=float(np.linalg.norm(Mf-Mn)/max(np.linalg.norm(Mf),1e-30));errP=abs(Pf-Pn)/max(abs(Pf),1e-30)
 assert max(errors)<1e-10 and max(errF,errM,errP)<1e-10
 # A rigid velocity field is reproduced exactly at the reference face centroids.
 omega=np.array([3.,-2.,1.]);rigid=np.array([.002,.001,-.003])+np.cross(omega,X-x0);rigiderr=float(np.max(abs(W.T@rigid-(np.array([.002,.001,-.003])+np.cross(omega,cent-x0)))))
 np.savez_compressed(O/'V3_MAPPED_FORCES.npz',material_atom_id=geom['material_ids'],reference_coordinates_m=X,force_N=fn,face_force_N=ff,reference_point_m=x0,test_node_velocity_m_s=vel,test_face_velocity_m_s=vf)
 dump(O/'FORCE_MAPPING_AUDIT.json',dict(status='PASS',coupling='FROZEN_SURFACE_ONE_WAY',force_relative_error=errF,moment_relative_error=errM,power_relative_error=errP,face_power_W=Pf,node_power_W=Pn,rigid_velocity_max_error_m_s=rigiderr,max_affine_constraint_error=max(errors),max_condition_number=max(conditions),negative_weight_fraction=float(np.mean(W.data<0)),min_weight=float(W.data.min()),max_weight=float(W.data.max()),max_column_absolute_sum=float(abs(W).sum(0).max()),support='Accepted historical local/geodesic support; per-column weighted affine correction',minimum_interface_change='Scalar weights now reproduce reference centroid, instead of load-dependent global force/moment correction; same weights transpose interpolate velocity',power_scope='Exact algebraic virtual power in reference geometry, not physical moving-boundary power or dynamic FSI',units='Pa*m2=N; m cross N=Nm; N dot m/s=W',normal_sign='n_clot points into fluid; positive upper drive gives positive lower-wall tangential force',source_hashes={p.name:sha(p) for p in [MAP/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz',CFD/'refined/CLOT_SURFACE_TRACTIONS.npz']}))
if __name__=='__main__':build()
