"""Actual top-lattice triangles and conservative barycentric load transfer.

Frozen reference surface: no unearned current-surface/FSI claim.
The high-amplitude case scales traction only; it never rescales a claimed
bubble field or reports an implied high-speed Stokes solution as physical.
"""
from common import *
from acoustic_streaming_shear import *
from scipy.sparse import csr_matrix,save_npz

class ShearLoad:
 def __init__(self,ref,name):
  assert name in ['stokes','inverse'];self.name=name;self.ref=ref;self.flow=StokesMode();self.N=len(ref['X']);X=ref['X'];top=ref['top_ids']-1
  xs=np.unique(X[top,0]);ys=np.unique(X[top,1]);lookup={(round(X[i,0]/DX,6),round(X[i,1]/DX,6)):i for i in top};tri=[]
  for i in range(len(xs)-1):
   for j in range(len(ys)-1):
    q=[lookup[(round(x/DX,6),round(y/DX,6))] for x,y in [(xs[i],ys[j]),(xs[i+1],ys[j]),(xs[i+1],ys[j+1]),(xs[i],ys[j+1])]]
    tri.extend([[q[0],q[1],q[2]],[q[0],q[2],q[3]]] if (i+j)%2==0 else [[q[0],q[1],q[3]],[q[1],q[2],q[3]]])
  self.tri=np.asarray(tri);vertices=X[self.tri];cross=np.cross(vertices[:,1]-vertices[:,0],vertices[:,2]-vertices[:,0]);self.area=.5*np.linalg.norm(cross,axis=1);self.normal=cross/(2*self.area[:,None]);self.cent=vertices.mean(1)
  self.W=csr_matrix((np.full(self.tri.size,1/3),(self.tri.ravel(),np.repeat(np.arange(len(tri)),3))),shape=(self.N,len(tri)))
  self.mode=np.sin(WAVENUMBER*self.cent[:,1]);self.unit_face_force=self.area[:,None]*self.mode[:,None]*np.array([1.,0,0]);self.unit_node_force=self.W@self.unit_face_force
  dense=np.linspace(0,END,14401);self.physical_peak=float(self.flow.amplitude(dense).max()*np.max(abs(self.mode)))
  self.scale=1. if name=='stokes' else 1000/self.physical_peak
  self.boundaries=[RISE,HOLD,OFF,COHORT,END];self.loaded_ids=np.flatnonzero(np.linalg.norm(self.unit_node_force,axis=1)>1e-30)
  self.label='PRESCRIBED_STOKES_SHEAR_TEMPLATE' if name=='stokes' else 'INVERSE_SHEAR_REQUIREMENT'
 def amplitude(self,t):return self.scale*self.flow.amplitude(t)
 def force(self,t):return float(self.amplitude(t))*self.unit_node_force
 def faces(self,t):return float(self.amplitude(t))*self.unit_face_force
 def peak(self,t):return float(abs(self.amplitude(t))*np.max(abs(self.mode)))
 def phase(self,t):return 'RISE' if t<RISE else 'HOLD' if t<HOLD else 'DRIVE_DOWN' if t<OFF else 'VISCOUS_DECAY_OBSERVATION'
 def mean_force(self,t0,t1,order=4):
  q,w=np.polynomial.legendre.leggauss(order);t=(t1+t0)/2+(t1-t0)*q/2
  return .5*float(w@self.amplitude(t))*self.unit_node_force
 def work(self,x0,x1,t0,t1,order=4):return float(np.sum(self.mean_force(t0,t1,order)*(x1-x0)))
 def geometry_metrics(self,x):
  moved=x[self.tri];cross=np.cross(moved[:,1]-moved[:,0],moved[:,2]-moved[:,0]);area=.5*np.linalg.norm(cross,axis=1);normal=cross/np.maximum(2*area[:,None],1e-30)
  return dict(max_loaded_displacement_m=float(np.linalg.norm(x[self.loaded_ids]-self.ref['X'][self.loaded_ids],axis=1).max()),max_surface_normal_rotation_deg=float(np.degrees(np.arccos(np.clip(np.sum(normal*self.normal,axis=1),-1,1))).max()),max_surface_area_ratio=float((area/self.area).max()),min_surface_area_ratio=float((area/self.area).min()))
 def audit(self):
  X=self.ref['X'];ff=self.faces(HOLD);fn=self.force(HOLD);v=np.random.default_rng(20261017).normal(size=X.shape)*.001;vf=self.W.T@v
  force_error=float(np.linalg.norm(fn.sum(0)-ff.sum(0)));moment_error=float(np.linalg.norm(np.cross(X,fn).sum(0)-np.cross(self.cent,ff).sum(0)));power_error=float(abs(np.sum(fn*v)-np.sum(ff*vf)))
  assert force_error<1e-18 and moment_error<1e-21 and power_error<1e-20 and np.all(self.normal[:,2]>1-1e-14)
  return dict(status='PASS',triangles=len(self.tri),loaded_material_nodes=len(self.loaded_ids),area_m2=float(self.area.sum()),min_weight=1/3,max_weight=1/3,negative_weights=0,force_error_N=force_error,moment_error_Nm=moment_error,virtual_power_error_W=power_error,max_node_force_N=float(np.linalg.norm(fn,axis=1).max()),sum_nodal_force_norm_N=float(np.linalg.norm(fn,axis=1).sum()),sum_face_force_norm_N=float(np.linalg.norm(ff,axis=1).sum()),net_force_N=fn.sum(0),torque_Nm=np.cross(X,fn).sum(0),local_force_amplification_bound='Convex 1/3 weights; no extrapolation or sign reversal introduced by mapping',coupling='FROZEN_SURFACE_ONE_WAY',surface_definition='Top boundary lattice squares split into actual original-node triangles, alternating diagonals; not an invented medical shell',force_quadrature='Piecewise-constant centroid stress per actual triangle; barycentric 1/3 distribution',normal_load_N=0,physical_baseline_peak_face_Pa=self.physical_peak,traction_multiplier=self.scale)

def freeze():
 ref=load_npz(O/'provenance/reference.npz');a=ShearLoad(ref,'stokes');b=ShearLoad(ref,'inverse');times=np.linspace(0,END,145)
 save_npz(O/'SMALL_SURFACE_MAPPING.npz',a.W)
 np.savez_compressed(O/'SHEAR_TRACTION_FIELD.npz',time_s=times,original_atom_ids=np.arange(1,len(ref['X'])+1),loaded_atom_ids=a.loaded_ids+1,reference_coordinates_m=ref['X'],triangle_atom_ids=a.tri+1,triangle_vertices_m=ref['X'][a.tri],centroid_m=a.cent,area_m2=a.area,n_clot=a.normal,stokes_shear_traction_Pa=np.array([a.faces(t)/a.area[:,None] for t in times]),normal_traction_Pa=np.zeros((len(times),len(a.tri),3)),stokes_node_force_N=np.array([a.force(t) for t in times]),inverse_traction_multiplier=np.array(b.scale))
 dump(O/'SHEAR_PROTOCOL.json',dict(cases={'stokes':dict(label=a.label,peak_face_shear_Pa=a.physical_peak,upper_speed_m_s=SPEED),'inverse':dict(label=b.label,evidence='MECHANICAL_REQUIREMENT_ONLY',peak_face_shear_Pa=1000,traction_multiplier=b.scale,gate='Only after A completed with no meaningful fragmentation; known historical 1000 Pa uniform shear detached 882 nodes, but counter-shear and unloading differ so target threshold is unknown.')},same_geometry_and_space_time_mode=True,max_formal_cases=2,native_base_dt_s=BASE,min_dt_s=MIN_H,max_dd=.1,relative_energy_tolerance=.01,begin_s=0,end_s=END,post_upper_drive_stop_s=OFF,near_unloaded_cohort_s=COHORT,cohort_remaining_load_fraction=float(a.peak(COHORT)/a.physical_peak),comparison='B/A is a controlled multiplier between these templates, NOT a bubble capacity deficit or success threshold',surface_limitation='Frozen top nodes remain loaded after motion or detachment. No updated fluid field; flag first displacement>0.1gap or rotation>5deg as conservative diagnostic, not a validated cutoff.'))
 dump(O/'SMALL_SURFACE_MAPPING_AUDIT.json',{'stokes':a.audit(),'inverse':b.audit()});print('SHEAR_PROTOCOL_FROZEN',a.physical_peak,b.scale)
if __name__=='__main__':freeze()
