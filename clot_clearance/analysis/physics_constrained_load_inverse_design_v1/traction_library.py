"""Unsteady plane Couette Stokes solution, a prescribed-flow surrogate for streaming.

Fluid occupies 0<y<H: u(0,t)=0, u(H,t)=U(t), u(y,0)=0; dp/dx=0.
u_t=nu*u_yy. The sine expansion solves each diffusion mode by quadrature.
This is NOT a bubble-generated streaming calculation or a global V3 CFD solve.
"""
from common import *
from scipy.integrate import solve_ivp
MU=.001; FLUID_RHO=1000.; GAP=25e-6; RAMP=.1e-3
BOUNDS={'drive_speed_m_s':(.0005,.004),'hold_until_s':(.0004,.0008)}
def physical_admission(params):
 valid=all(k in params and np.isfinite(params[k]) and lo<=params[k]<=hi for k,(lo,hi) in BOUNDS.items())
 return dict(status='VALID_LOCAL_STOKES_APPROXIMATION' if valid else 'OUTSIDE_VALIDATED_RANGE',valid=bool(valid),reason='Both physical parameters must remain in the frozen domain; no traction extrapolation')
def drive(t,duration):
 t=np.asarray(t);q=np.zeros_like(t,dtype=float);dq=np.zeros_like(t,dtype=float)
 a=(t>=0)&(t<RAMP);q[a]=.5*(1-np.cos(np.pi*t[a]/RAMP));dq[a]=.5*np.pi/RAMP*np.sin(np.pi*t[a]/RAMP)
 a=(t>=RAMP)&(t<duration);q[a]=1
 a=(t>=duration)&(t<duration+RAMP);s=(t[a]-duration)/RAMP;q[a]=.5*(1+np.cos(np.pi*s));dq[a]=-.5*np.pi/RAMP*np.sin(np.pi*s)
 return q,dq
class Couette:
 def __init__(self,speed,duration,modes=64):
  self.speed=float(speed);self.duration=float(duration);self.end=duration+.4e-3
  self.re=FLUID_RHO*abs(speed)*GAP/MU
  if not physical_admission(dict(drive_speed_m_s=speed,hold_until_s=duration))['valid']:raise ValueError('OUTSIDE_PHYSICS_DOMAIN')
  self.n=np.arange(1,modes+1);self.lam=(MU/FLUID_RHO)*(self.n*np.pi/GAP)**2
  c=2*(-1.)**(self.n+1)/(self.n*np.pi)
  # Normalized velocity coefficients; pressure is identically zero gauge.
  self.sol=solve_ivp(lambda t,b:-self.lam*b-c*float(drive(t,duration)[1]),(0,self.end),np.zeros(modes),method='BDF',jac=-np.diag(self.lam),rtol=2e-9,atol=1e-12,max_step=5e-6,dense_output=True)
  assert self.sol.success
 def shear(self,t):
  q=drive(t,self.duration)[0];b=self.sol.sol(t)
  return MU*self.speed/GAP*(q+np.einsum('i,i...->...',self.n*np.pi,b))
 def velocity(self,y,t):
  y=np.asarray(y);b=self.sol.sol(t)
  return self.speed*(drive(t,self.duration)[0]*y/GAP+np.sin(np.outer(y/GAP,self.n*np.pi))@b)
 def provenance(self):
  return dict(mechanism='STREAMING_SHEAR',status='PHYSICS_APPROXIMATION',fluid_model='Unsteady plane Stokes Couette; prescribed mean-flow drive, not solved acoustic forcing',mu_Pa_s=MU,rho_kg_m3=FLUID_RHO,gap_m=GAP,drive_speed_m_s=self.speed,Re_H=self.re,Re_admission_max=.1,pressure_increment_Pa=0,drive_onset_s=0,ramp_s=RAMP,hold_until_s=self.duration,drive_off_s=self.duration+RAMP,mechanical_end_s=self.end,diffusion_time_H2_over_nu_s=GAP**2*FLUID_RHO/MU,slowest_mode_time_s=GAP**2*FLUID_RHO/MU/np.pi**2,five_mode_times_s=5*GAP**2*FLUID_RHO/MU/np.pi**2,time_scope='Resolved idealized startup; not Cho experimental streaming formation. No pre-existing steady shear at t=0.',boundary_conditions='stationary bottom no-slip; moving upper plane U(t); zero pressure gradient; infinite lateral planes',v3_scope='Independent tangent-plane patches, projected global +x drive; curvature and cross-patch continuity unverified; not a global V3 fluid solution',coupling='FROZEN_SURFACE_ONE_WAY',increment='Additional prescribed local fluid drive. Historical background CFD is reference-only, not reapplied or rescaled.',limitations=['No bubble source or streaming conversion efficiency','Finite specimen edge effects omitted','No global curved no-slip solution','No deformation feedback or newly exposed fracture faces'])
def check():
 c=Couette(.004,.0008);t=np.linspace(0,c.end,1201);y=np.linspace(0,GAP,51)
 # Independent exact step-start wall shear convolution for a smooth moving plate.
 # Laplace transfer lower-wall shear: sqrt(s/nu)/sinh(H sqrt(s/nu)).
 # An independently integrated finite-difference diffusion model verifies startup.
 N=80;dy=GAP/N;nu=MU/FLUID_RHO
 from scipy.sparse import diags
 mat=diags([np.ones(N-2),-2*np.ones(N-1),np.ones(N-2)],[-1,0,1]).tocsc()*nu/dy**2
 def rhs(t,u):
  z=mat@u;z[-1]+=nu/dy**2*float(drive(t,c.duration)[0]);return z
 fd=solve_ivp(rhs,(0,c.end),np.zeros(N-1),method='BDF',jac=mat,rtol=1e-8,atol=1e-11,t_eval=t,max_step=5e-6)
 tau_fd=MU*c.speed*(4*fd.y[0]-fd.y[1])/(2*dy)
 rel=float(np.max(abs(c.shear(t)-tau_fd))/(MU*c.speed/GAP))
 wall_error=max(float(abs(c.velocity(np.array([0,GAP]),tt)-np.array([0,c.speed*float(drive(tt,c.duration)[0])])).max()) for tt in [0,.0001,.0004,.0009])
 assert rel<.003 and wall_error<1e-12 and abs(float(c.shear(0)))<1e-14
 z=dict(status='PASS',unit='Pa',steady_formula_mu_U_over_H_Pa=MU*c.speed/GAP,finite_difference_vs_series_max_scaled_difference=rel,noslip_error_m_s=wall_error,initial_shear_Pa=float(c.shear(0)),direction='+x at lower wall with n_clot=+y',pressure_Pa=0,convection='identically zero for infinite parallel planes; Re<=0.1 additionally imposed for patch use')
 dump(O/'FLUID_TRACTION_CHECK.json',z);return z
def build(speed=.004,duration=.0008):
 c=Couette(speed,duration);z=load_npz(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz');n=z['n_clot'];direction=np.array([1.,0,0])-n*n[:,0,None]
 t=np.linspace(0,c.end,61);tau=c.shear(t);traction=tau[:,None,None]*direction[None]
 np.savez_compressed(O/'PHYSICAL_TRACTION_LIBRARY.npz',time_s=t,face_id=np.arange(len(n))+1,triangle_mesh_node_ids=z['triangle_mesh_node_ids'],triangles_xyz_m=z['triangles_xyz_m'],face_centroid_m=z['surface_coordinates_m'],area_m2=z['triangle_area_m2'],n_clot=n,pressure_increment_Pa=np.zeros((len(t),len(n))),viscous_traction_increment_Pa=traction,total_traction_increment_Pa=traction,channel_y_m=np.linspace(0,GAP,101),channel_velocity_m_s=np.array([c.velocity(np.linspace(0,GAP,101),tt) for tt in t]),lower_wall_shear_Pa=tau)
 p=c.provenance();p.update(source_surface_sha256=sha(CFD/'refined/CLOT_SURFACE_TRACTIONS.npz'),source_face_id='1-based unchanged row in source surface plus original triangle mesh node IDs',traction_formula='mu*(U/H+sum_n b_n*n*pi/H)*(I-n*n.T)@ex; p_increment=0',background_cfd_status='ACTUAL_SOLVE historical steady 5 mm/s Y flow, deliberately not relabelled streaming')
 dump(O/'TRACTION_PROVENANCE.json',p);check()
if __name__=='__main__':build()
