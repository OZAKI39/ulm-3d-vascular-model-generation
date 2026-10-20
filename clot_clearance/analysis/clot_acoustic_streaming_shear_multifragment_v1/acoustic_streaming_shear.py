"""SI-only flow entrypoint. Never infer streaming from pressure or bubble speed.

Selected fallback is an exact unidirectional periodic Stokes model:
u=(sin(k*y)*v(z,t),0,0), v_t=nu*(v_zz-k^2*v), v(0)=0,
v(H)=U*q(t). Divergence and nonlinear convection vanish identically.
This prescribed upper boundary is not a solved bubble boundary.
"""
from common import *
from scipy.integrate import solve_ivp
MU=.001; FLUID_RHO=1000.; GAP=25e-6; SPEED=.004
WAVELENGTH=700e-6; WAVENUMBER=2*np.pi/WAVELENGTH
RISE=.002; HOLD=.0024; OFF=.0028; COHORT=.0034

def envelope(t):
 t=np.asarray(t);q=np.zeros_like(t,float);dq=q.copy()
 a=(t>=0)&(t<RISE);q[a]=.5*(1-np.cos(np.pi*t[a]/RISE));dq[a]=.5*np.pi/RISE*np.sin(np.pi*t[a]/RISE)
 a=(t>=RISE)&(t<HOLD);q[a]=1
 a=(t>=HOLD)&(t<OFF);z=(t[a]-HOLD)/(OFF-HOLD);q[a]=.5*(1+np.cos(np.pi*z));dq[a]=-.5*np.pi/(OFF-HOLD)*np.sin(np.pi*z)
 return q,dq

def bubble_model_input_gate(inputs):
 required=['frequency_Hz','pressure_field','bubble_positions_m','radius_time_series_m','radius_velocity_m_s','bubble_wall_gap_m','rho_kg_m3','mu_Pa_s','surface_geometry','time_s','first_order_velocity_field','oscillatory_boundary_conditions']
 missing=[k for k in required if inputs.get(k) is None]
 return dict(status='STREAMING_LOAD_NOT_YET_DETERMINED',missing=missing,reason='No boundary-consistent second-order solver is implemented; even complete input is not a solved streaming output.',flow=None)

class StokesMode:
 def __init__(self,k=WAVENUMBER,speed=SPEED,modes=64):
  self.k=k;self.speed=speed;self.n=np.arange(1,modes+1);self.nu=MU/FLUID_RHO
  self.lam=self.nu*((self.n*np.pi/GAP)**2+k*k)
  c=2*(-1.)**(self.n+1)*self.n*np.pi/((self.n*np.pi)**2+(k*GAP)**2)
  self.sol=solve_ivp(lambda t,b:-self.lam*b-c*float(envelope(t)[1]),(0,END),np.zeros(modes),method='BDF',jac=-np.diag(self.lam),rtol=2e-10,atol=1e-13,max_step=5e-6,dense_output=True)
  assert self.sol.success
 def g(self,z):return np.asarray(z)/GAP if self.k==0 else np.sinh(self.k*np.asarray(z))/np.sinh(self.k*GAP)
 def gp(self,z):return np.ones_like(np.asarray(z))/GAP if self.k==0 else self.k*np.cosh(self.k*np.asarray(z))/np.sinh(self.k*GAP)
 def modal_velocity(self,z,t):
  z=np.atleast_1d(z);return self.speed*(self.g(z)*float(envelope(t)[0])+np.sin(np.outer(z/GAP,self.n*np.pi))@self.sol.sol(t))
 def modal_gradient(self,z,t):
  z=np.atleast_1d(z);return self.speed*(self.gp(z)*float(envelope(t)[0])+np.cos(np.outer(z/GAP,self.n*np.pi))@(self.sol.sol(t)*self.n*np.pi/GAP))
 def amplitude(self,t):
  return MU*self.speed*(float(self.gp(0))*envelope(t)[0]+np.einsum('i,i...->...',self.n*np.pi/GAP,self.sol.sol(t)))
 def field(self,y,z,t):
  y,z=np.broadcast_arrays(y,z);v=self.modal_velocity(z.ravel(),t).reshape(z.shape);dv=self.modal_gradient(z.ravel(),t).reshape(z.shape)
  sine=np.sin(self.k*y) if self.k else np.ones_like(y);cosine=np.cos(self.k*y)
  u=np.zeros(y.shape+(3,));grad=np.zeros(y.shape+(3,3));u[...,0]=sine*v
  grad[...,0,1]=self.k*cosine*v;grad[...,0,2]=sine*dv
  return dict(velocity_m_s=u,gradient_s_inv=grad,pressure_increment_Pa=np.zeros(y.shape))

def traction(gradient,pressure,normal):
 sigma=-np.asarray(pressure)[...,None,None]*np.eye(3)+MU*(gradient+np.swapaxes(gradient,-1,-2))
 total=np.einsum('...ij,...j->...i',sigma,normal);normal_t=np.sum(total*normal,axis=-1)[...,None]*normal
 return total-normal_t,normal_t

def spec():
 return dict(model_label='PRESCRIBED_STOKES_SHEAR_TEMPLATE',evidence='PHYSICS_APPROXIMATION',flow_order='PRESCRIBED_NONOSCILLATORY_VISCOUS_FLOW',not_first_order_acoustics=True,not_second_order_streaming=True,
  equation='u_x=sin(k*y)*v(z,t); v_t=nu*(v_zz-k^2*v); p_increment=0',boundary='v(0,t)=0; v(H,t)=0.004*q(t); y-periodic 700 um, infinite x; initial rest',
  g='sinh(k*z)/sinh(k*H)',modal_ODE='b_n_dot=-nu*((n*pi/H)^2+k^2)*b_n-c_n*q_dot; c_n=2*(-1)^(n+1)*n*pi/((n*pi)^2+(k*H)^2)',
  reconstruction='v=U*(q*g+sum b_n*sin(n*pi*z/H))',shear='tau_xz=mu*sin(k*y)*v_z(0,t)',
  rho_kg_m3=FLUID_RHO,mu_Pa_s=MU,gap_m=GAP,upper_speed_m_s=SPEED,wavelength_m=WAVELENGTH,rise_s=RISE,hold_until_s=HOLD,upper_drive_off_s=OFF,end_s=END,cohort_s=COHORT,
  mode_reason='One lowest smooth counter-shear mode over actual 350 um transverse span: opposite x-directed lanes and zero centerline. Selected from specimen dimensions before damage outcomes, not prescribed crack or bond targets.',
  nonlinear_convection='Exactly zero because u=(u_x(y,z,t),0,0), independent of x; this does not validate finite edges or curved surfaces.',
  Re_H=FLUID_RHO*SPEED*GAP/MU,coupling='FROZEN_SURFACE_ONE_WAY',
  limitations=['Prescribed patterned boundary, no physical bubble-to-wall realization','Finite edge effects and moving clot feedback omitted','No new crack-face fluid stress','Reference loading persists on displaced/detached original nodes; only a mechanical dead-traction surrogate once deformation is appreciable','V3 projection is independent local-patch extension, not a global V3 solution'],
  references=[dict(url='https://joss.theoj.org/papers/10.21105/joss.05435',scope='APECSS bubble/acoustic-emission capabilities, not a streaming closure'),dict(url='https://arxiv.org/abs/1509.02554',scope='First/second-order streaming separation and temporal buildup; abstract read, not a reproduction of its full model'),dict(url='https://docs.lammps.org/fix_external.html',scope='Native callback force entrypoint')])

def checks(flow):
 from scipy.sparse import diags
 N=100;dz=GAP/N;nu=MU/FLUID_RHO;k=flow.k;times=np.linspace(0,END,801)
 mat=diags([np.ones(N-2),-2*np.ones(N-1),np.ones(N-2)],[-1,0,1]).tocsc()*nu/dz**2-diags([np.full(N-1,nu*k*k)],[0])
 def rhs(t,u):
  z=mat@u;z[-1]+=nu/dz**2*float(envelope(t)[0]);return z
 sol=solve_ivp(rhs,(0,END),np.zeros(N-1),method='BDF',jac=mat,t_eval=times,rtol=2e-10,atol=1e-13,max_step=5e-6)
 fd=MU*SPEED*(4*sol.y[0]-sol.y[1])/(2*dz);err=float(np.max(abs(fd-flow.amplitude(times)))/(MU*SPEED/GAP))
 wall=max(float(np.max(abs(flow.modal_velocity([0,GAP],t)-[0,SPEED*float(envelope(t)[0])]))) for t in [0,.001,.0024,.0028,.0036])
 steady=MU*SPEED*float(flow.gp(0));near=float(flow.amplitude(.0024));assert err<.001 and wall<1e-13 and abs(near/steady-1)<.002
 f=flow.field(np.array([.000175,-.000175]),np.zeros(2),.0024);ts,tn=traction(f['gradient_s_inv'],f['pressure_increment_Pa'],np.tile([0,0,1.],(2,1)))
 assert np.max(abs(ts[:,0]-np.array([near,-near])))<1e-14 and np.max(abs(tn))==0
 return dict(status='PASS',independent_FD_max_error_over_muU_H=err,no_slip_error_m_s=wall,steady_analytic_mu_U_k_over_sinh_kH_Pa=steady,plateau_shear_Pa=near,plateau_vs_exact_relative=abs(near/steady-1),pure_pressure_shear_N=0,convective_acceleration_exact=0,units='m/s, 1/s, Pa; normals point from solid into fluid')

if __name__=='__main__':
 flow=StokesMode();dump(O/'STREAMING_MODEL_SPEC.json',spec());dump(O/'FLUID_SHEAR_CHECK.json',checks(flow))
 t=np.linspace(0,END,145);y=np.linspace(-WAVELENGTH/2,WAVELENGTH/2,81);z=np.linspace(0,GAP,31);Y,Z=np.meshgrid(y,z,indexing='ij');fields=[flow.field(Y,Z,tt) for tt in t]
 np.savez_compressed(O/'FLOW_FIELD_DATA.npz',time_s=t,y_m=y,z_m=z,velocity_m_s=np.array([f['velocity_m_s'] for f in fields]),velocity_gradient_s_inv=np.array([f['gradient_s_inv'] for f in fields]),pressure_increment_Pa=np.array([f['pressure_increment_Pa'] for f in fields]),wall_shear_amplitude_Pa=flow.amplitude(t))
 print('FLOW_DONE',float(flow.amplitude(t).max()))
