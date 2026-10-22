"""Direct 3D F1 equilibria and drift Jacobians, never particle-pattern claims."""
import time
import numpy as np,h5py
from scipy.optimize import brentq
from mevo_mb.acoustics.water_disc_proxy import WaterDiscs,CachedWaterBasis
from scripts.straight_tube_exp1.f10af.field import response_context,direct_force
from .common import *
def force_context(f,deg=0):
 c=response_context(f);model=WaterDiscs(c['discs'].elements,f);a=[0,0,1000,1000];p=[0,0,0,np.radians(deg)]
 def force(x):
  P,G,F,ok=direct_force(model,c['response'],np.atleast_2d(x),a,p)
  # Invalid low-pressure nodes remain NaN, never clipped or used as equilibria.
  return F
 return c,model,a,p,force

def roots_x(force,y=0.,z=0.):
 xx=np.linspace(-.55e-3,.55e-3,2201);q=np.c_[xx,np.full(len(xx),y),np.full(len(xx),z)];F=force(q)[:,0];roots=[]
 for j in np.flatnonzero((F[:-1]>0)&(F[1:]<=0)):
  roots.append(brentq(lambda u:force([[u,y,z]])[0,0],xx[j],xx[j+1],xtol=1e-13))
 return np.array(roots)
def main():
 start=time.perf_counter();verify();b=read(BASE/'BASELINE_RESULTS.json');assert b['PASS'] and len(b['rows'])==5;assert (BASE/'BUDGET_FREEZE.json').exists()
 records=[];planes=[];validations=[];zeta=6*np.pi*.001*R;cfg=read(CONFIG)['stability']
 for f in [2500000.,2400000.,2600000.]:
  for deg in [0,90,180,270]:
   c,model,a,p,force=force_context(f,deg);roots=roots_x(force);rr=[]
   for x in roots:
    q=np.array([x,0,0]);Js=[]
    for eps in cfg['derivative_steps_m']:
     J=np.column_stack([(force(q+np.eye(3)[j]*eps)[0]-force(q-np.eye(3)[j]*eps)[0])/(2*eps*zeta) for j in range(3)]);Js.append(J)
    ev,vec=np.linalg.eig(Js[-1]);error=np.linalg.norm(Js[-1]-Js[-2])/max(np.linalg.norm(Js[-1]),1e-30);F=force(q)[0];zero=np.linalg.norm(F)<cfg['force_zero_N'];tol=cfg['neutral_eigen_abs_s1']
    cls='NUMERICALLY_UNRESOLVED' if error>cfg['jacobian_relative_limit'] else ('FULL_3D_ATTRACTOR' if zero and max(ev.real)<-tol else ('TRANSVERSE_ATTRACTING_PLANE' if Js[-1][0,0]<-tol and max(abs(ev.real[np.argsort(abs(ev.real))[:2]]))<=tol else ('SADDLE_OR_UNSTABLE' if max(ev.real)>tol else 'X_TRANSVERSE_RESTORING')))
    rr.append(dict(position_m=q,force_N=F,drift_m_s=F/zeta,jacobians_s1=Js,eigenvalues_real_s1=ev.real,eigenvalues_imag_s1=ev.imag,eigenvectors=vec.real,derivative_relative_difference=error,force_zero=bool(zero),classification=cls,X_TRANSVERSE_RESTORING=bool(Js[-1][0,0]<0)))
   # Off-axis local X roots are not full force zeros.
   off=[]
   for y in [-350e-6,0,350e-6]:
    for z in [-60e-6,0,60e-6]:
     roots2=roots_x(force,y,z);off.append(dict(y_m=y,z_m=z,x_roots_m=roots2,force_N=force(np.c_[roots2,np.full(len(roots2),y),np.full(len(roots2),z)])))
   records.append(dict(frequency_Hz=f,phase_deg=deg,roots=rr,spacing_m=np.diff(roots),off_axis=off,source_amplitudes=a,phases=p))
   print('STATIC',f,deg,[r['classification'] for r in rr],flush=True)
  # Three-dimensional independent direct pressure/gradient validation.
  c,model,a,p,force=force_context(f);rng=np.random.default_rng(463010+int(f));q=rng.uniform([-.00055,-.00045,-.00007],[.00055,.00045,.00007],(96,3));P,G,F,ok=direct_force(model,c['response'],q,a,p);P48,G48,F48,_=direct_force(model,c['response'],q,a,p,48)
  eps=1e-9;fd=np.stack([(direct_force(model,c['response'],q+np.eye(3)[j]*eps,a,p)[0]-direct_force(model,c['response'],q-np.eye(3)[j]*eps,a,p)[0])/(2*eps) for j in range(3)],-1)
  v=dict(frequency=f,pressure_gradient_FD_rel=float(np.linalg.norm(fd-G)/np.linalg.norm(G)),quadrature_force_rel=float(np.linalg.norm(F48-F)/np.linalg.norm(F)),all_response_valid=bool(ok.all()))
  for cache in ['coarse','fine']:
   bs=CachedWaterBasis(BASE/'fields'/str(int(f))/(cache+'.h5'));pp,gg=bs.query_basis(q);u=np.array(a)*np.exp(1j*np.array(p));cp=u@pp;cg=np.einsum('i,inj->nj',u,gg);cf=-.5*np.real(c['response'].query(cp)[:,0,None]*cg.conj());v[cache+'_force_rel']=float(np.linalg.norm(cf-F)/np.linalg.norm(F))
  validations.append(v)
  yy,zz=np.meshgrid(np.linspace(-450e-6,450e-6,51),np.linspace(-70e-6,70e-6,25),indexing='ij');qyz=np.c_[np.zeros(yy.size),yy.ravel(),zz.ravel()];axis=[np.linspace(-550e-6,550e-6,25),np.linspace(-450e-6,450e-6,17),np.linspace(-70e-6,70e-6,7)];q3=np.stack(np.meshgrid(*axis,indexing='ij'),-1).reshape(-1,3)
  with h5py.File(BASE/f'FIELD_3D_{int(f)}.h5','x') as h:
   for key,q in [('grid',q3),('yz_center',qyz)]:
    P,G,F,ok=direct_force(model,c['response'],q,a,p)
    for name,vv in dict(positions_m=q,pressure_Pa=P,gradient_Pa_m=G,force_N=F,drift_m_s=F/zeta,response_valid=ok).items():h[key+'/'+name]=vv
   h.attrs.update(scope='SINGLE_BUBBLE_F1_STOKES_DIAGNOSTIC',identity_json=__import__('json').dumps(identity()))
 write(BASE/'STATIC_RESULTS.json',dict(rows=records,validation=validations,wall_s=time.perf_counter()-start,identity=identity(),scope='Single bubble F1 / 6 pi mu R, no F2, no streaming or buoyancy; local roots and eigenvalues only'))
if __name__=='__main__':main()
