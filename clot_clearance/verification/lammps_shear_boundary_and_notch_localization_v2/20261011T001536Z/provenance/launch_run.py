"""Two registered boundary variants of the existing laboratory, with passive records."""
from pathlib import Path
import sys,csv,json,math
import numpy as np
from pd_clot.lammps_fragmentation.macro_lab_v1 import run,DELTA
from audit import regions

class Observer:
 def __init__(self,out,g,mask,sign,plan):
  self.out=out;self.g=g;self.X=g['X'];self.b=g['bonds'];self.mask=mask;self.sign=sign;self.plan=plan;self.previous=None;self.first=False
  self.low=self.X[:,2]<-.00035-1e-12;self.high=self.X[:,2]>.00035+1e-12;self.reg,self.dtip,self.dgrip,self.band=regions(self.X,self.b)
  d=self.X[self.b[:,1]]-self.X[self.b[:,0]]
  if plan['transverse_periodic_fixed_cell']:d[:,:2]-=.0006*np.rint(d[:,:2]/.0006)
  self.r=np.linalg.norm(d,axis=1);self.fitsteps={round(math.acos(1-2*v)*.01/math.pi/plan['dt_s']):'fit'+str(round(v*100)) for v in [.2,.5,.7]}
  self.fs=(out/'boundary_audit.csv').open('w');self.w=csv.writer(self.fs);self.w.writerow(['step','time_s','top_Rx_N','top_Ry_N','top_Rz_N','bottom_Rx_N','bottom_Ry_N','bottom_Rz_N','generalized_external_Q_N',*[f'F{i}{j}' for i in range(3) for j in range(3)],'boundary_displacement_error_m','boundary_velocity_error_m_s'])
  self.fr=(out/'regional_utilization.csv').open('w');self.wr=csv.writer(self.fr);self.wr.writerow(['step','time_s','notch_live_max','grip_live_max','surface_live_max','interior_live_max','notch_including_newbreak_max','grip_including_newbreak_max','surface_including_newbreak_max','interior_including_newbreak_max','broken_count'])
 def state(self,label,step,x,v,alive):np.savez_compressed(self.out/(label+'.npz'),step=step,time_s=step*self.plan['dt_s'],x_m=x,v_full_m_s=v,alive=alive)
 def record(self,step,x,v,fm,fd,alive,changed,F):
  t=step*self.plan['dt_s'];amp=self.plan['amplitude_m'];ramp=self.plan['ramp_s'];d=amp*.5*(1-math.cos(math.pi*t/ramp));vel=amp*.5*math.pi/ramp*math.sin(math.pi*t/ramp)
  disp_error=float(np.max(abs((x-self.X-d*self.sign)[self.mask])));vel_error=float(np.max(abs((v-vel*self.sign)[self.mask])));assert disp_error<1e-12 and vel_error<1e-10
  if step%20==0 or changed or step==self.plan['steps']:
   self.w.writerow([step,t,*(-fm[self.high].sum(0)),*(-fm[self.low].sum(0)),-np.sum(fm*self.sign*self.mask),*F.ravel(),disp_error,vel_error])
  if self.plan['case']=='M5B' and (step%20==0 or changed or step==self.plan['steps']):
   util=(np.linalg.norm(x[self.b[:,1]]-x[self.b[:,0]],axis=1)/self.r-1)/.15;eligible=alive.astype(bool);prior=self.previous[3].astype(bool) if self.previous is not None else eligible
   maxima=lambda mask:[float(util[(self.reg==j)&mask].max()) if np.any((self.reg==j)&mask) else float('nan') for j in range(4)]
   self.wr.writerow([step,t,*maxima(eligible),*maxima(prior),int((1-alive).sum())])
  if changed and not self.first:
   assert self.previous is not None;self.state('pre_first_break',*self.previous);self.first=True
  if self.plan['case']=='M2P' and step in self.fitsteps:self.state(self.fitsteps[step],step,x,v,alive)
  if self.plan['case']=='M5B' and step%2000==0:self.state(f'frame_{step:06d}',step,x,v,alive)
  self.previous=(step,x.copy(),v.copy(),alive.copy())
  if step%1000==0:self.fs.flush();self.fr.flush()
 def close(self):self.fs.close();self.fr.close()

if __name__=='__main__':
 case=sys.argv[1];out=Path(sys.argv[2]).resolve();assert case in ['M2P','M5B'];run(out,case,static='--static' in sys.argv,observer_factory=Observer)
