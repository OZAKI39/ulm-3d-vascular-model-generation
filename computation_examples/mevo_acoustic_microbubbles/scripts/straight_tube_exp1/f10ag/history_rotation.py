"""Recompute historical 50 ms negative rotation from every original accepted step."""
import time
import numpy as np,h5py
from .common import *
def main():
 begin=time.perf_counter();verify();old=read(AF/'SCIENTIFIC_RESULTS.json')['rows'];rows=[]
 for name in ['ROT_PLUS','ROT_MINUS']:
  id='long50_'+name;r=read(AF/'hf'/id/'SUMMARY.json');p=AF/r['trajectory_path'];assert sha(p)==r['trajectory_sha256'];ref=next(v for v in old if v['id']==id);omega=[];rigid=[]
  with h5py.File(p) as h:
   times=h['time_s'][:]
   for k in range(0,len(times)-1,512):
    x=h['positions_m'][k:min(k+512,len(times)-1)];v=h['force_samples/velocity_m_s'][k:min(k+512,len(times)-1)];q=x-x.mean(1)[:,None,:];w=v-v.mean(1)[:,None,:];cross=np.cross(q,w)[:,:,2].sum(1);om=cross/np.sum(q*q,axis=(1,2));omxy=cross/np.sum(q[:,:,:2]**2,axis=(1,2));fit=np.stack([-omxy[:,None]*q[:,:,1],omxy[:,None]*q[:,:,0],np.zeros(q.shape[:2])],-1);energy=np.sum(w*w,axis=(1,2));fraction=1-np.sum((w-fit)**2,axis=(1,2))/energy;omega.extend(om);rigid.extend(fraction)
   ix=np.rint(np.linspace(0,len(times)-1,401)).astype(int);positions=h['positions_m'][ix]
  omega=np.array(omega);rigid=np.array(rigid);angle=float(np.degrees(np.sum(omega*np.diff(times))));mean=float(rigid.mean());assert abs(angle-ref['cumulative_rotation_deg'])<1e-10 and abs(mean-ref['mean_rigid_fraction'])<1e-10
  with h5py.File(BASE/(id+'_historical_reanalysis.h5'),'x') as z:
   z['time_s']=times[:-1];z['Omega_z_s']=omega;z['rigid_fraction']=rigid;z['sample_time_s']=times[ix];z['sample_positions_m']=positions;z.attrs.update(scope='HISTORICAL_F10AF_DATA_REANALYSIS; not a new dynamics trajectory',raw_sha256=r['trajectory_sha256'])
  rows.append(dict(action=name,N=r['N'],duration_s=r['completed_s'],recomputed_angle_deg=angle,recomputed_mean_rigid_fraction=mean,historical_angle_deg=ref['cumulative_rotation_deg'],historical_rigid_fraction=ref['mean_rigid_fraction'],PASS=True,raw_path=str(p),raw_sha256=r['trajectory_sha256']))
 write(BASE/'HISTORICAL_ROTATION_REANALYSIS.json',dict(rows=rows,wall_s=time.perf_counter()-begin,scope='Independent all-step reanalysis, not new 50 ms dynamics. New frozen short ROT trajectories were separately reproduced bitwise.',identity=identity()))
if __name__=='__main__':main()
