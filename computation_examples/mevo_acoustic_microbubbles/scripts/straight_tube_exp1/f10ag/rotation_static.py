"""Ring force versus actual constrained velocity, no injected streaming."""
import time,json
import numpy as np,h5py
from mevo_mb.fast_sim.accelerated_water import WaterAccelerated
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from scripts.straight_tube_exp1.f10af.field import response_context
from scripts.straight_tube_exp1.f10ae.schedules import schedule
from scripts.straight_tube_exp1.f10af.metrics import rotation
from .common import *
def main():
 begin=time.perf_counter();verify();c=response_context();bs=CachedWaterBasis(BASE/'fields/2500000/coarse.h5');actions=read(BASE/'ACTIONS.json');rows=[];theta=np.arange(64)*2*np.pi/64;tan=np.c_[-np.sin(theta),np.cos(theta),np.zeros(len(theta))];rad=np.c_[np.cos(theta),np.sin(theta),np.zeros(len(theta))]
 with h5py.File(BASE/'ROTATION_STATIC.h5','x') as h:
  for a in ['ROT_PLUS','ROT_MINUS']:
   for radius in [50e-6,150e-6,300e-6,450e-6]:
    for z in [-60e-6,0,60e-6]:
     x=radius*rad+np.array([0,0,z]);sim=WaterAccelerated(x,c['physics'],c['response'],bs,c['cavity'],policy=policy());t=dict(actions[a],duration=.002);sim.control(schedule(t).segments[0].event());step=sim.step(1e-6);f=step['fields'];v=step['result']['velocity_m_s'];key=f'{a}_r{radius*1e6:.0f}_z{z*1e6:.0f}'
     for name,value in dict(positions_m=x,primary_force_N=f['primary_force_N'],secondary_force_N=f['secondary_force_N'],total_force_N=f['total_force_N'],constrained_velocity_m_s=v,pressure_Pa=f['pressure_Pa'],contact_force_N=step['result']['contact_force_N']).items():h[key+'/'+name]=value
     pp=np.r_[f['pressure_Pa'],f['pressure_Pa'][0]];F1t=np.sum(f['primary_force_N']*tan,axis=1);vt=np.sum((v-v.mean(0))*tan,axis=1);vr=np.sum((v-v.mean(0))*rad,axis=1)
     rows.append(dict(id=key,action=a,radius_m=radius,z_m=z,phases=t['phases'],mean_F1_tangent_N=float(F1t.mean()),rms_F1_tangent_N=float(np.sqrt(np.mean(F1t**2))),mean_Ftotal_tangent_N=float(np.sum(f['total_force_N']*tan,axis=1).mean()),mean_vtheta_m_s=float(vt.mean()),rms_vtheta_m_s=float(np.sqrt(np.mean(vt**2))),mean_vradial_m_s=float(vr.mean()),rms_vradial_m_s=float(np.sqrt(np.mean(vr**2))),positive_vtheta_fraction=float((vt>0).mean()),phase_winding=float(np.diff(np.unwrap(np.angle(pp))).sum()/(2*np.pi)),rotation=rotation(x,v),gates=step['gates'],accepted=step['accepted'],N=64))
     print('RING',key,rows[-1]['mean_vtheta_m_s'],flush=True)
 write(BASE/'ROTATION_STATIC.json',dict(rows=rows,wall_s=time.perf_counter()-begin,identity=identity(),scope='PRESCRIBED_RING_INITIAL_CONDITION; one actual full constrained step, not sustained rotation'))
if __name__=='__main__':main()
