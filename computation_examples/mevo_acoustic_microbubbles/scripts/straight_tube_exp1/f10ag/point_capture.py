"""Explicit 3D distance-to-point confirmation, distinct from X-band membership."""
import time
import numpy as np,h5py
from .common import *
from .single import integrate
from .capture import capture

def main():
 begin=time.perf_counter();verify();rootrow=next(r for r in read(BASE/'STATIC_RESULTS.json')['rows'] if r['frequency_Hz']==2500000 and r['phase_deg']==0);roots=np.array([r['position_m'][0] for r in rootrow['roots']]);x=np.array([[root+sx*25e-6,sy*8e-6,sz*8e-6] for root in roots for sx in [-1,1] for sy in [-1,1] for sz in [-1,1]])
 write(BASE/'POINT_CAPTURE_PLAN.json',dict(initial_positions_m=x,duration_s=.25,point_distance_thresholds_m=[10e-6,20e-6,30e-6],hold_s=.01,scope='PRESCRIBED_LOCAL_3D_INITIAL_CONDITION; single bubble F1, no F2',identity=identity()));r=integrate('local_3d_point',2500000,0,x);rr=[]
 for name in ['atlas','local_3d_point']:
  with h5py.File(BASE/'single'/(name+'.h5')) as h:
   t=h['time_s'][:];q=h['positions_m'][:];rows={str(th):capture(t,q,roots,th*1e-6,full_3d=True) for th in [10,20,30]}
  rr.append(dict(id=name,capture=rows,transverse_displacement_m=(q[-1]-q[0])[:,1:],duration_s=float(t[-1])))
 write(BASE/'POINT_CAPTURE_RESULTS.json',dict(rows=rr,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
