"""Independent field-cache check for single-bubble basin and point capture."""
import time
import numpy as np,h5py
from .common import *
from .single import integrate
from .capture import capture

def main():
 begin=time.perf_counter();verify();roots=np.array([r['position_m'][0] for v in read(BASE/'STATIC_RESULTS.json')['rows'] if v['frequency_Hz']==2500000 and v['phase_deg']==0 for r in v['roots']]);rows=[]
 write(BASE/'SINGLE_CACHE_PLAN.json',dict(sources=['atlas','local_3d_point'],cache_comparison=['coarse','fine'],position_limit_m=50e-9,capture_agreement_required=True,identity=identity()))
 for name in ['atlas','local_3d_point']:
  with h5py.File(BASE/'single'/(name+'.h5')) as h:x=h['positions_m'][0];times=h['time_s'][:];fine=h['positions_m'][:]
  result=integrate(name+'_coarse',2500000,0,x,cache='coarse')
  with h5py.File(BASE/result['trajectory_path']) as h:ct=h['time_s'][:];coarse=h['positions_m'][:]
  assert np.array_equal(times,ct);error=float(np.linalg.norm(fine-coarse,axis=2).max());outcomes=[]
  for full in [False,True]:
   for th in [10e-6,20e-6,30e-6]:
    a=capture(times,fine,roots,th,full_3d=full);b=capture(times,coarse,roots,th,full_3d=full);outcomes.append(dict(full_3d=full,threshold_m=th,disagreements=sum(p['right_censored']!=q['right_censored'] or p['final_band']!=q['final_band'] for p,q in zip(a,b))))
  rows.append(dict(id=name,max_cache_difference_m=error,capture=outcomes,PASS=error<=50e-9 and all(r['disagreements']==0 for r in outcomes)))
 write(BASE/'SINGLE_CACHE_RESULTS.json',dict(rows=rows,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
