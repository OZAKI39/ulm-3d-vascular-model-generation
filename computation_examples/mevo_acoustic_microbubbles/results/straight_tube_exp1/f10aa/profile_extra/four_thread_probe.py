import faulthandler;faulthandler.enable()
import os,time,json,h5py,numpy as np
from threadpoolctl import threadpool_info
from scripts.straight_tube_exp1.f10aa.common import *
from mevo_mb.fast_sim.accelerated_sequential import sequential_step
from mevo_mb.fast_sim.sequential_hard_sphere import sequential_step as original
from mevo_mb.fast_sim.accelerated_contact import resistance
n=int(os.environ['CASE_N']);backend=os.environ['CASE_BACKEND']
with h5py.File(BASE/'baseline/REFERENCE.h5') as h:
 q=h[f'N{n}_cubic_ON'];x=q['x'][:];F=q['force_N'][:];A=q['A'][:];z=q.attrs['zeta'];fn=original if backend=='B0' else sequential_step;tt=[]
 for k in range(11):
  t=time.perf_counter();r=fn(x,F,R,z,1e-6,A,policy());elapsed=time.perf_counter()-t
  assert r['converged'] and all(r['gates'].values())
  if k>=2:tt.append(elapsed)
 print(json.dumps(dict(N=n,backend=backend,threads=threadpool_info(),wall=stats(tt),status='PASS')))
