import faulthandler;faulthandler.enable()
import os,h5py,numpy as np
from threadpoolctl import threadpool_limits
from scripts.straight_tube_exp1.f10aa.common import *
from mevo_mb.fast_sim import f10g_contact as old,accelerated_contact as opt
from mevo_mb.fast_sim.sequential_hard_sphere import sequential_step as original
from mevo_mb.fast_sim.accelerated_sequential import sequential_step as fast
backend=os.environ['BACKEND'];threads=int(os.environ['THREADS'])
with h5py.File(BASE/'baseline/REFERENCE.h5') as h:
 q=h['N128_cubic_ON'];x,F,A=q['x'][:],q['force_N'][:],q['A'][:];z=q.attrs['zeta']
 with threadpool_limits(limits=threads):
  for k in range(3):
   print(backend,threads,k,flush=True)
   a=old.resistance(x,R)['A'] if backend=='B0' else opt.resistance(x,R)['A']
   r=original(x,F,R,z,1e-6,a,policy()) if backend=='B0' else fast(x,F,R,z,1e-6,a,policy(),backend='cupy' if backend=='B3_CuPy' else 'cpu')
   print(r['converged'],r['active_count'],flush=True)
