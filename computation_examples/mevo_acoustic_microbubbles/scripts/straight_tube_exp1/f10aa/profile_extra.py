"""Supplementary allocation/transfer accounting and isolated four-thread failures."""
import os,time,subprocess,sys,tracemalloc,psutil
import numpy as np,h5py
from mevo_mb.fast_sim.accelerated_water import WaterAccelerated
from mevo_mb.fast_sim.accelerated_fields import WaterPacked
from mevo_mb.fast_sim.accelerated_contact import GPUFactor
from .profile import Timers
from .trajectory import Campaign
from .common import *

def main():
    out=BASE/'profile_extra';out.mkdir(exist_ok=False);rows=[];start=time.perf_counter()
    for n in [32,64,128]:
      for backend in ['B1','PACKED','B3_CuPy']:
        c=Campaign(out/'scratch',backend);s=c.engine(n);s.control(c.schedules['PAIR'].segments[0].event());x=s.x.copy();s.step(1e-6);s.x=x.copy()
        transfers=dict(H2D_bytes=0,D2H_bytes=0,H2D_calls=0,D2H_calls=0);init=GPUFactor.__init__;solve=GPUFactor.solve
        def wrapped_init(self,A):
            transfers['H2D_bytes']+=A.nbytes;transfers['H2D_calls']+=1;return init(self,A)
        def wrapped_solve(self,b):
            transfers['H2D_bytes']+=b.nbytes;transfers['H2D_calls']+=1
            r=solve(self,b);transfers['D2H_bytes']+=r.nbytes;transfers['D2H_calls']+=1;return r
        GPUFactor.__init__=wrapped_init;GPUFactor.solve=wrapped_solve
        try:
            with Timers() as timer:r=s.step(1e-6)
        finally:GPUFactor.__init__=init;GPUFactor.solve=solve
        s.x=x.copy();tracemalloc.start();t=time.perf_counter();cpu=time.process_time();s.step(1e-6);wall=time.perf_counter()-t;cpu=time.process_time()-cpu;current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
        rows.append(dict(N=n,backend=backend,warp_transfer_bytes=timer.bytes,warp_transfer_calls=timer.calls,cupy_transfer=transfers,
            transfer_inclusive_wall_s={k:float(np.sum(v,axis=0)[0]) for k,v in timer.times.items()},tracemalloc_current_bytes=current,tracemalloc_peak_bytes=peak,
            instrumented_wall_s=wall,cpu_s=cpu,RSS_MiB=psutil.Process().memory_info().rss/2**20,actual_process_threads=psutil.Process().num_threads(),scope='Allocations and wrapper times instrumented; compare uninstrumented benchmark for speed'))
    write(out/'TRANSFER_ALLOCATION.json',rows)
    code='''import faulthandler;faulthandler.enable()
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
'''
    path=out/'four_thread_probe.py';path.write_text(code);isolated=[]
    for n in [32,64,128]:
      for backend in (['B1','B0'] if n==128 else ['B1']):
        env=dict(os.environ,OPENBLAS_NUM_THREADS='4',OMP_NUM_THREADS='4',CASE_N=str(n),CASE_BACKEND=backend)
        t=time.perf_counter();r=subprocess.run([sys.executable,str(path)],cwd=ROOT,env=env,capture_output=True,text=True,timeout=45)
        log=out/f'threads4_N{n}_{backend}.log';log.write_text(r.stdout+r.stderr)
        isolated.append(dict(N=n,backend=backend,exitcode=r.returncode,status='PASS' if r.returncode==0 else 'FAIL',elapsed_s=time.perf_counter()-t,log=log.name,stdout=r.stdout))
        write(out/'FOUR_THREADS.json',isolated)
    write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
