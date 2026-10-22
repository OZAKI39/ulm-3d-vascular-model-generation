"""Fixed-startup BLAS threads: safely test the alternative to crashing runtime resize."""
import os,sys,time,subprocess,json
import numpy as np,h5py
from threadpoolctl import threadpool_info
from mevo_mb.fast_sim.accelerated_contact import resistance
from mevo_mb.fast_sim.accelerated_sequential import sequential_step
from .trajectory import Campaign
from .common import *

def worker(threads):
    assert all(v['num_threads']==threads for v in threadpool_info())
    out=BASE/f'threaded/B2_t{threads}';out.mkdir(parents=True,exist_ok=False);rows=[];start=time.perf_counter()
    with h5py.File(BASE/'baseline/REFERENCE.h5') as h:
      for name,q in h.items():
        x=q['x'][:];F=q['force_N'][:];A=q['A'][:];z=q.attrs['zeta'];times=[]
        for k in range(11):
            t=time.perf_counter();r=sequential_step(x,F,R,z,1e-6,resistance(x,R)['A'],policy());dt=time.perf_counter()-t
            assert r['converged'] and all(r['gates'].values())
            err=float(np.linalg.norm(r['positions_m']-q['reference/positions_m'][:],axis=1).max());assert err<=policy()['fixed_id_tolerance_m']
            if k>=2:times.append(dt)
        rows.append(dict(case=name,threads=threads,status='PASS',wall=stats(times),max_fixed_ID_m=err));write(out/'FIXED.json',rows)
    c=Campaign(out,'B1')
    for n in [32,64,128]:
        s=c.engine(n);s.control(c.schedules['PAIR'].segments[0].event());assert s.step(1e-6)['accepted']
        r=c.run(n,'PAIR',read(CONFIG)['full_windows_s'][str(n)])
    write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,BLAS=threadpool_info(),identity=identity(),backend='B1',threads=threads))

def main():
    from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
    out=BASE/'threaded';out.mkdir(exist_ok=False);start=time.perf_counter();rows=[];checks={}
    for threads in [2,4]:
        env=dict(os.environ,OPENBLAS_NUM_THREADS=str(threads),OMP_NUM_THREADS=str(threads),MKL_NUM_THREADS=str(threads))
        with (out/f't{threads}.log').open('w') as f:
            r=subprocess.run([sys.executable,'-m','scripts.straight_tube_exp1.f10aa.threaded','--worker',str(threads)],env=env,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,timeout=240)
        rows.append(dict(threads=threads,exitcode=r.returncode));write(out/'STATUS.json',rows)
        if r.returncode==0:
            for path in (out/f'B2_t{threads}/runs').glob('*/RESULT.json'):
                checks[f't{threads}_'+path.parent.name]=compare(BASE/'full/B0/runs'/path.parent.name,path.parent)
    write(out/'COMPARISONS.json',checks);write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':
    if '--worker' in sys.argv:worker(int(sys.argv[-1]))
    else:main()
