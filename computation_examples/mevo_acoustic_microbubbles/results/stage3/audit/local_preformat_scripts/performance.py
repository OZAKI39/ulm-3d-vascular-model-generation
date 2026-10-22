"""Same-machine implementation timings, not a physical validation criterion."""
import json
from pathlib import Path
import time
import numpy as np
from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from mevo_mb.bubbles.contact import contact_candidates
from mevo_mb.bubbles.auto_backend import AutoPairBackend
from scripts.stage3.common import new_run, finish

out=new_run('performance',dict(seed=20261002,particle_counts=[64,256,1024],contact_count=10000),['Stage2 same-server measured AutoPairBackend evidence'])
rng=np.random.default_rng(20261002);rows=[]
for n in [64,256,1024]:
    x=rng.uniform(-.001,.001,(n,3));r=np.full(n,1e-6);v=(rng.normal(size=(n,8))+1j*rng.normal(size=(n,8)))*1e-19
    values={}
    for backend in ['numpy','cupy']:
        model=ClassicalFarFieldBjerknes(backend,enforce_far_field=False)
        elapsed=[]
        for repeat in range(4):
            start=time.perf_counter();value=model.forces(x,r,v,2.5e6,1000)
            if hasattr(value,'get'):value=value.get()
            elapsed.append(time.perf_counter()-start)
        values[backend]=value
        rows.append(dict(N=n,backend=backend,median_wall_s=float(np.median(elapsed[1:])),includes='host input and returned host output; first call excluded'))
    error=float(np.linalg.norm(values['numpy']-values['cupy'])/np.linalg.norm(values['numpy']))
    rows[-1]['relative_L2_against_numpy']=error
    if error>1e-11:raise RuntimeError('CPU/GPU implementation mismatch')
x=rng.uniform(-.0001,.0001,(10000,3));r=rng.uniform(.5e-6,2e-6,10000);contacts={};timings={}
for backend in ['brute','ckdtree']:
    start=time.perf_counter();contacts[backend]=contact_candidates(x,r,1e-11,backend);timings[backend]=time.perf_counter()-start
np.testing.assert_array_equal(contacts['brute'],contacts['ckdtree'])
auto=AutoPairBackend(benchmark_path='results/stage2/performance/STATE.json')
finish(out,dict(status='PASS',pair_timings=rows,contact_timings_s=timings,contact_count=len(contacts['brute']),
    contact_missed=0,contact_extra=0,retained_auto_crossover=auto.crossover,
    auto_evidence='results/stage2/performance/STATE.json',scope='This server only; timings do not certify physical accuracy'))
