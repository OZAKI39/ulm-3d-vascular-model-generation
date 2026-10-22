"""Bounded 50/100 ms extensions, with per-action cost and numerical admission."""
import time
from .common import *
from .run import task,batch
from .slow import batch as slow_batch
from .paired import compare
from .extension import compare_hf

def main():
    begin=time.perf_counter();assert read(BASE/'EXTENSION_RESULTS.json')['PASS'];allrows=[];allchecks=[];allfine=[];allconv=[]
    tasks=[task('long50_'+a,a,duration=.05) for a in ['X_PAIR','XY_QUADRATURE','ROT_PLUS','ROT_MINUS']]
    write(BASE/'LONG50_PLAN.json',dict(tasks=tasks,identity=identity()));rows=batch(tasks);slow_batch(tasks);checks=[compare(t) for t in tasks];allrows+=rows;allchecks+=checks
    refs=[r for r in rows if r['status']=='PASS' and r['task']['action']=='X_PAIR'];tt=[dict(r['task'],id=r['id']+'_'+v,variant=v) for r in refs for v in ['DT05','FINE']];write(BASE/'LONG50_REFINEMENT_PLAN.json',dict(tasks=tt,identity=identity()));fine=batch(tt);conv=[compare_hf(r,z) for r in refs for z in fine];allfine+=fine;allconv+=conv
    allowed=[r for r in rows if r['status']=='PASS' and r['task']['action'] in ['X_PAIR','ROT_PLUS'] and next(c['PASS'] for c in checks if c['id']==r['id']) and r['wall_s']*2.2<read(CONFIG)['budget']['per_HF_wall_s']]
    # Require room for the independently refined longest reference, not just DT1.
    dtref=next((r for r in fine if r['id'].endswith('DT05')),None)
    if not(dtref and dtref['status']=='PASS' and dtref['wall_s']*2.2<read(CONFIG)['budget']['per_HF_wall_s'] and all(z['PASS'] for z in conv)):allowed=[]
    tasks=[task('long100_'+r['task']['action'],r['task']['action'],duration=.1) for r in allowed];write(BASE/'LONG100_PLAN.json',dict(tasks=tasks,reason='Cost and 50 ms convergence gate',identity=identity()))
    if tasks:
        rows=batch(tasks);slow_batch(tasks);checks=[compare(t) for t in tasks];allrows+=rows;allchecks+=checks;refs=[r for r in rows if r['status']=='PASS' and r['task']['action']=='X_PAIR'];tt=[dict(r['task'],id=r['id']+'_'+v,variant=v) for r in refs for v in ['DT05','FINE']];write(BASE/'LONG100_REFINEMENT_PLAN.json',dict(tasks=tt,identity=identity()));fine=batch(tt);conv=[compare_hf(r,z) for r in refs for z in fine];allfine+=fine;allconv+=conv
    write(BASE/'LONG_TIME_RESULTS.json',dict(rows=allrows,paired=allchecks,refinements=allfine,convergence=allconv,wall_s=time.perf_counter()-begin,identity=identity(),not_run=[.5,1.],scope='Restricted zero-flow force-only model; no buoyancy, streaming or experimental lifetime calibration'))
if __name__=='__main__':main()
