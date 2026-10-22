"""Independent states, frequencies and particle counts; no optimizer decisions."""
import time
from .common import *
from .run import task,batch
from .slow import batch as slow_batch
from .paired import compare
from .extension import compare_hf

def main():
    begin=time.perf_counter();assert read(BASE/'RESTRICTED_LF_ADMISSION.json')['PASS']
    tasks=[task(f'general_v{j}_{a}',a,state=f'N64_wide_v{j}',duration=.02) for j in range(3) for a in ['X_PAIR','ROT_PLUS','ROT_MINUS']]
    tasks += [task(f'scale_N{n}','X_PAIR',f'N{n}_wide_r0',duration=.002) for n in [32,128,200]]
    tasks += [task(f'frequency_{f}',action='X_PAIR',duration=.02,frequency=f) for f in [2400000.,2600000.]]
    tasks += [task('position_'+a,a,duration=.02,offset_m=[0,50e-6,0]) for a in ['ROT_PLUS','ROT_MINUS']]
    write(BASE/'GENERALIZATION_PLAN.json',dict(tasks=tasks,identity=identity(),scope='Independent fixed hypotheses; no automatic LF pruning'))
    rr=batch(tasks);slow_batch(tasks);checks=[compare(t) for t in tasks]
    sensitivity=[dict(t,id=t['id']+'_'+v,variant=v) for t in tasks if t['id'].startswith('frequency_') for v in ['DT05','FINE']]
    write(BASE/'FREQUENCY_REFINEMENT_PLAN.json',dict(tasks=sensitivity,identity=identity()));fine=batch(sensitivity);convergence=[compare_hf(r,z) for r in rr for z in fine if z['id'].startswith(r['id']+'_')]
    write(BASE/'GENERALIZATION_RESULTS.json',dict(rows=rr,paired=checks,refinements=fine,convergence=convergence,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
