"""Time ladder with same-duration paired gates; no unvalidated extrapolation."""
import time
import numpy as np,h5py
from .common import *
from .run import task,batch
from .slow import batch as slow_batch
from .paired import compare

def compare_hf(a,b):
    if a['status']!='PASS' or b['status']!='PASS':return dict(reference=a['id'],comparison=b['id'],PASS=False,reason='Incomplete trajectory')
    with h5py.File(BASE/a['trajectory_path']) as h,h5py.File(BASE/b['trajectory_path']) as z:
        t=h['time_s'][:];u=z['time_s'][:];ix=np.unique(np.rint(np.linspace(0,len(t)-1,201)).astype(int));errors=[]
        for j in ix:
            k=int(np.argmin(abs(u-t[j])));assert abs(u[k]-t[j])<1e-12;errors.append(np.linalg.norm(h['positions_m'][j]-z['positions_m'][k],axis=1).max())
    return dict(reference=a['id'],comparison=b['id'],max_fixed_ID_m=float(max(errors)),PASS=max(errors)<=50e-9)

def main():
    begin=time.perf_counter();old=read(BASE/'PAIRED_VALIDATION.json');valid=[r for r in old['rows'] if r['id']!='pair_compact'];assert all(r['PASS'] for r in valid)
    # Matching failure in compact state is evidence of safe fallback, not admission.
    failure=[]
    with h5py.File(BASE/'hf/pair_compact/runs/N64_CONTROL_2000us_DT1/trajectory.h5') as h,h5py.File(BASE/'lf/pair_compact/trajectory.h5') as z:
        diff=float(np.max(abs(h['positions_m'][:]-z['positions_m'][:])))
    write(BASE/'RESTRICTED_LF_ADMISSION.json',dict(PASS=True,allowed_states='Registered dispersed states in central cache; contact fallback checked separately',excluded='Dense compact state fails original reaction-balance gate; not admitted',compact_failure_same_prefix=diff,successful_pair_count=len(valid),maximum_difference_m=max(r['max_fixed_ID_m'] for r in valid),screening_for_pattern_success='INCONCLUSIVE_NO_HF_SUCCESSES',identity=identity()))
    results=[];alltasks=[];checks=[]
    for ms in [5,10,20]:
        tasks=[task(f'extend_{ms}_{a}',a,duration=ms/1000) for a in ['X_PAIR','XY_QUADRATURE','ROT_PLUS','ROT_MINUS']]
        alltasks+=tasks;write(BASE/f'EXTENSION_{ms}_PLAN.json',dict(tasks=tasks,gate='Same-duration LF/HF; stopped at first failed level',identity=identity()))
        rr=batch(tasks);slow_batch(tasks);cc=[compare(t) for t in tasks];results+=rr;checks+=cc;write(BASE/'EXTENSION_PROGRESS.json',dict(rows=results,paired=checks,wall_s=time.perf_counter()-begin,identity=identity()))
        if not all(r['PASS'] for r in cc):break
    refs=[r for r in results if r['task']['duration']==.02 and r['task']['action'] in ['X_PAIR','ROT_MINUS']]
    sensitivity=[]
    for r in refs:
        for v in ['DT05','FINE']:
            t=dict(r['task'],id=r['id']+'_'+v,variant=v);sensitivity.append(t)
    write(BASE/'EXTENSION_REFINEMENT_PLAN.json',dict(tasks=sensitivity,identity=identity()));fine=batch(sensitivity);comparisons=[compare_hf(r,z) for r in refs for z in fine if z['id'].startswith(r['id']+'_')]
    write(BASE/'EXTENSION_RESULTS.json',dict(rows=results,paired=checks,refinements=fine,convergence=comparisons,PASS=all(r['PASS'] for r in checks+comparisons),wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
