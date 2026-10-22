"""Small paired pilot before admitting larger batches."""
import time
import numpy as np,h5py
from .common import *
from .run import task,perform
from .slow import integrate
from .states import main as states

def main():
    begin=time.perf_counter();verify()
    if not (BASE/'INITIAL_STATES.h5').exists():states()
    hist=read(AE/'runs/Ynode_-5e-06/SUMMARY.json');t=dict(hist['task']);t.update(id='historical_Ynode_replay',state='N64_compact_reference',historical_cache=True)
    plan=[t,task('pilot32_OFF',action='OFF',state='N32_wide_r0',duration=.0002),task('pilot64_X',duration=.002)]
    write(BASE/'PILOT_PLAN.json',dict(tasks=plan,identity=identity()))
    rows=[]
    for item in plan:
        r=perform(item);rows.append(r)
        if r['status']!='PASS':write(BASE/'PILOT_RESULTS.json',rows);raise RuntimeError(r.get('reason',r.get('reasons')))
    with h5py.File(AE/hist['trajectory_path']) as h, h5py.File(BASE/rows[0]['trajectory_path']) as z:
        same=np.array_equal(h['positions_m'][:],z['positions_m'][:]);difference=float(np.max(abs(h['positions_m'][:]-z['positions_m'][:])))
    write(BASE/'HISTORICAL_REPLAY.json',dict(bitwise_positions_equal=same,max_difference_m=difference,history_task=hist['task'],historical_sha256=hist['trajectory_sha256'],replay_sha256=rows[0]['trajectory_sha256'],identity=identity()))
    assert same
    for item in plan[1:]:rows.append(integrate(item))
    write(BASE/'PILOT_RESULTS.json',dict(rows=rows,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
