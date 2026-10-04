"""Final idle CPU constitutive/q benchmark; no stationary solve or new law."""
import time
from types import SimpleNamespace
import numpy as np
from scripts import c7_early_branch as a
r=a.read(a.OUT/'GPU_BACKEND.json');weak=a.read(a.OUT/'activation/WEAK.json');data=a.arrays(a.OUT/f'intact/{weak["label"]}/STATE.npz');p=a.problem(7,weak['strain'])
F=p.state(data['x'],data['s']).Ftilde
class Probe:
    yc=p.yc
    def state(self,x,s):return SimpleNamespace(W0=a.old.spectral.svk_split_response(F,p.material).Psi_plus)
probe=Probe()
for _ in range(2):a.old.q_fields(probe,None,data['s'])
times=[]
for _ in range(10):
    start=time.perf_counter();a.old.q_fields(probe,None,data['s']);times.append(time.perf_counter()-start)
r['CPU_initial_concurrent_samples_seconds']=r['CPU_full_bond_batch_seconds'];r['CPU_full_bond_batch_seconds']=times;r['CPU_batch_median_seconds']=float(np.median(times));r['CPU_timing_context']='After all scientific jobs finish, no competing task worker pool; two warmups plus ten repetitions. Input F already assembled; original full spectral response and original q_fields code executed.';r['CPU_benchmark_state']=weak['strain']
a.save(a.OUT/'GPU_BACKEND.json',r);print('CPU batch median seconds',r['CPU_batch_median_seconds'])
