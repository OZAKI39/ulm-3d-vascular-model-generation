"""Explicit sub-0.1 Pa dynamics and Boolean-safe final extension admission.

Discovery runs retain their recorded historical domain. This does not retune
actions or relabel old trajectories. All extra trajectories receive new IDs.
"""
import time,numpy as np
from .common import *
from .run import task,batch
from .baseline import compare_rows
def main():
    verify();start=time.perf_counter();old=read(BASE/'LOW_PRESSURE_FINAL.json')
    plans=[task('low_below_'+v,amps=[0,0,.003,.003],phase=[0]*4,variant=v,joined_diagnostic=True) for v in ['DT1','DT05','FINE']]
    write(BASE/'LOW_ADMISSION_PLAN.json',dict(plans=plans,reason='Earlier amplitude .1 Pa-eq gave incident pressure .37 Pa; require an explicit incident-pressure trajectory below .1 Pa. Preserve earlier records. NumPy scalar Boolean identity is not a numerical rejection.',identity=identity()))
    rows=batch(plans,3);checks=[compare_rows(rows[0],r) for r in rows[1:]]
    pressure_pass=all(r['status']=='PASS' and r['pressure_min_Pa']>=.001 and r['max_metrics']['max_pressure_Pa']<.1 for r in rows)
    gates={k:bool(v) for k,v in old['gates'].items()};gates['SUB_0P1_INCIDENT_PRESSURE_DYNAMIC_PASS']=pressure_pass and all(r['admitted'] for r in checks)
    result=dict(gates=gates,production_extension_admitted=all(gates.values()),admitted_min_Pa=.001,admitted_max_Pa=40000.,frequency_Hz=2.5e6,shell_multiplier=1.,rows=rows,comparisons=checks,prior_record_sha256=sha(BASE/'LOW_PRESSURE_FINAL.json'),scope='F10AE opt-in only; historical queries/files untouched. Discovery comparison retained the historical table. 0.0001 Pa excluded after timeout.',wall_s=time.perf_counter()-start,identity=identity())
    write(BASE/'LOW_PRESSURE_ADMISSION.json',result)
    print('LOW_ADMISSION',gates,flush=True)
if __name__=='__main__':main()
