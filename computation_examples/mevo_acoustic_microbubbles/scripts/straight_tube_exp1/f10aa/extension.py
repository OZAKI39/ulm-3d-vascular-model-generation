"""Conditional scale preflight and duration extension; no enlarged-N long trajectory."""
import time,traceback
import numpy as np,h5py
from threadpoolctl import threadpool_limits
from mevo_mb.fast_sim.accelerated_contact import contact_geometry
from scripts.straight_tube_exp1.f10z.prepare import cloud
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from scripts.straight_tube_exp1.f10z.run import Monitor
from .trajectory import Campaign
from .common import *

def main():
    verify();out=BASE/'extension';out.mkdir(exist_ok=False);start=time.perf_counter();cfg=read(CONFIG)
    gates=read(BASE/'full/COMPARISONS.json')
    admitted=all(v['admitted'] and v['event_topology_agrees'] for v in gates.values())
    results=read(BASE/'full/RUNS.json');probe=next(r for r in results if r['N']==128 and r['backend']=='B1' and r['variant']=='DT1' and r['state']=='dispersed' and r['action']=='PAIR')
    projected=probe['mean_step_s']*cfg['max_extension_s']/cfg['dt_s']*4
    rule=dict(admitted=admitted,projected_triple_s=projected,limit_s=120,N256_seed=18256,N256_minimum_separation_m=5e-6,N256_half_width_m=36e-6,N256_rule='Same fixed nominal density random cloud; preflight only, not scientific trajectory',extension_s=.002)
    write(out/'PRE_RUN_DECISION.json',rule)
    if not admitted or projected>=120:
        write(out/'COMPLETE.json',dict(status='NOT_RUN',reason='Registered eligibility gate',wall_s=time.perf_counter()-start));return
    mon=Monitor();mon.thread.start();records=[];checks={}
    with threadpool_limits(limits=1):
        c=Campaign(out/'preflight','B1');x=cloud(18256,36e-6,5e-6,n=256);x-=x.mean(0)
        with h5py.File(out/'N256_INITIAL.h5','w') as h:h['positions_m']=x
        c.states[256]=x;t=time.perf_counter();sim=c.engine(256);cold=time.perf_counter()-t;sim.control(c.schedules['PAIR'].segments[0].event())
        times=[];r=None
        try:
            for k in range(5):
                sim.x=x.copy();t=time.perf_counter();r=sim.step(1e-6);times.append(time.perf_counter()-t);assert r['accepted']
            write(out/'N256_PREFLIGHT.json',dict(status='PASS',N=256,pairs=32640,A_shape=[768,768],A_bytes=768**2*8,B_dense_bytes=32640*768*8,B_matrix_free_bytes=contact_geometry(x,R)[0].nbytes,cold_s=cold,wall=stats(times[1:]),gates=r['gates'],active=r['result']['active_count'],pressure_min_Pa=float(abs(r['fields']['pressure_Pa']).min()),pressure_max_Pa=float(abs(r['fields']['pressure_Pa']).max()),initial_sha256=sha(out/'N256_INITIAL.h5'),scope='Five reset steps; no forward long integration'))
        except Exception as e:
            write(out/'N256_PREFLIGHT.json',dict(status='FAIL',reason=repr(e)));(out/'N256_EXCEPTION.txt').write_text(traceback.format_exc())
        paths={}
        for backend,variants in [('B1',['DT1','DT05','FINE']),('B0',['DT1'])]:
            c=Campaign(out/backend,backend)
            for variant in variants:
                r=c.run(128,'PAIR',.002,variant);records.append(r);paths[backend+'_'+variant]=c.out/'runs'/r['name'];write(out/'RUNS.json',records)
                if r['status']!='PASS':break
        for other in ['B1_DT05','B1_FINE','B0_DT1']:
            checks[other]=compare(paths['B1_DT1'],paths[other]) if other in paths else dict(admitted=False,reason='NOT_RUN')
    mon.done.set();mon.thread.join(4);write(out/'resources.json',mon.rows);write(out/'COMPARISONS.json',checks)
    write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
