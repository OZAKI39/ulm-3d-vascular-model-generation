"""Bounded, output-matched complete trajectories, followed by independent refinements."""
import time,traceback,resource
import numpy as np
from threadpoolctl import threadpool_limits
from scripts.straight_tube_exp1.f10z.run import Monitor
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from .trajectory import Campaign
from .common import *

def main():
    verify();cfg=read(CONFIG);bench=BASE/'benchmark_v3'
    assert read(bench/'COMPLETE.json')['status']=='COMPLETE'
    fixed=read(bench/'FIXED_BENCHMARK.json')
    assert all(r['status']=='PASS' for r in fixed if r['backend']=='B1')
    out=BASE/'full';out.mkdir(exist_ok=False);start=time.perf_counter();mon=Monitor();mon.thread.start()
    records=[];comparisons={};status='EXCEPTION'
    try:
      with threadpool_limits(limits=1):
        for backend in ['B0','B1','PACKED','B3_CuPy']:
          if backend=='B3_CuPy' and not all(r['status']=='PASS' for r in fixed if r['backend']==backend):
            write(out/'B3_NOT_RUN.json',dict(status='NOT_RUN',reason='Fixed-state numerical failure'));continue
          c=Campaign(out/backend,backend)
          # Compile/warm each N with a reset frozen step outside formal trajectory timing.
          warm_start=time.perf_counter()
          for n in cfg['N']:
            sim=c.engine(n);sim.control(c.schedules['PAIR'].segments[0].event());assert sim.step(cfg['dt_s'])['accepted']
          write(c.out/'WARMUP.json',dict(wall_s=time.perf_counter()-warm_start))
          for n in cfg['N']:
            r=c.run(n,'PAIR',cfg['full_windows_s'][str(n)]);records.append(r)
            if backend!='B0':
              comparisons[f'{n}_{backend}']=compare(out/'B0/runs'/r['name'],c.out/'runs'/r['name'])
            write(out/'COMPARISONS.json',comparisons);write(out/'RUNS.json',records)
        c=Campaign(out/'B1','B1')
        for n in cfg['N']:
          duration=cfg['full_windows_s'][str(n)]
          records.append(c.run(n,'OFF',duration))
          for variant in ['DT05','FINE']:
            r=c.run(n,'PAIR',duration,variant);records.append(r)
            comparisons[f'{n}_{variant}']=compare(c.out/'runs'/r['name'].replace(variant,'DT1'),c.out/'runs'/r['name'])
            write(out/'COMPARISONS.json',comparisons);write(out/'RUNS.json',records)
        for state,n in [('original',32),('historical_end',64)]:
          paths={}
          for backend in ['B0','B1']:
            c=Campaign(out/(backend+'_'+state),backend,state)
            r=c.run(n,'PAIR',cfg['near_contact_window_s']);records.append(r);paths[backend]=c.out/'runs'/r['name']
          comparisons[f'{n}_{state}']=compare(paths['B0'],paths['B1'])
        status='COMPLETE'
    except Exception:
        (out/'EXCEPTION.txt').write_text(traceback.format_exc());raise
    finally:
        mon.done.set();mon.thread.join(4)
        write(out/'resources.json',mon.rows);write(out/'COMPARISONS.json',comparisons);write(out/'RUNS.json',records)
        write(out/'COMPLETE.json',dict(status=status,wall_s=time.perf_counter()-start,identity=identity(),cpu_user_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime))
if __name__=='__main__':main()
