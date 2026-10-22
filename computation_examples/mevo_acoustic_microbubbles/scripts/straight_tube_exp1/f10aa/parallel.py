"""Four independent physical jobs, bounded one/two/four-process throughput experiment."""
import os,time,multiprocessing,traceback
from concurrent.futures import ProcessPoolExecutor
from threadpoolctl import threadpool_limits
from .common import *

def worker(item):
    # BLAS variables are set in the launcher before Python/NumPy import.
    from .trajectory import Campaign
    index,workers,round_id,root=item
    begin=time.perf_counter();out=Path(root)/f'w{workers}_{round_id}_job{index}'
    with threadpool_limits(limits=1):
        state='dispersed' if index<2 else 'historical_end'
        action='OFF' if index%2==0 else 'PAIR'
        c=Campaign(out,'B1',state)
        sim=c.engine(64);sim.control(c.schedules[action].segments[0].event());assert sim.step(1e-6)['accepted']
        warm=time.perf_counter()-begin
        r=c.run(64,action,read(CONFIG)['parallel_window_s'])
    return dict(job=index,workers=workers,round=round_id,PID=os.getpid(),warm_s=warm,latency_s=time.perf_counter()-begin,result=r,path=str(out.relative_to(BASE)))

def main():
    from scripts.straight_tube_exp1.f10z.run import Monitor
    from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
    assert read(BASE/'full/COMPLETE.json')['status']=='COMPLETE'
    verify();cfg=read(CONFIG);out=BASE/'parallel';out.mkdir(exist_ok=False);allrows=[];start=time.perf_counter()
    mon=Monitor();mon.thread.start();comparisons=[]
    # Spawn, never fork an initialized CUDA context. Same jobs, same outputs at every concurrency.
    for workers in cfg['concurrent_workers']:
        launched=time.perf_counter()
        with ProcessPoolExecutor(workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            for round_id in ['cold','warm']:
                t=time.perf_counter();jobs=list(pool.map(worker,[(j,workers,round_id,str(out)) for j in range(cfg['parallel_jobs'])]));elapsed=time.perf_counter()-t
                row=dict(workers=workers,round=round_id,batch_wall_s=elapsed,total_physical_s=sum(j['result']['completed_s'] for j in jobs),jobs=jobs,
                    throughput_physical_s_per_wall_s=sum(j['result']['completed_s'] for j in jobs)/elapsed,
                    all_pass=all(j['result']['status']=='PASS' for j in jobs))
                allrows.append(row);write(out/'BATCHES.json',allrows);print(workers,round_id,elapsed,row['all_pass'],flush=True)
        if workers>1:
            for j in range(cfg['parallel_jobs']):
                name=allrows[-1]['jobs'][j]['result']['name']
                comparisons.append(dict(workers=workers,job=j,comparison=compare(out/f'w1_warm_job{j}/runs'/name,out/f'w{workers}_warm_job{j}/runs'/name)))
        write(out/'COMPARISONS.json',comparisons)
    mon.done.set();mon.thread.join(4);write(out/'resources.json',mon.rows)
    write(out/'COMPLETE.json',dict(status='COMPLETE',wall_s=time.perf_counter()-start,identity=identity(),CPU_monitor_scope='Parent process only; worker CPU time and RSS are in each RESULT; nvidia-smi is whole device'))
if __name__=='__main__':main()
