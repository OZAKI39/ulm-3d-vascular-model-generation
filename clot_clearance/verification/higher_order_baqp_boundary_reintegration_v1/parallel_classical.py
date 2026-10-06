"""Resource-proven overlap of the last two independent requested solves.

The original stdin dispatcher is paused (only its parent PID, never the solver
process group) so it cannot submit a duplicate classical case. It is retired
only after both scientific workers have finished. Solver configuration is not
changed. This script creates no additional equilibrium case.
"""
import json,subprocess,sys,time
from pathlib import Path
import psutil

V=Path(__file__).resolve().parent;ROOT=V.parents[1]
def read(n):return json.loads((V/n).read_text())
def save(n,x):(V/n).write_text(json.dumps(x,indent=2)+'\n')

def main():
    matches=[]
    for p in psutil.process_iter(['pid','cmdline']):
        a=p.info['cmdline'] or []
        if len(a)>=2 and a[-1]=='nonaffine_shear' and a[-2].endswith('higher_order_baqp_boundary_reintegration_v1/run_audit.py') and p.cwd()==str(ROOT):matches.append(p)
    assert len(matches)==1,'Expected exactly the active fourth worker'
    fourth=matches[0];dispatcher=fourth.parent();args=dispatcher.cmdline()
    assert args[-1]=='-' and 'python' in Path(args[0]).name and dispatcher.cwd()==str(ROOT),'Unexpected dispatcher identity'
    assert not (V/'CLASSICAL_SHEAR_SOLVE.log').exists(),'Fifth case already submitted'
    peak=max(read('RESOURCES_'+k+'.json')['peak_process_tree_RSS_bytes'] for k in ['normal_affine','shear_affine','nonaffine_axial'])
    mem=psutil.virtual_memory();reserve=1.5*peak;headroom=(mem.available-reserve)/mem.total
    assert headroom>=.2,'Insufficient proven headroom'
    decision=dict(fourth_pid=fourth.pid,paused_dispatcher_pid=dispatcher.pid,previous_measured_max_RSS_bytes=peak,
        additional_worker_reserve_factor=1.5,additional_reserved_bytes=reserve,available_before_bytes=mem.available,
        RAM_total_bytes=mem.total,predicted_headroom_fraction=headroom,minimum_required_headroom_fraction=.2,
        reason='Task section 42 permits two independent equilibrium solves after memory proof',solver_configuration_modified=False)
    dispatcher.suspend();save('PARALLEL_EXECUTION_DECISION.json',decision)
    start=time.time();minimum=1.;max_rss=0;fourth_finished=None
    with (V/'CLASSICAL_SHEAR_SOLVE.log').open('x') as log:
        proc=subprocess.Popen([sys.executable,'-B',str(V/'run_audit.py'),'classical_shear'],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        fifth=psutil.Process(proc.pid)
        while True:
            fourth_live=fourth.is_running() and fourth.status()!=psutil.STATUS_ZOMBIE
            if not fourth_live and fourth_finished is None:fourth_finished=time.time()
            fifth_live=proc.poll() is None
            rss=0
            for p in [fourth,fifth]:
                try:rss+=p.memory_info().rss
                except psutil.Error:pass
            max_rss=max(max_rss,rss);minimum=min(minimum,psutil.virtual_memory().available/mem.total)
            if not fourth_live and not fifth_live:break
            time.sleep(1.)
    # Both workers are finished before retiring the paused dispatcher; no live
    # scientific computation is killed and no duplicate case is launched.
    dispatcher.kill()
    summary=dict(simultaneous_equilibrium_workers=2,peak_combined_workers_RSS_bytes=max_rss,minimum_available_RAM_fraction=minimum,
        fifth_returncode=proc.returncode,coordination_wall_seconds=time.time()-start,
        paused_dispatcher_retired_after_both_workers_finished=True,scientific_workers_terminated=False,
        exactly_five_equilibrium_cases=True,solver_configuration_modified=False)
    save('PARALLEL_EXECUTION_SUMMARY.json',summary)
    assert proc.returncode==0 and minimum>=.2
    assert read('PD2_NONAFFINE_SHEAR.json')['PASS'] and read('PD2_CLASSICAL_SHEAR.json')['PASS']
    print('Both final shear cases completed; five total solves; no configuration changes',flush=True)

if __name__=='__main__':main()
