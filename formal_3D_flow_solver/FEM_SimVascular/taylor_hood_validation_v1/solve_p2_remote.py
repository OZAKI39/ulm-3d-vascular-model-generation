"""Native P2 GPU run with the unchanged baseline steady and health gates.

Existing scientific helpers are imported read-only. P2 measurements integrate
all ten velocity nodes on authoritative double-precision input geometry.
"""
from pathlib import Path
import sys,os,json,time,subprocess,hashlib,re,signal,argparse,xml.etree.ElementTree as ET
sys.dont_write_bytecode=True
sys.path[:0]=['/workspace/flow_mean_2p0_mmps_20260922/src','/workspace/flow_mean_2p0_mmps_20260922/scripts/sv13q']
from sv_validation.sv13 import SteadyStopMonitor
from sv_validation.sv13n import checkpoint_one_rank
from sv_validation.sv11 import linear_gate,nonlinear_gate
from flow_parser import parse_solver_log
from p2_measure import P2Measurements

def write(path,data):
    tmp=path.with_suffix(path.suffix+'.pending');tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n');tmp.replace(path)

def resource_sample(start):
    processes={}
    for d in Path('/proc').iterdir():
        if not d.name.isdigit():continue
        try:
            s=(d/'status').read_text();ppid=int(re.search(r'^PPid:\s*(\d+)',s,re.M)[1])
            rss=re.search(r'^VmRSS:\s*(\d+)',s,re.M);hwm=re.search(r'^VmHWM:\s*(\d+)',s,re.M)
            processes[int(d.name)]=dict(pid=int(d.name),ppid=ppid,rss_kib=int(rss[1]) if rss else 0,hwm_kib=int(hwm[1]) if hwm else 0)
        except (OSError,TypeError):pass
    ids={os.getpid()}
    for _ in range(8):ids|={pid for pid,v in processes.items() if v['ppid'] in ids}
    tree=[v for pid,v in processes.items() if pid in ids]
    gpu=subprocess.run(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],capture_output=True,text=True).stdout.strip()
    return dict(elapsed_s=time.time()-start,gpu=gpu,tree_RSS_MiB=sum(v['rss_kib'] for v in tree)/1024,
                max_process_HWM_MiB=max((v['hwm_kib'] for v in tree),default=0)/1024,processes=tree)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('root');parser.add_argument('--shift',action='store_true');parser.add_argument('--label');parser.add_argument('--solver-extra',default='');parser.add_argument('--inlet-value',type=float);a=parser.parse_args()
    root=Path(a.root);case=root/'case_int32';run=case/(a.label or ('production_shift' if a.shift else 'production_baseline'));run.mkdir()
    build=json.loads((root/'audit/remote_protection_manifest.json').read_text())['build'];w=build['PETSc_build'];exe=build['executable']
    assert hashlib.sha256(Path(exe).read_bytes()).hexdigest()=='0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7'
    (run/'solver.xml').write_bytes((case/'run/solver.xml').read_bytes())
    inlet_change=None
    if a.inlet_value is not None:
        tree=ET.parse(run/'solver.xml');node=tree.find(".//Add_BC[@name='INLET']/Value");old=float(node.text)
        assert a.inlet_value<0 and abs(a.inlet_value/old-1)<1e-6
        node.text=format(a.inlet_value,'.17g');tree.write(run/'solver.xml',encoding='utf-8',xml_declaration=True)
        inlet_change=dict(old=old,new=a.inlet_value,ratio=a.inlet_value/old,reason='Single allowed scalar correction of native nodal-normal flux; unchanged Flat/rim-zero profile')
    options=(root/'baseline_PETSC_OPTIONS.txt').read_text().strip()
    extra=' -sub_pc_factor_shift_type nonzero -sub_pc_factor_shift_amount 1e-6' if a.shift else ''
    extra+=' '+a.solver_extra if a.solver_extra else ''
    options+=extra;(run/'PETSC_OPTIONS.txt').write_text(options+'\n')
    write(run/'configuration_diff.json',dict(solver_only_options_added=extra,allowed_inlet_scalar_normalization=inlet_change,
        unchanged=['geometry','rho','mu','outlet pressures','WALL','dt','steady criteria','output cadence','tolerances']))
    policy=json.loads((case/'policy.json').read_text());dt=policy['dt_s'];measure=P2Measurements(case);monitor=SteadyStopMonitor(measure,policy)
    env={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
    env.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',LD_LIBRARY_PATH=build['runtime_library_path'],PETSC_OPTIONS=options)
    cmd=[w['candidate_wrapper'],w['launcher'],'-n','1',w['candidate_wrapper'],exe,'solver.xml']
    start=time.time();seen=set();stop=None;retries=[];samples=[];log=run/'solver.log';hard_failure=None
    with log.open('w') as f:
        proc=subprocess.Popen(cmd,cwd=run,env=env,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            for path in sorted((run/'1-procs').glob('result_*.vtu')):
                step=int(path.stem.rsplit('_',1)[1]);cp=path.with_name('stFile_%03d.bin'%step)
                if not step or step%policy['save_interval_steps'] or step in seen:continue
                if not cp.exists() or time.time()-max(cp.stat().st_mtime,path.stat().st_mtime)<3:continue
                try:
                    checkpoint_one_rank(cp,step,dt);u,p=measure.read(path);state=dict(measure.measure(u,p),step=step,time_s=step*dt)
                except Exception as e:
                    retries.append(dict(step=step,error=repr(e)));continue
                seen.add(step);history=parse_solver_log(log.read_text(errors='replace'),dt)
                failed=history['failed_linear_solves']+len(history['nonlinear_failure_messages'])
                qualified=monitor.observe(state,u,linear_failures=failed)
                if qualified and stop is None:stop=dict(reason='five consecutive steady intervals and physical gates',step=step,time_unix=time.time())
            if time.time()-start>policy['maximum_wall_time_s'] and stop is None:
                stop=dict(reason='wall time budget exceeded: FAIL',time_unix=time.time())
            text=log.read_text(errors='replace')
            if 'SV13Q_FATAL' in text or '[svMultiPhysics] ERROR' in text:
                hard_failure='native solver fatal error';stop=stop or dict(reason=hard_failure,time_unix=time.time())
            if stop and not (run/'STOP_SIM').exists():(run/'STOP_SIM').write_text('0\n')
            samples.append(resource_sample(start));write(run/'resource_samples.json',samples)
            write(run/'progress.json',dict(elapsed_s=time.time()-start,pid=proc.pid,states=monitor.states,intervals=monitor.intervals,stop=stop,read_retries=retries))
            # Give native stop a full step; bound unresponsive fatal/timeout cleanup.
            if stop and time.time()-stop['time_unix']>1800:
                hard_failure='native stop did not finish within 1800 seconds';os.killpg(proc.pid,signal.SIGTERM);break
            time.sleep(5)
        try:code=proc.wait(timeout=30)
        except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);code=proc.wait()
    text=log.read_text(errors='replace');history=parse_solver_log(text,dt)
    record=dict(exit_code=code,wall_time_s=time.time()-start,command=cmd,PETSC_OPTIONS=options,
                initial_state='zero',MPI_ranks=1,OMP_NUM_THREADS=1,history=history,states=monitor.states,intervals=monitor.intervals,
                stop=stop,hard_failure=hard_failure,read_retries=retries,resources=samples,
                peak_tree_RSS_MiB_sampled=max((s['tree_RSS_MiB'] for s in samples),default=0),
                peak_single_process_HWM_MiB=max((s['max_process_HWM_MiB'] for s in samples),default=0),
                GPU_matrix_verified='seqaijcusparse' in text,GPU_vector_verified='seqcuda' in text)
    failures=[]
    for name,fn in [('linear',lambda:linear_gate(record)),('nonlinear',lambda:nonlinear_gate(history))]:
        try:fn()
        except Exception as e:failures.append(dict(gate=name,error=str(e)))
    outputs=list((run/'1-procs').glob('result_*.vtu'))
    if outputs:
        path=max(outputs,key=lambda p:int(p.stem.rsplit('_',1)[1]));step=int(path.stem.rsplit('_',1)[1])
        try:
            u,p=measure.read(path);record.update(final_step=step,final_vtu=str(path),final=measure.measure(u,p),
                checkpoint=checkpoint_one_rank(path.with_name('stFile_%03d.bin'%step),step,dt))
        except Exception as e:failures.append(dict(gate='final_output',error=repr(e)))
    else:failures.append(dict(gate='final_output',error='No saved solution'))
    record['health_failures']=failures
    record['status']='PASS' if code==0 and not failures and not hard_failure and stop and 'step' in stop and record['GPU_matrix_verified'] and record['GPU_vector_verified'] else 'FAIL'
    write(run/'execution.json',record)
    print(json.dumps({k:v for k,v in record.items() if k not in ('history','resources','states','intervals')},indent=2),flush=True)

if __name__=='__main__':main()
