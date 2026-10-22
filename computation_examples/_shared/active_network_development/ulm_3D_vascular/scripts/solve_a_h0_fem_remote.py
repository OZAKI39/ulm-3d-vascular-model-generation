"""H0 case runner adapted from the unchanged production runner.
Source SHA256: d458c0c97a9f6d61a3f98e1d12e5e40ac6877b70aca1173be87bda3d19417ed8
Numerical command, options, policy and gates are reused without changes.
"""
from pathlib import Path
import sys, os, json, time, subprocess, socket, hashlib, traceback, re, argparse, resource
parser=argparse.ArgumentParser()
parser.add_argument('--case', type=Path, required=True)
parser.add_argument('--reference-root', type=Path, required=True)
args=parser.parse_args()
ROOT = args.reference_root
sys.path[:0] = [str(ROOT/'src'), str(ROOT/'scripts/sv13q')]
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.sv13 import SteadyStopMonitor
from sv_validation.sv13n import checkpoint_one_rank
from sv_validation.sv11 import linear_gate, nonlinear_gate
from flow_parser import parse_solver_log

def write(path, data):
    tmp = path.with_suffix(path.suffix+'.pending')
    tmp.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n'); tmp.replace(path)
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def capture(args): return subprocess.check_output(args,text=True).strip()

def main():
    case = args.case; run = case/'run'; reports = case/'reports'
    assert not (run/'solver.log').exists() and not (run/'STOP_SIM').exists(), 'No overwrite permitted'
    for name, expected in json.loads((case/'input_hashes.json').read_text()).items():
        assert sha(case/name)==expected, name
    policy=json.loads((case/'policy.json').read_text()); dt=policy['dt_s']
    stack=Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy')
    build=json.loads((stack/'sv1_3q/reports/svmp_reuse_build.json').read_text())
    exe=build['executable']; w=build['PETSc_build']
    assert sha(exe)=='0e509fe21424b5b1f731c4b4f519839b0104bfb2aab0817b01556cc4054f7fc7'
    assert sha(Path(w['prefix'])/'lib/libpetsc.so')=='b79e98d278e105aa8a1848bd08d66ef36a3e2425458610f97213d3fc8e45c5ef'
    options=(run/'PETSC_OPTIONS.txt').read_text().strip()
    env={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
    env.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',
        CUDA_VISIBLE_DEVICES='0', LD_LIBRARY_PATH=build['runtime_library_path'], PETSC_OPTIONS=options)
    command=[w['candidate_wrapper'],w['launcher'],'-n','1',w['candidate_wrapper'],exe,'solver.xml']
    host=dict(hostname=socket.gethostname(),remote=True,cpu=capture(['lscpu']),ram=capture(['free','-b']),
        gpu=capture(['nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader']),
        cgroup={p:Path(p).read_text().strip() for p in ('/sys/fs/cgroup/cpu.max','/sys/fs/cgroup/memory.max')},
        output_path=str(case),python=sys.version)
    write(reports/'host.json',host)
    measure=SolutionMeasurements(case/'SV_MESH/mesh_arrays.npz',policy['Q_target_m3_s'],policy['Umean_m_s'])
    monitor=SteadyStopMonitor(measure,policy); seen=set(); stop=None; errors=[]; samples=[]
    start=time.time(); log=run/'solver.log'; peak_tree_rss_kib=0; peak_process_hwm_kib=0
    with log.open('w') as out:
        proc=subprocess.Popen(command,cwd=run,env=env,stdout=out,stderr=subprocess.STDOUT,start_new_session=True)
        while proc.poll() is None:
            for p in sorted((run/'1-procs').glob('result_*.vtu')):
                step=int(p.stem.rsplit('_',1)[1])
                if step==0 or step%10 or step in seen: continue
                cp=p.with_name('stFile_%03d.bin'%step)
                if not cp.exists() or time.time()-max(p.stat().st_mtime,cp.stat().st_mtime)<2: continue
                try:
                    checkpoint_one_rank(cp,step,dt); u,pressure=measure.read(p)
                except Exception as exc:
                    errors.append(dict(step=step,error=str(exc),retry=True)); continue
                state=dict(measure.measure(u,pressure),step=step,time_s=step*dt)
                seen.add(step); qualified=monitor.observe(state,u)
                if qualified and stop is None:
                    stop=dict(reason='five consecutive steady intervals and physical gates',step=step,time_unix=time.time())
                    tmp=run/'STOP_SIM.pending'; tmp.write_text('0\n'); tmp.replace(run/'STOP_SIM')
            if time.time()-start>policy['maximum_wall_time_s'] and stop is None:
                stop=dict(reason='wall time budget exceeded: FAIL',time_unix=time.time())
                tmp=run/'STOP_SIM.pending';tmp.write_text('0\n');tmp.replace(run/'STOP_SIM')
            # /proc memory sampling: sum resident memory in this runner's descendant tree.
            info={}
            for status_path in Path('/proc').glob('[0-9]*/status'):
                try:
                    fields=dict(line.split(':',1) for line in status_path.read_text().splitlines() if ':' in line)
                    info[int(status_path.parent.name)]=(int(fields['PPid']),int(fields.get('VmRSS','0 kB').split()[0]),int(fields.get('VmHWM','0 kB').split()[0]))
                except (OSError,ValueError,KeyError): pass
            descendants={os.getpid()}
            while True:
                more={pid for pid,(ppid,_,_) in info.items() if ppid in descendants}
                if more<=descendants:break
                descendants|=more
            peak_tree_rss_kib=max(peak_tree_rss_kib,sum(info.get(pid,(0,0,0))[1] for pid in descendants))
            peak_process_hwm_kib=max(peak_process_hwm_kib,max(info.get(pid,(0,0,0))[2] for pid in descendants))
            sample=capture(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'])
            samples.append(dict(elapsed_s=time.time()-start,gpu=sample))
            write(reports/'progress.json',dict(elapsed_s=time.time()-start,pid=proc.pid,states=monitor.states,
                  intervals=monitor.intervals,stop=stop,read_retries=errors))
            time.sleep(5)
        code=proc.wait()
    history=parse_solver_log(log.read_text(errors='replace'),dt)
    record=dict(exit_code=code,wall_time_s=time.time()-start,start_unix=start,end_unix=time.time(),
        command=command,solver_sha256=sha(exe),PETSc_library_sha256=sha(Path(w['prefix'])/'lib/libpetsc.so'),
        PETSC_OPTIONS=options,initial_state='zero',MPI_ranks=1,OMP_NUM_THREADS=1,history=history,
        stop=stop,states=monitor.states,intervals=monitor.intervals,GPU_samples=samples,log_sha256=sha(log), peak_tree_RSS_MiB_sampled=peak_tree_rss_kib/1024, peak_single_process_HWM_MiB=peak_process_hwm_kib/1024, child_rusage_maxrss_MiB=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024, memory_method='5-second /proc sampling; tree RSS includes runner and descendants; process HWM and wait4 child high-water also recorded')
    failures=[]
    for name,fn in [('linear',lambda:linear_gate(record)),('nonlinear',lambda:nonlinear_gate(history))]:
        try: fn()
        except Exception as exc: failures.append(dict(gate=name,error=str(exc)))
    outputs=sorted((run/'1-procs').glob('result_*.vtu'))
    final=max(outputs,key=lambda p:int(p.stem.rsplit('_',1)[1]))
    step=int(final.stem.rsplit('_',1)[1]); u,p=measure.read(final)
    record.update(final_step=step,final_vtu=str(final.relative_to(case)),final=measure.measure(u,p),
        checkpoint=checkpoint_one_rank(final.with_name('stFile_%03d.bin'%step),step,dt),health_failures=failures)
    record['status']='PASS' if code==0 and not failures and stop and 'step' in stop else 'FAIL'
    write(reports/'execution.json',record)
    print(json.dumps({k:record[k] for k in ('status','exit_code','wall_time_s','final_step','final','health_failures')},indent=2),flush=True)

if __name__=='__main__':
    try: main()
    except Exception:
        traceback.print_exc(); raise
