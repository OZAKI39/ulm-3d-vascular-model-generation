"""One-shot, sequential, durable stage controller. No automatic run retries."""
from pathlib import Path
import argparse,csv,datetime,hashlib,json,os,re,signal,subprocess,sys,time,traceback

def put(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def tree_memory(pid):
    todo=[pid];rss=0;seen=set()
    while todo:
        x=todo.pop()
        if x in seen:continue
        seen.add(x)
        try:
            text=Path(f'/proc/{x}/status').read_text();m=re.search(r'^VmRSS:\s+(\d+)',text,re.M)
            if m:rss+=int(m[1])
            for task in Path(f'/proc/{x}/task').iterdir():
                todo.extend(int(v) for v in (task/'children').read_text().split())
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    return rss

def main(root,work):
    sys.path.insert(0,str(root/'src'))
    from finalize_lammps_particle_engine import integrity,evaluate_case,evaluate_all,kernel_summary
    state=root/'STAGE_STATE.json'
    assert not state.exists(),'Stage already exists; do not launch again'
    contract=json.loads((root/'contracts/MINIMAL_LAMMPS_PARTICLE_ENGINE_CONTRACT.json').read_text());g=contract['gates']
    assert integrity(root)['status']=='PASS'
    build=json.loads((root/'build_provenance/LAMMPS_BUILD_PROVENANCE.json').read_text());assert build['status']=='PASS'
    for k in ['cpu','gpu']:
        binary=work/('build_'+k)/'lmp';assert hashlib.sha256(binary.read_bytes()).hexdigest()==build['builds'][k]['sha256']
    completed=[];put(state,{'state':'RUNNING','completed':completed,'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    try:
        for sp in contract['cases']:
            key=sp['case'];case=root/'cases'/key
            assert not (case/'RUN_STARTED.json').exists(),f'{key} already launched'
            assert integrity(root)['status']=='PASS'
            cmd=['mpirun','--allow-run-as-root','--bind-to','none','-np',str(sp['mpi_ranks'])]
            if key=='case_f_kokkos_gpu':cmd+=['nsys','profile','--trace=cuda,nvtx','--sample=none','--cpuctxsw=none','--force-overwrite=false','--output',str(case/'gpu_trace')]
            cmd+=[str(work/('build_'+sp['engine'])/'lmp')]
            if sp['engine']=='gpu':cmd+=contract['gpu_flags']
            cmd+=['-in','in.lammps','-log','log.lammps']
            put(case/'RUN_STARTED.json',{'case':key,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'command':cmd,'cwd':str(case),'retry':False,'before_run_frozen_input_status':'PASS'})
            put(state,{'state':'RUNNING','current':key,'completed':completed})
            env=os.environ.copy();env['OMP_NUM_THREADS']='1';env['CUDA_VISIBLE_DEVICES']='0';env['HWLOC_COMPONENTS']='-gl'
            start=time.monotonic();peak=0;timeout=False;poll_count=0
            tele=None
            if sp['engine']=='gpu':
                # Persistent low-overhead telemetry process. No profiler install.
                telefile=(case/'NVIDIA_SMI_RAW.csv').open('w')
                tele=subprocess.Popen(['nvidia-smi','--query-gpu=timestamp,utilization.gpu,memory.used','--format=csv,noheader,nounits','-lms','100'],stdout=telefile,stderr=subprocess.DEVNULL,text=True)
            with (case/'stdout.log').open('w') as log,(case/'RESOURCE_TRACE.csv').open('w') as rf:
                writer=csv.DictWriter(rf,fieldnames=['elapsed_seconds','process_tree_rss_kib','gpu_utilization_percent','gpu_memory_mib']);writer.writeheader();rf.flush()
                p=subprocess.Popen(cmd,cwd=case,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
                while True:
                    now=time.monotonic()-start;rss=tree_memory(p.pid);peak=max(peak,rss);ut='';mem=''
                    if tele:
                        try:
                            ls=(case/'NVIDIA_SMI_RAW.csv').read_text().strip().splitlines()
                            if ls: _,ut,mem=[x.strip() for x in ls[-1].split(',')]
                        except (ValueError,FileNotFoundError):pass
                    writer.writerow({'elapsed_seconds':round(now,6),'process_tree_rss_kib':rss,'gpu_utilization_percent':ut,'gpu_memory_mib':mem});rf.flush();poll_count+=1
                    if p.poll() is not None:break
                    if now>g['run_wall_time_limit_s']:
                        timeout=True;os.killpg(p.pid,signal.SIGTERM)
                        try:p.wait(timeout=10)
                        except subprocess.TimeoutExpired:os.killpg(p.pid,signal.SIGKILL);p.wait()
                        break
                    time.sleep(contract['resource_sampling_seconds'])
            elapsed=time.monotonic()-start
            if tele:
                tele.terminate();tele.wait(timeout=5);telefile.close()
            metrics={'case':key,'returncode':p.returncode,'timed_out':timeout,'requested_steps':sp['steps'],'requested_particle_count':sp['N'],'wall_seconds':elapsed,'end_to_end_steps_per_s':sp['steps']/elapsed,'process_tree_peak_rss_kib':peak,'resource_sample_count':poll_count,'wall_time_scope':'mpirun command including initialization, sparse output, and profiler export-to-report overhead if enabled','profiled':key=='case_f_kokkos_gpu'}
            put(case/'RUN_METRICS.json',metrics)
            if p.returncode!=0 or timeout:raise RuntimeError(f'{key} process failed, no dependent run launched')
            if key=='case_f_kokkos_gpu':
                with (case/'nsys_export.log').open('w') as f:
                    r=subprocess.run(['nsys','export','--type','sqlite','--output',str(case/'gpu_trace.sqlite'),str(case/'gpu_trace.nsys-rep')],stdout=f,stderr=subprocess.STDOUT,timeout=g['nsys_export_time_limit_s'])
                assert r.returncode==0,'Nsight evidence export failed'
                summary=kernel_summary(case/'gpu_trace.sqlite');put(root/'validation/GPU_KERNEL_ACTIVITY.json',summary)
                with (root/'validation/GPU_KERNEL_SUMMARY.csv').open('w') as f:
                    wr=csv.DictWriter(f,fieldnames=['name','launches','total_device_seconds','first_ns','last_ns']);wr.writeheader();wr.writerows(summary.get('kernels',[]))
            evaluated=evaluate_case(root,key);put(case/'INDEPENDENT_CASE_VALIDATION.json',evaluated)
            if evaluated['status']!='PASS':raise RuntimeError(f'{key} independent gate failed: '+repr([k for k,v in evaluated['checks'].items() if not v])+' '+evaluated.get('exception',''))
            completed.append(key);print(key,'PASS',flush=True)
        final=evaluate_all(root);put(root/'validation/INDEPENDENT_FINALIZER.json',final)
        assert final['status']=='PASS','Final independent evaluation failed'
        put(state,{'state':'PASS','completed':completed,'end_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
    except Exception as e:
        put(state,{'state':'FAIL','completed':completed,'error':repr(e),'traceback':traceback.format_exc(),'automatic_retry':False,'dependent_cases_stopped':True});raise

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',required=True);ap.add_argument('--work',required=True);a=ap.parse_args();main(Path(a.root),Path(a.work))
