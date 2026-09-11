"""One job, sequential bounded tests, closed-file manifest, no automatic retry or deletion."""
from common import *
import subprocess,signal,resource,shutil,sys,fcntl,csv,math,traceback

def main():
    cfg=read(sys.argv[1]);root=Path(cfg['output_directory']);job=cfg['job_id'];start=time.monotonic()
    assert root.parent==Path('/workspace/bloodflow/cloud_runs') and root.name==job and root.resolve()==root
    assert not (root/'RESULTS_READY').exists(),'A completed job can never rerun'
    resource.setrlimit(resource.RLIMIT_CORE,(0,0))
    lock=Path('/workspace/bloodflow/metadata/cloud_smoke.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    env=dict(os.environ,HWLOC_COMPONENTS='-opencl,-gl',PYTHONDONTWRITEBYTECODE='1',CUDA_VISIBLE_DEVICES='0',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
        PATH='/workspace/bloodflow/.venv/bin:/usr/local/cuda/bin:/usr/bin:/bin',LD_LIBRARY_PATH='/usr/local/cuda/lib64',
        PYTHONPATH=cfg['cloud_project_root']+':'+str(root/'tools'),CLOUD_JOB_AUTHORIZED=cfg['cloud_plan_sha256'])
    env.pop('SINGLE_RBC_AUTHORIZED',None);env.pop('RBC_REPAIR_NATIVE_TRACE',None)
    ledger_path=Path(cfg['budget_ledger']);ledger=read(ledger_path) if ledger_path.exists() else dict(limit_s=180.,used_s=cfg['previous_checks_elapsed_s'],attempts=cfg['previous_checks'])
    state=dict(job_id=job,status='RUNNING',pid=os.getpid(),started_at=now(),tests={},source_hash=cfg['source_hash'],runtime_binary_hash=cfg['runtime_binary_hash'],physical_model_acceptance='UNCHANGED_NOT_MATCHED')
    write(root/'run_status.json',state);write(root/'budget_before.json',ledger)
    def run_test(name,argv,limit):
        remaining=ledger['limit_s']-ledger['used_s'];limit=min(limit,remaining)
        if limit<=3:raise RuntimeError('TOTAL_TEST_TIME_LIMIT')
        if not disk_allowed(shutil.disk_usage('/workspace').free,'run'):raise RuntimeError('DISK_LIMIT')
        d=root/name;d.mkdir();test_start=time.monotonic();record=dict(name=name,argv=argv,started_at=now(),limit_s=limit,status='RUNNING',gpus=1)
        state['tests'][name]=record;write(root/'run_status.json',state)
        with (d/'console.log').open('wb') as log:
            p=subprocess.Popen(argv,cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            record['pid']=p.pid;record['pgid']=p.pid;write(root/'run_status.json',state)
            waited=wait_guarded(p,max(.01,limit-(time.monotonic()-test_start)),lambda:shutil.disk_usage('/workspace').free)
            reason=waited['stop_reason']
        elapsed=time.monotonic()-test_start
        record.update(exit_code=p.returncode,elapsed_s=elapsed,status=reason or ('PASS' if p.returncode==0 else 'FAILED'),finished_at=now())
        ledger['attempts'].append(dict(job_id=job,**record));ledger['used_s']+=elapsed;write(ledger_path,ledger);write(root/'run_status.json',state)
        if p.returncode or reason:raise RuntimeError(name+':'+record['status'])
    terminal='FAILED'
    try:
        assert sha(cfg['runtime_binary_path'])==cfg['runtime_binary_hash'],'RUNTIME_BINARY_CHANGED'
        for relative,h in cfg['tool_hashes'].items():assert sha(root/'tools'/relative)==h,'CLOUD_TOOL_CHANGED'
        before=shutil.disk_usage('/workspace');write(root/'disk_before.json',dict(free=before.free,used=before.used))
        mpi=['/usr/bin/mpirun','--allow-run-as-root','--bind-to','none','-np','2',cfg['python'],'-B']
        run_test('A_loader',mpi+[str(root/'tools/smoke_worker.py'),str(root/'cloud_runtime.json'),'A'],15)
        loader=read(root/'A_loader/loader.json');assert loader['mpi_sum']==2.0
        run_test('B_liquid',mpi+[str(root/'tools/smoke_worker.py'),str(root/'cloud_runtime.json'),'B'],60)
        result=read(root/'B_liquid/completion.json');assert result['actual_steps']==cfg['liquid_steps'] and result['finite']
        if cfg['run_repair_path']:
            run_test('C_repair',mpi+[str(root/'tools/cloud_repair_worker.py'),'--spec',str(root/'repair_spec.json')],60)
            c=read(root/'C_repair/completion.json');assert c['completed'] and c['actual_steps']==250 and c['wall_hidden_steps']==4000 and c['shear_steps']==0 and c['bouncer']=='bounce_back' and c['membrane_updates']==250
            with (root/'C_repair/moments.csv').open() as f:rows=list(csv.DictReader(f))
            assert len(rows)>0 and all(math.isfinite(float(r[k])) for r in rows for k in ['temperature','max_speed'])
            assert int(rows[-1]['N'])==110592,'UNEXPECTED_FLUID_PARTICLE_COUNT'
        state['status']='CLOUD_SMOKE_PASS';terminal='COMPLETED'
    except Exception as ex:
        state['status']='CLOUD_SMOKE_FAILED';state['error']=type(ex).__name__+': '+str(ex)
        (root/'failure.txt').write_text(traceback.format_exc())
        if 'TIMEOUT' in str(ex):terminal='TIMEOUT'
        elif 'DISK_LIMIT' in str(ex):terminal='DISK_LIMIT'
    finally:
        state.update(finished_at=now(),elapsed_s=time.monotonic()-start,final_state=terminal,total_test_process_wall_s=ledger['used_s'],billing_note='Process wall time and deployment elapsed time are not Vast rental billing amounts')
        processes=[]
        groups={x.get('pgid') for x in state['tests'].values()}
        for p in Path('/proc').iterdir():
            if not p.name.isdigit():continue
            try:
                raw=(p/'stat').read_text();fields=raw[raw.rfind(')')+2:].split()
                if int(fields[2]) in groups and fields[0]!='Z':processes.append(dict(pid=int(p.name),state=fields[0]))
            except (FileNotFoundError,ProcessLookupError,PermissionError):pass
        state['remaining_test_processes']=processes
        if processes:state['status']='CLOUD_SMOKE_FAILED';state['error']='LINGERING_TEST_PROCESSES';terminal='FAILED'
        state['final_state']=terminal
        write(root/'run_status.json',state);write(root/'budget_after.json',ledger)
        disk=shutil.disk_usage('/workspace');write(root/'disk_after.json',dict(free=disk.free,used=disk.used,result_bytes=sum(x['size'] for x in files(root)),largest20=sorted(files(root),key=lambda x:x['size'],reverse=True)[:20]))
        seal(root,job,terminal)
    return 0 if terminal=='COMPLETED' else 2
if __name__=='__main__':sys.exit(main())
