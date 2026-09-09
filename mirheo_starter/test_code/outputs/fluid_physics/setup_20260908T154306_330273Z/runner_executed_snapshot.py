"""Persistent one-GPU campaign ledger and monotonic, process-group-local watchdog.

An interrupted runner leaves the full reservation charged. It cannot silently
restart/reset a campaign. Only scheduler state is mutable, every attempt is new.
"""
from contextlib import contextmanager
import ctypes
import fcntl
import os
from pathlib import Path
import signal
import subprocess
import time
from .common import output, read_json, write_json, atomic_state, now, resource_snapshot, fingerprint


@contextmanager
def exclusive_lock(path):
    path=output(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('a') as lock:
        try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError as exc: raise RuntimeError('GPU_OR_CAMPAIGN_ALREADY_RUNNING') from exc
        try: yield
        finally: fcntl.flock(lock,fcntl.LOCK_UN)


def group_members(pgid):
    result=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            # comm may contain spaces/parentheses; parse after its final ')'.
            fields=(p/'stat').read_text().rsplit(')',1)[1].split()
            if int(fields[2])==pgid and fields[0]!='Z':result.append(int(p.name))
        except (FileNotFoundError,ProcessLookupError,PermissionError):pass
    return result


def kill_group(pgid,sig):
    if pgid==os.getpgrp():raise RuntimeError('refusing to signal controller group')
    try:os.killpg(pgid,sig)
    except ProcessLookupError:pass


def bounded_process(command,directory,allocation_s,budget,*,monitor_gpu=True):
    """Timer includes creation, MPI/CUDA init, chunks, output and group cleanup."""
    directory=output(directory)
    started=time.monotonic(); deadline=started+allocation_s
    margin=min(float(budget['graceful_margin_s']),allocation_s*.35)
    stop_at=deadline-margin; term_at=deadline-max(2.0,margin/2); kill_at=deadline-1.0
    if allocation_s<2:stop_at=started+allocation_s*.45;term_at=started+allocation_s*.65;kill_at=started+allocation_s*.8
    process=None;timed_out=False;reason='NORMAL_EXIT';exit_code=None; samples=[]; signals=[];remaining=[]
    # Adopt and reap our MPI descendants if their immediate parent exits early.
    ctypes.CDLL(None).prctl(36,1,0,0,0)  # Linux PR_SET_CHILD_SUBREAPER
    next_sample=started
    try:
        with (directory/'console.log').open('xb') as log:
            process=subprocess.Popen(command,cwd=directory,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            atomic_state(directory/'process_identity.json',{'pid':process.pid,'pgid':process.pid,'started_at':now(),'controller_pid':os.getpid()})
            while True:
                current=time.monotonic()
                rc=process.poll()
                if rc is not None:
                    exit_code=rc
                    if not group_members(process.pid):break
                    # A failed mpirun must not leave its compute rank behind.
                    reason='DESCENDANT_CLEANUP';kill_group(process.pid,signal.SIGTERM)
                    if current>=kill_at:kill_group(process.pid,signal.SIGKILL)
                if current>=stop_at and not timed_out:
                    timed_out=True;reason='BUDGET_STOP_REQUESTED'
                    (directory/'STOP_REQUESTED').touch(exist_ok=True)
                    signals.append({'signal':'COOPERATIVE_FILE','elapsed_s':current-started})
                if current>=term_at and not any(s['signal']=='SIGTERM' for s in signals):
                    kill_group(process.pid,signal.SIGTERM);signals.append({'signal':'SIGTERM','elapsed_s':current-started})
                if current>=kill_at:
                    kill_group(process.pid,signal.SIGKILL)
                    if not any(s['signal']=='SIGKILL' for s in signals):signals.append({'signal':'SIGKILL','elapsed_s':current-started})
                # Keep sampling calls away from the stop/kill deadline.
                if monitor_gpu and current>=next_sample and current<stop_at-6:
                    try:samples.append({'elapsed_s':time.monotonic()-started,**resource_snapshot()})
                    except (RuntimeError,subprocess.TimeoutExpired) as exc:samples.append({'elapsed_s':time.monotonic()-started,'error':str(exc)})
                    next_sample=time.monotonic()+budget.get('gpu_sample_interval_s',2)
                if current>=deadline:
                    kill_group(process.pid,signal.SIGKILL)
                    process.wait(timeout=1)
                    break
                time.sleep(min(.1,max(.001,deadline-time.monotonic())))
    except BaseException:
        reason='CONTROLLER_EXCEPTION'
        if process is not None:
            kill_group(process.pid,signal.SIGTERM)
            try:process.wait(timeout=min(1,max(.05,deadline-time.monotonic())))
            except subprocess.TimeoutExpired:kill_group(process.pid,signal.SIGKILL);process.wait(timeout=1)
        raise
    finally:
        if process is not None:
            if process.poll() is None:kill_group(process.pid,signal.SIGKILL);process.wait(timeout=1)
            if group_members(process.pid):kill_group(process.pid,signal.SIGKILL)
            # Child exits already awaited; reap adopted children of this job only.
            try:
                while True:
                    pid,_=os.waitpid(-process.pid,os.WNOHANG)
                    if pid==0:break
            except ChildProcessError:pass
            remaining=group_members(process.pid)
            exit_code=process.returncode
        elapsed=time.monotonic()-started
        good=[x for x in samples if x.get('gpus')]
        peak=max((g['used_MiB'] for x in good for g in x['gpus']),default=None)
        record={'started_at_utc':now(),'elapsed_monotonic_s':elapsed,'allocation_s':allocation_s,
                'exit_code':exit_code,'timeout':timed_out,'reason':reason,'signals':signals,
                'remaining_own_group_processes':remaining,'device_sampled_peak_used_MiB':peak,
                'memory_scope':'sampled device-wide usage including other applications; not exact process peak',
                'resource_samples':samples,'status':'COMPLETED' if exit_code==0 and not timed_out and not remaining else 'STOPPED_OR_FAILED'}
        write_json(directory/'execution.json',record)
    return record


def ledger_read(campaign, campaign_id, budget):
    path=output(campaign)/'budget_ledger.json'
    if path.exists():
        ledger=read_json(path)
        if ledger['campaign_id']!=campaign_id or ledger['campaign_limit_s']!=budget['campaign_limit_s']:
            raise ValueError('CAMPAIGN_ID_OR_BUDGET_MISMATCH; no implicit reset')
        return ledger
    return {'schema_version':1,'campaign_id':campaign_id,'campaign_limit_s':budget['campaign_limit_s'],
            'task_limit_s':budget['task_limit_s'],'created_at':now(),'attempts':[],
            'accounting':'Actual monotonic job wall time including failures; RUNNING/crashed jobs retain full reservation.'}


def used_budget(ledger):
    return sum(a.get('charged_s',a['reserved_s']) for a in ledger['attempts'])


def run_attempt(campaign,c,task_id,cache_identity,create_command,allocation,*,retry_failed=False,monitor_gpu=True):
    campaign=output(campaign);campaign.mkdir(parents=True,exist_ok=True)
    budget=c['budget'];key=fingerprint(cache_identity)
    with exclusive_lock(campaign.parent/'.gpu_exclusive.lock'),exclusive_lock(campaign/'.campaign.lock'):
        ledger=ledger_read(campaign,c['campaign_id'],budget)
        for old in ledger['attempts']:
            if old['status']=='RUNNING':
                raise RuntimeError('UNRECONCILED_RUNNING_RESERVATION: retain charge, check recorded process before explicit reconciliation')
        previous=[a for a in ledger['attempts'] if a['task_id']==task_id]
        if previous and previous[-1]['cache_key']==key and previous[-1]['status']=='COMPLETED':
            d=Path(previous[-1]['directory'])
            # A cache hit also needs the exact immutable output checksums.
            from .common import sha256_file
            hashes=read_json(d/'output_sha256.json')
            if any(not (d/p).is_file() or sha256_file(d/p)!=v for p,v in hashes.items()):
                raise ValueError('CACHE_OUTPUT_HASH_MISMATCH')
            return d,read_json(d/'execution.json'),True
        if previous and not retry_failed:
            raise RuntimeError('EXISTING_DIFFERENT_OR_FAILED_ATTEMPT: explicit --retry-failed required; budget is preserved')
        remaining=budget['campaign_limit_s']-used_budget(ledger)
        allocation=min(float(allocation),float(budget['task_limit_s']),remaining)
        if allocation<budget['minimum_launch_budget_s']:
            return None,{'status':'BUDGET_EXHAUSTED','remaining_s':remaining},False
        if monitor_gpu:
            resources=resource_snapshot()
            if resources['gpus'][0]['free_MiB']<budget['vram_reserve_MiB']+budget['estimated_task_vram_MiB']:
                return None,{'status':'RESOURCE_RESERVE_BLOCKED','resources':resources},False
            if resources['host_memory_bytes']['MemAvailable']<budget['host_reserve_MiB']*1024**2:
                return None,{'status':'HOST_MEMORY_RESERVE_BLOCKED','resources':resources},False
        name=task_id if not previous else task_id+f'_attempt_{len(previous)+1:02d}'
        directory=output(campaign/name);directory.mkdir(exist_ok=False)
        entry={'task_id':task_id,'directory':str(directory),'cache_key':key,'cache_identity':cache_identity,
               'reserved_s':allocation,'charged_s':allocation,'status':'RUNNING','created_at':now()}
        ledger['attempts'].append(entry);atomic_state(campaign/'budget_ledger.json',ledger)
        try:
            command=create_command(directory)
            record=bounded_process(command,directory,allocation,budget,monitor_gpu=monitor_gpu)
        finally:
            # If the controller dies before this, reservation stays charged.
            ep=directory/'execution.json'
            if ep.exists():
                record=read_json(ep);entry['charged_s']=record['elapsed_monotonic_s'];entry['status']=record['status']
                entry['execution_file']=str(ep);atomic_state(campaign/'budget_ledger.json',ledger)
        from .common import sha256_file
        hashes={str(p.relative_to(directory)):sha256_file(p) for p in directory.rglob('*') if p.is_file() and p.name!='output_sha256.json'}
        write_json(directory/'output_sha256.json',hashes)
        return directory,record,False
