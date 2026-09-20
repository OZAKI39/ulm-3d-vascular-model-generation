"""Stage-local native execution with a clean environment and captured provenance."""
import json,os,signal,subprocess,time
from pathlib import Path
from environment_remote import BASE,digest,env
R=BASE/'reports';L=BASE/'logs';B=BASE/'benchmarks'
for p in (R,L,B):p.mkdir(exist_ok=True)
def write(name,d):(R/(name+'.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def load(name):return json.loads((R/(name+'.json')).read_text())
def run(args,name,cwd=BASE,timeout=60,cuda=True,extra_env=None):
    command=[str(BASE/'scripts/use_cuda12_gpu_env.sh')]+list(map(str,args)) if cuda else list(map(str,args))
    runtime_env=env()
    if extra_env:runtime_env.update(extra_env)
    start=time.monotonic();timed_out=False
    log=L/(name+'.log');assert not log.exists(), 'Do not overwrite native execution evidence: '+name
    with log.open('x') as f:
        p=subprocess.Popen(command,cwd=cwd,env=runtime_env,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
        try:p.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timed_out=True;os.killpg(p.pid,signal.SIGKILL);p.wait(timeout=2)
    record={'name':name,'command':command,'cwd':str(cwd),'exit_code':p.returncode,'timeout':timed_out,'timeout_s':timeout,'wall_time_s':time.monotonic()-start,'log':'logs/'+log.name,'sha256':digest(log)}
    write(name,record);print(json.dumps(record),flush=True)
    return record
def okay(d):return d['exit_code']==0 and not d['timeout']
def text(d):return (BASE/d['log']).read_text(errors='replace')
