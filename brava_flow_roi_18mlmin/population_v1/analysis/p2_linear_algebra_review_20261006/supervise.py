from common import *
import subprocess,time
stage=sys.argv[1];assert stage in {'preflight','prepare','build','graph','p1_baseline','p1_asm0','p1_schur','parity','p2_system','finalize','regression','run_pipe','analyze_pipe','output_audit'}
m=managed();m.ROOT=HERE
# Extend the existing supervisor's sampling with GPU and rank RSS evidence.
original=m.cgroup_memory;last=[0,{}]
def sample(*a,**kw):
 r=original(*a,**kw)
 if time.monotonic()-last[0]>4:
  try:
   gpu=subprocess.check_output(['nvidia-smi','--query-gpu=memory.used,memory.free,utilization.gpu','--format=csv,noheader,nounits'],text=True,timeout=4).strip().split(',')
   last[:]=[time.monotonic(),dict(GPU_total_used_MiB=float(gpu[0]),GPU_free_MiB=float(gpu[1]),GPU_utilization_percent=float(gpu[2]))]
  except Exception:pass
 r.update(last[1]);return r
m.cgroup_memory=sample
old_exceeded=m.exceeded
def exceeded(task,sample,initial):
 reason=old_exceeded(task,sample,initial)
 if reason:return reason
 if stage in {'p2_system','run_pipe'} and sample.get('GPU_free_MiB',999999)<2048:return 'GPU_FREE_MEMORY_GUARD'
 return None
m.exceeded=exceeded
name=stage;i=1
while (HERE/'monitoring'/name).exists():i+=1;name=f'{stage}_attempt{i}'
limit=int(Path('/sys/fs/cgroup/memory.max').read_text())/2**30
put('TASK.json',dict(name=name,command=[sys.executable,str(HERE/(stage+'.py')),*sys.argv[2:]],cwd=str(HERE),max_seconds=14400 if stage=='run_pipe' else (1800 if stage.startswith('p1_') or stage=='p2_system' else 1200),max_tree_RSS_GiB=24,max_nonreclaimable_GiB=.8*limit,environment={'PYTHONUNBUFFERED':'1','PYTHONDONTWRITEBYTECODE':'1'}))
try:m.main()
finally:(HERE/'TASK.json').unlink(missing_ok=True)
r=read(HERE/'monitoring'/name/'EXECUTION.json');print(name,r['exit_code'],r['stop_reason'],r['elapsed_s'])
if r['exit_code']!=0 or r['stop_reason']:sys.exit(1)
