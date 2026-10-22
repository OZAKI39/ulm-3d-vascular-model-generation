from common import *
import subprocess,time
stage=sys.argv[1];assert stage in {'preflight','prepare','build','graph','p1_baseline','p1_asm0','p1_schur','p1_candidate','build_guard','graph_bound','parity','p2_system','finalize','regression','run_pipe','analyze_pipe','output_audit'}
m=managed();m.ROOT=HERE
# Extend the existing supervisor's sampling with GPU and rank RSS evidence.
original=m.cgroup_memory
# NVML avoids launching nvidia-smi per sample; one second host/GPU observations.
import ctypes
class Mem(ctypes.Structure):
 _fields_=[('total',ctypes.c_ulonglong),('free',ctypes.c_ulonglong),('used',ctypes.c_ulonglong)]
nvml=ctypes.CDLL('libnvidia-ml.so.1');assert nvml.nvmlInit_v2()==0
handle=ctypes.c_void_p();assert nvml.nvmlDeviceGetHandleByIndex_v2(0,ctypes.byref(handle))==0
def sample(*a,**kw):
 r=original(*a,**kw);v=Mem();code=nvml.nvmlDeviceGetMemoryInfo(handle,ctypes.byref(v))
 if code:raise RuntimeError('NVML monitor unavailable')
 r.update(GPU_total_used_MiB=v.used/2**20,GPU_free_MiB=v.free/2**20)
 return r
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
