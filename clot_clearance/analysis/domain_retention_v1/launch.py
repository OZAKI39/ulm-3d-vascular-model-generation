"""External process timing survives a Python/native child failure."""
from common import *
import sys,subprocess,time,resource
case=sys.argv[1];assert case in ['preflight','full']
target=A/('preflight.py' if case=='preflight' else 'run_protocol.py')
folder=O/case
with (folder/'PROCESS_STARTED.json').open('x') as f:json.dump(dict(command=[sys.executable,str(target)],parent_pid=os.getpid()),f)
start=time.perf_counter();before=resource.getrusage(resource.RUSAGE_CHILDREN);code=None
try:
 with (folder/'CONSOLE.log').open('w') as log:
  child=subprocess.Popen([sys.executable,'-u',str(target)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  for line in child.stdout:log.write(line);log.flush();print(line,end='',flush=True)
  code=child.wait()
finally:
 after=resource.getrusage(resource.RUSAGE_CHILDREN)
 dump(folder/'PROCESS_RESOURCES.json',dict(exit_code=code,wall_s=time.perf_counter()-start,user_CPU_s=after.ru_utime-before.ru_utime,system_CPU_s=after.ru_stime-before.ru_stime,child_peak_RSS_bytes=after.ru_maxrss*1024,source='external parent wait4 resource accounting',MPI_ranks=1,threads=1))
sys.exit(code or 0)
