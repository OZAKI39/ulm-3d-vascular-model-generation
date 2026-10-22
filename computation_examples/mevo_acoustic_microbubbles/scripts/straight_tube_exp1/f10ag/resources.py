"""Bounded read-only runtime sampling for only this campaign."""
import time,subprocess,json,os
from .common import *
def main():
 start=time.perf_counter();dest=BASE/'RESOURCE_SAMPLES.jsonl';samples=0
 with dest.open('x') as out:
  while time.perf_counter()-start<6000 and not (BASE/'RESOURCE_STOP').exists():
   r=dict(elapsed_s=time.perf_counter()-start,loadavg=os.getloadavg(),gpu=subprocess.check_output(['nvidia-smi','--query-gpu=timestamp,utilization.gpu,utilization.memory,memory.used,memory.total,power.draw','--format=csv,noheader,nounits'],text=True).strip(),identity=identity())
   try:
    import psutil
    r.update(cpu_percent=psutil.cpu_percent(),memory=psutil.virtual_memory()._asdict())
   except ImportError:pass
   out.write(json.dumps(r)+'\n');out.flush();samples+=1;time.sleep(10)
 write(BASE/'RESOURCE_SUMMARY.json',dict(samples=samples,wall_s=time.perf_counter()-start,scope='Read-only campaign monitoring, shared hardware; not scientific throughput denominator',identity=identity()))
if __name__=='__main__':main()
