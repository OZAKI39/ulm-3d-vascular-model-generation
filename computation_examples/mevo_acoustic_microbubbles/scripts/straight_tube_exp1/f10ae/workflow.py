"""Bounded phases; original complete/failed evidence retained on resumption."""
import time,threading,subprocess,os,json,sys
from .common import *
def monitor(stop):
    with (BASE/'RESOURCES.jsonl').open('a') as out:
        while not stop.is_set():
            try:
                r=dict(utc=time.time(),load=os.getloadavg(),cpu_count=os.cpu_count(),gpu=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,power.draw','--format=csv,noheader,nounits'],text=True).strip())
                out.write(json.dumps(r)+'\n');out.flush()
            except Exception:pass
            stop.wait(3)
def main():
    verify();start=time.perf_counter();stop=threading.Event();thread=threading.Thread(target=monitor,args=(stop,),daemon=True);thread.start()
    if not (BASE/'WORKFLOW_START.json').exists():write(BASE/'WORKFLOW_START.json',dict(identity=identity(),max_wall_s=5400))
    try:
        subprocess.run([sys.executable,'-m','pytest','-q','tests/f10ae'],check=True)
        from . import baseline,search,explore,validate
        for name,module,marker in [('baseline',baseline,'BASELINE_VALIDATION.json'),('search',search,'OPTIMIZER_COMPARISON.json'),('explore',explore,'SELECTED_FROZEN.json'),('validate',validate,'VALIDATION_COMPLETE.json')]:
            if (BASE/marker).exists():continue
            if time.perf_counter()-start>5400:raise RuntimeError('Overall scientific budget reached')
            write(BASE/'WORKFLOW_PROGRESS.json',dict(phase=name,identity=identity(),elapsed_s=time.perf_counter()-start))
            module.main()
        write(BASE/'WORKFLOW_COMPLETE.json',dict(wall_s=time.perf_counter()-start,identity=identity()))
    finally:stop.set();thread.join(5)
if __name__=='__main__':main()
