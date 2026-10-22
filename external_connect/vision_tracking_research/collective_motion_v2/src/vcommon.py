from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get('F10V_DATA_ROOT', ROOT / 'data'))
OUT = Path(os.environ.get('F10V_OUTPUT_ROOT', ROOT))
P = json.loads((ROOT / 'configs/protocol.json').read_text())
NAMES = P['parameter_names']
LIMITS = np.asarray(P['parameter_error_limits'], float)
SCALES = np.asarray(P['parameter_scales'], float)

def clean(x):
    if isinstance(x, dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x, (tuple,list,np.ndarray)): return [clean(v) for v in x]
    if isinstance(x, np.bool_): return bool(x)
    if isinstance(x, np.integer): return int(x)
    if isinstance(x, (float,np.floating)): return float(x) if np.isfinite(x) else None
    if isinstance(x, Path): return str(x)
    return x

def dump(path, obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(clean(obj),ensure_ascii=False,indent=2,allow_nan=False)+'\n');tmp.replace(path)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(2**20),b''): h.update(block)
    return h.hexdigest()

def identity():
    try:
        head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,stderr=subprocess.DEVNULL,text=True).strip()
        branch=subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip()
    except subprocess.CalledProcessError:
        record=json.loads((ROOT/'EXECUTION.json').read_text());head=record['git_head'];branch=record['git_branch']
    return dict(git_head=head,git_branch=branch,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for folder in ['src','scripts','configs'] for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts})

def stats(x):
    a=np.asarray(x,float);a=a[np.isfinite(a)]
    return dict(n=len(a),mean=float(a.mean()),p50=float(np.median(a)),p95=float(np.percentile(a,95)),max=float(a.max())) if len(a) else dict(n=0,mean=None,p50=None,p95=None,max=None)

class ResourceMonitor:
    def __init__(self, interval=1.0):
        import threading,psutil
        self.proc=psutil.Process();self.cpu0=self.proc.cpu_times();self.start=time.perf_counter();self.rows=[]
        self.stop_event=threading.Event();self.interval=interval
        self.thread=threading.Thread(target=self.run,daemon=True);self.thread.start()
    def run(self):
        while not self.stop_event.is_set():
            row={'elapsed_s':time.perf_counter()-self.start,'rss_mb':self.proc.memory_info().rss/2**20}
            try:
                s=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used,power.draw','--format=csv,noheader,nounits'],text=True,timeout=3).splitlines()[0].split(',')
                row.update(gpu_util_percent=float(s[0]),gpu_device_memory_mb=float(s[1]),gpu_power_w=float(s[2]))
            except Exception as e:row['gpu_query_error']=repr(e)
            self.rows.append(row);self.stop_event.wait(self.interval)
    def close(self):
        self.stop_event.set();self.thread.join(4);cpu=self.proc.cpu_times();dt=time.perf_counter()-self.start
        return dict(elapsed_s=dt,cpu_core_equivalents=(cpu.user+cpu.system-self.cpu0.user-self.cpu0.system)/dt,
            samples=self.rows,gpu_scope='device-wide; concurrent processes audited separately')
