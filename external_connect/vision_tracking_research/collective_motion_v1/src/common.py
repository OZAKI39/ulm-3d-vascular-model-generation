from pathlib import Path
import hashlib, json, subprocess, time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
P = json.loads((ROOT/'configs/protocol.json').read_text())
NAMES = ['translation_vx','translation_vy','divergence','vorticity','shear1','shear2']
VIDEO = Path('/home/lzy/projects/external_connect/sample_videos/Recording 2026-10-15 185901 截取视频.mp4')
SIM = Path('/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10s_20261015/computation_examples/mevo_acoustic_microbubbles')

def clean(x):
    if isinstance(x, dict): return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x, (tuple,list,np.ndarray)): return [clean(v) for v in x]
    if isinstance(x, (np.bool_,)): return bool(x)
    if isinstance(x, (np.integer,)): return int(x)
    if isinstance(x, (float,np.floating)): return float(x) if np.isfinite(x) else None
    return x

def dump(path, obj):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(clean(obj),indent=2,allow_nan=False)+'\n')

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''): h.update(chunk)
    return h.hexdigest()

def stats(a):
    a=np.array(a,float); a=a[np.isfinite(a)]
    return dict(zip(['n','p50','p95','max'],[len(a),*np.percentile(a,[50,95,100])])) if len(a) else None

def identity():
    return {'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            'git_branch':subprocess.check_output(['git','branch','--show-current'],cwd=ROOT,text=True).strip(),
            'source_sha256':{str(p.relative_to(ROOT)):sha(p) for d in ['src','scripts','configs'] for p in (ROOT/d).rglob('*') if p.is_file() and '__pycache__' not in str(p)},
            'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}

class ResourceMonitor:
    """Read-only process/GPU samples. GPU utilization includes other desktop users."""
    def __init__(self):
        import threading,psutil
        self.process=psutil.Process();self.samples=[];self.stop_event=threading.Event();self.start_time=time.monotonic();self.cpu0=self.process.cpu_times()
        self.thread=threading.Thread(target=self._run,daemon=True);self.thread.start()
    def _run(self):
        while not self.stop_event.is_set():
            row={'t_s':time.monotonic()-self.start_time,'rss_mb':self.process.memory_info().rss/2**20}
            try:
                s=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True,timeout=2).splitlines()[0].split(',')
                row.update(gpu_device_utilization_percent=float(s[0]),gpu_device_memory_mb=float(s[1]))
            except Exception as e:row['gpu_error']=str(e)
            self.samples.append(row);self.stop_event.wait(.2)
    def close(self):
        self.stop_event.set();self.thread.join(3);cpu=self.process.cpu_times();dt=time.monotonic()-self.start_time
        return {'elapsed_s':dt,'cpu_core_equivalents':((cpu.user+cpu.system)-(self.cpu0.user+self.cpu0.system))/dt,'samples':self.samples,'gpu_scope':'device-wide; desktop graphics are not isolated'}
