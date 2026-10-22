import sys,os,time,json,subprocess,platform,importlib.metadata as md
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import torch,psutil
from vcommon import OUT,dump,identity

def command(c):
 try:return subprocess.check_output(c,text=True,stderr=subprocess.STDOUT,timeout=15)
 except Exception as e:return repr(e)
start=time.perf_counter();f=OUT/'evidence/disk_probe.tmp';f.parent.mkdir(parents=True,exist_ok=True)
block=os.urandom(2**20)
with f.open('wb') as stream:
 for _ in range(64):stream.write(block)
 stream.flush();os.fsync(stream.fileno())
write_s=time.perf_counter()-start;t=time.perf_counter();b=f.read_bytes();read_s=time.perf_counter()-t;f.unlink()
record=dict(identity=identity(),platform=platform.platform(),python=sys.version,cpu_physical=psutil.cpu_count(False),cpu_logical=psutil.cpu_count(),memory=psutil.virtual_memory()._asdict(),disk=psutil.disk_usage(OUT)._asdict(),disk_probe=dict(bytes=len(b),write_fsync_s=write_s,warm_read_s=read_s,caveat='64 MiB bounded probe; warm read, not sustained storage benchmark'),lscpu=command(['lscpu']),nvidia_smi=command(['nvidia-smi']),gpu_processes=command(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv']),torch=torch.__version__,cuda_runtime=torch.version.cuda,bf16=torch.cuda.is_bf16_supported(),devices=[torch.cuda.get_device_properties(i).__dict__ if hasattr(torch.cuda.get_device_properties(i),'__dict__') else str(torch.cuda.get_device_properties(i)) for i in range(torch.cuda.device_count())],packages={d.metadata['Name']:d.version for d in md.distributions()})
x=torch.randn(128,128,device='cuda');y=x@x;torch.cuda.synchronize();record['cuda_operation_pass']=bool(torch.isfinite(y).all());dump(OUT/'evidence/environment.json',record);print(json.dumps({k:record[k] for k in ['torch','cuda_runtime','bf16','cuda_operation_pass']}))
