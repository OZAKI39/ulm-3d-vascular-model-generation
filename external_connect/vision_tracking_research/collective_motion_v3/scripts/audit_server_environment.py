"""Small isolated filesystem probe and current resource/environment inventory."""
import subprocess,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
script=r'''
import os,time,json,tempfile,subprocess,platform,hashlib
from pathlib import Path
import torch,psutil
root=Path('/workspace/f10w_vision_20261017/setup');blob=os.urandom(32*1024*1024)
fd,p=tempfile.mkstemp(prefix='f10w_io_probe_',dir=root);t=time.perf_counter()
with os.fdopen(fd,'wb') as f:f.write(blob);f.flush();os.fsync(f.fileno())
write_s=time.perf_counter()-t;t=time.perf_counter();read=Path(p).read_bytes();read_s=time.perf_counter()-t;assert read==blob;Path(p).unlink()
r=dict(utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),cpu=platform.processor(),logical_cores=psutil.cpu_count(),physical_cores=psutil.cpu_count(logical=False),memory=psutil.virtual_memory()._asdict(),disk=psutil.disk_usage(str(root))._asdict(),gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used,utilization.gpu','--format=csv'],text=True),gpu_processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv'],text=True),torch=torch.__version__,cuda=torch.version.cuda,cuda_available=torch.cuda.is_available(),io=dict(bytes=len(blob),write_fsync_s=write_s,read_s=read_s,write_mib_s=32/write_s,read_mib_s=32/read_s,scope='32 MiB bounded file; read may be page-cache; not sustained physical-disk benchmark'),supervisor_f10w=subprocess.run(['supervisorctl','status'],capture_output=True,text=True).stdout.splitlines())
r['supervisor_f10w']=[x for x in r['supervisor_f10w'] if x.startswith('f10w_')]
print(json.dumps(r))
'''
cmd=['ssh','-i','/home/lzy/.ssh/vast_step3b_ed25519','-o','BatchMode=yes','-p','15283','root@95.3.33.46','/workspace/f10w_vision_20261017/.venv/bin/python -']
r=subprocess.run(cmd,input=script,text=True,capture_output=True);print(r.stderr,file=__import__('sys').stderr);r.check_returncode();d=json.loads(r.stdout);(ROOT/'evidence/server_environment_final.json').write_text(json.dumps(d,indent=2)+'\n');print(d['gpu'],d['io'])
