import json, os, platform, shutil, subprocess
from pathlib import Path
from datetime import datetime, timezone

def run(args):
    try:
        p=subprocess.run(args,capture_output=True,text=True,timeout=45)
        return dict(returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
    except Exception as e:return dict(error=str(e))

r=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),hostname=platform.node())
r['os_release']=Path('/etc/os-release').read_text()
r['commands']={k:run(v) for k,v in {
    'gpu':['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used,utilization.gpu,compute_cap','--format=csv,noheader'],
    'disk':['df','-B1','/workspace'], 'cpu':['lscpu','-J'],
    'nvcc':['/usr/local/cuda/bin/nvcc','--version'],
    'gcc':['gcc','--version'],'gfortran':['gfortran','--version'],
    'cmake':['cmake','--version'],'python':['/usr/bin/python3','--version'],
    'packages':['dpkg-query','-W','-f=${binary:Package}\t${Version}\n'],
    'libraries':['ldconfig','-p'], 'services':['supervisorctl','status'],
}.items()}
r['cgroup']={k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in ['cpu.max','memory.max','memory.current','memory.events','cpu/cpu.cfs_quota_us','cpu/cpu.cfs_period_us','memory/memory.limit_in_bytes','memory/memory.usage_in_bytes','memory/memory.failcnt','memory/memory.oom_control','cpuset/cpuset.cpus'] if (Path('/sys/fs/cgroup')/k).is_file()}
r['affinity']=sorted(os.sched_getaffinity(0))
r['tools']={k:shutil.which(k) for k in ['rsync','uv','patchelf','mpicc','mpiexec','ninja','cmake','gcc','g++','gfortran']}
r['workspace']=[dict(name=p.name,symlink=p.is_symlink(),target=os.readlink(p) if p.is_symlink() else None) for p in Path('/workspace').iterdir()]
r['cuda_libs']=[str(p) for pattern in ['libcublas.so*','libcusparse.so*','libcudart.so*','libcusolver.so*'] for p in Path('/usr/local/cuda/lib64').glob(pattern)]
caps=run(['vast-capabilities'])
try:
    cap=json.loads(caps['stdout']);r['capabilities']={k:cap.get(k) for k in ['image','hardware']};r['workspace_is_volume']=cap.get('instance',{}).get('workspace_is_volume')
except Exception:r['capabilities_error']='Could not parse capabilities'
print(json.dumps(r,indent=2))
