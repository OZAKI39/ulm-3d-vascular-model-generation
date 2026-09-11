"""Read-only server probe. Executed over SSH stdin; emits one JSON record."""
from pathlib import Path
import subprocess,shutil,os,json,hashlib,stat

def run(argv):
    try:
        r=subprocess.run(argv,capture_output=True,text=True,timeout=20)
        return dict(argv=argv,exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr)
    except subprocess.TimeoutExpired:return dict(argv=argv,exit_code=124,stdout='',stderr='TIMEOUT')
def main():
    workspace=Path('/workspace/bloodflow');commands={}
    probes={'gpu':['nvidia-smi','--query-gpu=name,uuid,driver_version,memory.total,memory.used,utilization.gpu','--format=csv,noheader'],
        'gpu_processes':['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv,noheader'],
        'nvcc':['/usr/local/cuda/bin/nvcc','--version'],'architectures':['/usr/local/cuda/bin/nvcc','--list-gpu-code'],
        'gcc':['gcc','-dumpfullversion'],'gxx':['g++','-dumpfullversion'],'cmake':['cmake','--version'],'mpi':['mpirun','--version'],
        'hdf5':['h5pcc','-showconfig'],'hdf5_link':['h5pcc','-show'], 'df':['df','-h','/workspace'],
        'packages':[str(workspace/'.venv/bin/python'),'-B','-m','pip','list','--format=json'],
        'python':[str(workspace/'.venv/bin/python'),'-B','-c','import sys,sysconfig,json,os;print(json.dumps(dict(executable=sys.executable,version=sys.version,prefix=sys.prefix,base_prefix=sys.base_prefix,include=sysconfig.get_path("include"),python_h_exists=os.path.isfile(sysconfig.get_path("include")+"/Python.h"))))'],
        'processes':['ps','-eo','pid,ppid,comm'],'sizes':['du','-sb',str(workspace/'uploads/20260910T212439Z'),str(workspace/'.venv')]}
    for name,argv in probes.items():commands[name]=run(argv)
    disk=shutil.disk_usage('/workspace')
    return dict(recorded_at=__import__('datetime').datetime.now(__import__('datetime').timezone.utc).isoformat(),uid=os.getuid(),commands=commands,
        tools={x:shutil.which(x) for x in ['nvcc','gcc','g++','cmake','mpicc','mpicxx','mpirun','h5pcc','rsync','cuobjdump','tmux']},
        cgroup={x:(Path('/sys/fs/cgroup')/x).read_text().strip() for x in ['memory.max','memory.current','cpu.max','cpuset.cpus.effective']},
        disk=dict(total=disk.total,used=disk.used,free=disk.free),
        existing_builds=[p.name for p in (workspace/'build').iterdir()] if (workspace/'build').is_dir() else [],
        image_environment=dict(path='/venv/main',resolved=str(Path('/venv/main').resolve()),policy='Never write or install here'),
        python_environment=dict(path=str(workspace/'.venv'),exists=(workspace/'.venv').is_dir(),symlink=(workspace/'.venv').is_symlink()))
print(json.dumps(main()))
