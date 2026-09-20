"""Read-only environment fingerprint, including every existing CUDA13 file."""
import hashlib,json,os,subprocess,time
from pathlib import Path
BASE=Path(__file__).resolve().parent
def digest(p,algorithm='sha256'):
    h=hashlib.new(algorithm)
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
def inventory(root):
    root=Path(root);files={}
    for base,dirs,names in os.walk(root,followlinks=False):
        for name in dirs+names:
            p=Path(base)/name;s=p.lstat();k=str(p.relative_to(root))
            if p.is_symlink():files[k]={'symlink':os.readlink(p),'mode':s.st_mode}
            elif p.is_file():files[k]={'sha256':digest(p),'size':s.st_size,'mtime_ns':s.st_mtime_ns,'mode':s.st_mode}
            elif p.is_dir():files[k]={'directory':True,'mode':s.st_mode}
    return files
def env():
    d={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
    return dict(d,PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0')
def probe(args):
    start=time.monotonic();p=subprocess.run(args,env=env(),capture_output=True,text=True,timeout=20)
    return {'command':args,'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'wall_time_s':time.monotonic()-start}
def snapshot():
    link=Path('/usr/local/cuda');target=link.resolve()
    paths=set([Path('/usr/bin/nvidia-smi')])
    for folder in ('/usr/lib/x86_64-linux-gnu','/usr/local/nvidia/lib64'):
        for pattern in ('libcuda.so*','libnvidia*.so*'):
            paths.update(Path(folder).glob(pattern))
    drivers={str(p):{'realpath':str(p.resolve()),'sha256':digest(p)} for p in paths if p.is_file()}
    cmds=[['hostname'],['uname','-a'],['cat','/etc/os-release'],['gcc','--version'],['g++','--version'],['nvidia-smi'],['nvidia-smi','--query-gpu=name,driver_version,memory.total,compute_cap','--format=csv,noheader'],[str(target/'bin/nvcc'),'--version'],['dpkg-query','-W','-f=${Package}\t${Version}\t${Status}\n','*nvidia*','*cuda*'],['df','-h',str(BASE),'/tmp']]
    extra_cuda=BASE.parent/'sv1_3h/external/cuda-12.6'
    compilers={}
    for n in ('gcc','g++'):
        p=Path('/usr/bin')/n;compilers[n]={'path':str(p),'realpath':str(p.resolve()),'symlink':os.readlink(p) if p.is_symlink() else None,'sha256':digest(p)}
    historical={name:inventory(BASE.parent/name) for name in ('sv1_3j','sv1_3l','sv1_3m','sv1_3n','sv1_3g/external/gpu_mpi')}
    return {'historical_gpu_stacks':historical,'cuda126_prefix':str(extra_cuda),'cuda126_files':inventory(extra_cuda),'default_compilers':compilers,'timestamp':time.time(),'selected_environment':{k:os.environ.get(k) for k in ('PATH','LD_LIBRARY_PATH','CUDA_HOME','CUDA_VISIBLE_DEVICES')},'cuda_default_link':os.readlink(link) if link.is_symlink() else None,'cuda13_prefix':str(target),'cuda13_files':inventory(target),'driver_files':drivers,'ld_configuration':inventory('/etc/ld.so.conf.d'),'probes':[probe(a) for a in cmds]}
if __name__=='__main__':
    (BASE/'reports').mkdir(exist_ok=True)
    p=BASE/'reports/pre_install_environment.json';assert not p.exists()
    p.write_text(json.dumps(snapshot(),indent=2));print('Pre-install CUDA13 and driver fingerprint saved',flush=True)
