#!/usr/bin/env python3
import os,sys,time,json,subprocess,hashlib,fcntl,traceback,shutil
from pathlib import Path
W=Path(__file__).resolve().parents[1];A=Path('/workspace/hemocell_restore/archive/minimal_restore_20260915_143631');T=Path('/workspace/hemocell_restore/toolchains/nvhpc_26_5');phase='START'
def save(name,obj):
 p=W/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(p)
def run(name,cmd,cwd=None,env=None):
 global phase
 phase=name;save('PREPARATION_STATE.json',dict(state='RUNNING',phase=phase,pid=os.getpid(),unix=time.time()));print('START',name,flush=True)
 t=time.monotonic()
 with (W/'logs'/f'{name}.log').open('w') as f:p=subprocess.run(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT)
 save('provenance/'+name+'.json',dict(command=cmd,cwd=str(cwd),returncode=p.returncode,seconds=time.monotonic()-t))
 if p.returncode:raise RuntimeError(name+' failed '+str(p.returncode))
 print('PASS',name,flush=True)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
lock=(W/'preparation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
with (W/'PREPARATION_STARTED.json').open('x') as f:json.dump(dict(pid=os.getpid(),unix=time.time()),f)
try:
 run('archive_preverify',['python3','-B',str(A/'reports/verify_bundle.py'),str(A),'--readonly'])
 run('materialize_work',['python3','-B',str(W/'scripts/prepare_work.py'),str(A),str(W)])
 D=W/'toolchain_env/download';D.mkdir();X=W/'toolchain_env/installer';X.mkdir()
 package=D/'nvhpc_2026_265_Linux_x86_64_cuda_13.2.tar.gz'
 url='https://developer.download.nvidia.com/hpc-sdk/26.5/'+package.name
 run('nvhpc_download',['curl','--fail','--location','--connect-timeout','30','--output',str(package),url])
 assert package.stat().st_size==6851790769 and sha(package)=='8849b3909f7e75df4c32af4dc585bf0701d85495da43221b090a9af77df203b6','NVHPC SHA mismatch'
 run('nvhpc_md5_source',['curl','--fail','--location','--output',str(D/'official_md5sum.txt'),'https://developer.download.nvidia.com/hpc-sdk/26.5/md5sum.txt'])
 assert 'db12c569f852fd4789a1e3c20011a1db' in (D/'official_md5sum.txt').read_text()
 with package.open('rb') as f:md5=hashlib.file_digest(f,'md5').hexdigest()
 assert md5=='db12c569f852fd4789a1e3c20011a1db'
 save('provenance/NVHPC_DOWNLOAD_VERIFIED.json',dict(status='PASS',url=url,bytes=package.stat().st_size,sha256=sha(package),md5=md5,prefix=str(T)))
 run('nvhpc_extract',['tar','-xzf',str(package),'-C',str(X)])
 dirs=list(X.iterdir());assert len(dirs)==1 and (dirs[0]/'install').is_file()
 assert not T.exists(),'Unique exact toolchain prefix must not be overwritten'
 env=dict(os.environ,NVHPC_SILENT='true',NVHPC_INSTALL_DIR=str(T),NVHPC_INSTALL_TYPE='single',NVHPC_DEFAULT_CUDA='13.2');env.pop('NVHPC_STDPAR_CUDACC',None)
 run('nvhpc_install',[str(dirs[0]/'install')],cwd=dirs[0],env=env)
 nvc=T/'Linux_x86_64/26.5/compilers/bin/nvc++';assert nvc.is_file()
 run('nvc_version',[str(nvc),'--version'])
 missing=[]
 for pkg in ['openmpi-bin','libopenmpi-dev']:
  p=subprocess.run(['dpkg-query','-W','-f=${Status}',pkg],capture_output=True,text=True)
  if p.returncode or 'install ok installed' not in p.stdout:missing.append(pkg)
 if missing:
  run('apt_index',['apt-get','update'])
  run('minimal_mpi_install',['apt-get','install','-y','--no-install-recommends',*missing],env=dict(os.environ,DEBIAN_FRONTEND='noninteractive'))
 save('PREPARATION_STATE.json',dict(state='PASS',phase='READY_FOR_BUILD',pid=os.getpid(),unix=time.time(),minimal_packages_requested=missing,solver_runs=0))
except BaseException as e:
 traceback.print_exc();save('PREPARATION_STATE.json',dict(state='FAIL',phase=phase,error=repr(e),traceback=traceback.format_exc(),pid=os.getpid(),unix=time.time(),solver_runs=0));sys.exit(2)
