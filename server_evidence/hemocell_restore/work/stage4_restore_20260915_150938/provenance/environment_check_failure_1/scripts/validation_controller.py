#!/usr/bin/env python3
import os,sys,time,json,subprocess,hashlib,fcntl,traceback,shutil,socket,csv
from pathlib import Path
W=Path(__file__).resolve().parents[1];A=Path('/workspace/hemocell_restore/archive/minimal_restore_20260915_143631');T=Path('/workspace/hemocell_restore/toolchains/nvhpc_26_5/Linux_x86_64/26.5');phase='START';attempts=[]
def save(p,obj):
 p=W/p;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(obj,indent=2)+'\n');t.replace(p)
def capture(cmd):
 p=subprocess.run(cmd,capture_output=True,text=True);return dict(command=cmd,returncode=p.returncode,stdout=p.stdout,stderr=p.stderr)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(2**20),b''):h.update(b)
 return h.hexdigest()
def command(name,cmd,cwd=None):
 global phase
 phase=name;save('VALIDATION_STATE.json',dict(state='RUNNING',phase=phase,pid=os.getpid(),unix=time.time(),run_id=W.name));print('START',name,flush=True)
 t=time.monotonic()
 with (W/'logs'/f'{name}.log').open('w') as f:p=subprocess.run(cmd,cwd=cwd,stdout=f,stderr=subprocess.STDOUT)
 attempts.append(dict(stage=name,command=cmd,cwd=str(cwd),returncode=p.returncode,seconds=time.monotonic()-t,log='logs/'+name+'.log'))
 save('provenance/EXECUTION_COMMANDS.json',attempts)
 if p.returncode:raise RuntimeError(name+' failed '+str(p.returncode))
 print('PASS',name,flush=True)
def binary_receipt(exe,name):
 d=dict(path=str(exe),sha256=sha(exe),ldd=capture(['ldd',str(exe)]),old_binary_executed=False,host=socket.gethostname(),compiler='NVHPC26.5',gpu_arch='cc89')
 assert d['ldd']['returncode']==0 and 'not found' not in d['ldd']['stdout'],d
 save('provenance/'+name+'.json',d)
lock=(W/'validation.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
with (W/'VALIDATION_STARTED.json').open('x') as f:json.dump(dict(pid=os.getpid(),unix=time.time(),automatic_solver_retries=0),f)
try:
 assert json.loads((W/'PREPARATION_STATE.json').read_text())['state']=='PASS'
 for d in ['tmp','tmp/cache','tmp/cuda-cache','scripts/recipe']:(W/d).mkdir(parents=True,exist_ok=True)
 lib=W/'tbb/usr/lib/x86_64-linux-gnu'
 bundled_nsys=T/'profilers/13.2/Nsight_Systems/bin'
 profiler_path=(str(bundled_nsys)+':') if (bundled_nsys/'nsys').is_file() else ''
 os.environ.update(PATH=str(T/'compilers/bin')+':'+str(T/'cuda/13.2/bin')+':'+profiler_path+'/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin',LD_LIBRARY_PATH=str(lib)+':'+str(T/'compilers/lib')+':'+str(T/'cuda/13.2/lib64'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',HWLOC_COMPONENTS='-opencl',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1',TMPDIR=str(W/'tmp'),CUDA_CACHE_PATH=str(W/'tmp/cuda-cache'),XDG_CACHE_HOME=str(W/'tmp/cache'))
 environment={k:os.environ[k] for k in ['PATH','LD_LIBRARY_PATH','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','HWLOC_COMPONENTS','CUDA_VISIBLE_DEVICES','TMPDIR','CUDA_CACHE_PATH','XDG_CACHE_HOME']};save('toolchain_env/ENV.json',environment)
 commands={'nvc':['nvc++','--version'],'nvcc':['nvcc','--version'],'mpi':['/usr/bin/mpirun','--version'],'cmake':['cmake','--version'],'gcc':['g++','--version'],'make':['make','--version'],'python':['python3','--version'],'numpy':['python3','-c','import numpy;print(numpy.__version__)'],'nsys':['nsys','--version'],'gpu_arch':['nvidia-smi','--query-gpu=name,compute_cap,driver_version','--format=csv,noheader,nounits'],'compiler_arch_help':['nvc++','-help','-gpu']}
 deps={k:capture(c) for k,c in commands.items()};save('toolchain_env/DEPENDENCY_COMMANDS.json',deps)
 assert all(x['returncode']==0 for x in deps.values()),'Dependency command failed'
 assert '26.5' in deps['nvc']['stdout']+deps['nvc']['stderr']
 assert 'release 13.2' in deps['nvcc']['stdout']
 assert 'RTX 4090' in deps['gpu_arch']['stdout'] and '8.9' in deps['gpu_arch']['stdout']
 assert 'cc89' in deps['compiler_arch_help']['stdout']+deps['compiler_arch_help']['stderr']
 with (W/'SYSTEM_DEPENDENCIES.tsv').open('w',newline='') as f:
  out=csv.writer(f,delimiter='\t');out.writerow(['dependency','status','executable','version_evidence'])
  for k,v in deps.items():out.writerow([k,'PASS',shutil.which(commands[k][0]),(v['stdout']+v['stderr']).strip().replace('\n',' | ')])
  out.writerow(['TBB','PASS',str(lib), 'Archived exact dependency; hashes in WORKING_COPY_SHA256.json'])
 save('toolchain_env/NVHPC_READY.json',dict(status='PASS',version='26.5',cuda_component='13.2',prefix=str(T),arch='cc89',actual_compute_capability='8.9',driver_changed=False,global_CUDA_changed=False))
 # Verify actual physical-core binding without running any numerical solver.
 probe="import os,json,pathlib; a=sorted(os.sched_getaffinity(0)); c=sorted(set((pathlib.Path(f'/sys/devices/system/cpu/cpu{x}/topology/physical_package_id').read_text().strip(),pathlib.Path(f'/sys/devices/system/cpu/cpu{x}/topology/core_id').read_text().strip()) for x in a)); print(json.dumps(dict(affinity=a,physical_cores=c))); assert len(c)==1"
 command('mpi_binding_preflight',['/usr/bin/mpirun','--allow-run-as-root','-np','1','--bind-to','core','--map-by','core','--report-bindings','/usr/bin/python3','-c',probe])
 save('provenance/SCRIPT_EXECUTION_SHA256.json',{str(p.relative_to(W)):sha(p) for p in sorted((W/'scripts').rglob('*')) if p.is_file()})
 linker='-L'+str(lib);tbb=lib/'cmake/TBB';B=W/'build'
 command('official_configure',['cmake','-S',str(W/'official_palabos/examples/gpuExamples/cavity3d'),'-B',str(B/'official'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_C_COMPILER='+str(T/'compilers/bin/nvc'),'-DCMAKE_CXX_COMPILER='+str(T/'compilers/bin/nvc++'),'-DCMAKE_CXX_COMPILER_ARG1=-gpu=cc89','-DTBB_DIR='+str(tbb),'-DENABLE_MPI=ON','-DCMAKE_EXE_LINKER_FLAGS='+linker])
 command('official_build',['cmake','--build',str(B/'official'),'--parallel','8','--verbose'])
 binary_receipt(B/'cavity3d','OFFICIAL_BINARY')
 command('official_smoke',['python3','-u','-B',str(W/'scripts/run_official.py')])
 assert json.loads((W/'OFFICIAL_GPU_SMOKE.json').read_text())['status']=='PASS'
 command('stage4_configure',['cmake','-S',str(W/'scripts/recipe'),'-B',str(B/'stage4'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_C_COMPILER='+str(T/'compilers/bin/nvc'),'-DCMAKE_CXX_COMPILER='+str(T/'compilers/bin/nvc++'),'-DTBB_DIR='+str(tbb),'-DRESTORE_GPU_ARCH=89','-DRESTORE_GPU_BUILD_ROOT='+str(B),'-DCMAKE_EXE_LINKER_FLAGS='+linker])
 command('stage4_build',['cmake','--build',str(B/'stage4'),'--parallel','1','--verbose'])
 binary_receipt(B/'stage4/vascularPoC','STAGE4_BINARY')
 save('BUILD_PROVENANCE.json',dict(status='PASS',compiler=deps['nvc'],cuda=deps['nvcc'],MPI=deps['mpi'],cmake=deps['cmake'],gpu_architecture='cc89',actual_compute_capability='8.9',compiler_flags=(B/'stage4/CMakeFiles/vascularPoC.dir/flags.make').read_text(),link_command=(B/'stage4/CMakeFiles/vascularPoC.dir/link.txt').read_text(),source_sha256=sha(W/'source/vascularPoC.cpp'),binary=json.loads((W/'provenance/STAGE4_BINARY.json').read_text()),palabos_commit='4127697e90169bbef982295f1d1c933cf6e90caa',source_diff='NONE',platform_adaptation_only=True,old_RTX5090_binary='PROVENANCE_ONLY',recipes_sha256={str(p.relative_to(W)):sha(p) for p in [W/'scripts/recipe/CMakeLists.txt',W/'source/CMakeLists.txt',W/'official_palabos/examples/gpuExamples/cavity3d/CMakeLists.txt']},attempts=attempts))
 for n in (200,1000,5000):
  command(f'stage4_{n}_run',['python3','-u','-B',str(W/'scripts/run_stage4.py'),str(n)])
  command(f'stage4_{n}_compare',['python3','-u','-B',str(W/'scripts/compare_restore.py'),str(n)])
 result=dict(state='PASS',phase='ALL_REQUESTED_RUNS_TERMINAL',unix=time.time(),pid=os.getpid(),run_id=W.name,automatic_solver_retries=0)
except BaseException as e:
 traceback.print_exc();result=dict(state='FAIL',phase=phase,error=repr(e),traceback=traceback.format_exc(),unix=time.time(),pid=os.getpid(),run_id=W.name,automatic_solver_retries=0)
save('VALIDATION_STATE.json',result);save('provenance/VALIDATION_TERMINAL.json',result)
print('VALIDATION_TERMINAL',json.dumps(result),flush=True)
sys.exit(0 if result['state']=='PASS' else 2)
