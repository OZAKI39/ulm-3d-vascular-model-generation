from pathlib import Path
import os,sys,json,hashlib,subprocess,time,fcntl,traceback,math
N=Path(__file__).resolve().parents[1];B=Path('/workspace/hemocell_restore/work/stage4_restore_20260915_150938');phase='STATIC_CHECK'
def save(name,d):
 p=N/name;t=p.with_suffix(p.suffix+'.tmp');t.write_text(json.dumps(d,indent=2)+'\n');t.replace(p)
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def command(name,cmd):
 global phase
 phase=name;save('BUILD_STATE.json',dict(state='RUNNING',phase=phase,pid=os.getpid(),unix=time.time()))
 with (N/'logs'/f'{name}.log').open('w') as f:p=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT)
 save('provenance/command_'+name+'.json',dict(command=cmd,returncode=p.returncode))
 if p.returncode:raise RuntimeError(name+' failed '+str(p.returncode))
 print('PASS',name,flush=True)
lock=(N/'build.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
with (N/'BUILD_STARTED.json').open('x') as f:json.dump(dict(pid=os.getpid(),unix=time.time()),f)
try:
 env=json.loads((B/'toolchain_env/ENV.json').read_text());env.update(TMPDIR=str(N/'tmp'),CUDA_CACHE_PATH=str(N/'tmp/cuda-cache'),XDG_CACHE_HOME=str(N/'tmp/cache'))
 for p in ['tmp/cuda-cache','tmp/cache']:(N/p).mkdir(parents=True,exist_ok=True)
 os.environ.update(env);save('provenance/EXECUTION_ENV.json',env)
 c=json.loads((N/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json').read_text());x=json.loads((N/'inputs/NUMERICS_INPUT.json').read_text());g=json.loads((N/'provenance/NUMERICS_GENERATION_RECEIPT.json').read_text())
 assert sha(N/'contracts/NEW_MEDIUM_NUMERICS_CONTRACT.json')==g['output_sha256'] and sha(N/'scripts/prepare_numerics.py')==c['script_sha256'] and sha(N/'inputs/NUMERICS_INPUT.json')==c['input_sha256']
 expected={'dt_s':6.65994708146172e-9,'outlet_01_rho_lu':1.0000484343922278,'outlet_02_rho_lu':1.0004402376508774,'outlet_03_rho_lu':.9999543772756865}
 assert all(math.isfinite(c[k]) and abs(c[k]-v)<=8*math.ulp(v) for k,v in expected.items())
 assert c['tau']==1. and c['tau']>.5 and c['dt_s']>0 and c['command_inlet_mach_estimate']<.01
 assert all(.99<=c[k]<=1.01 for k in expected if k!='dt_s')
 old=(N/'inputs/OLD_SOLVER_PARAMETERS.txt').read_text().splitlines();new=(N/'contracts/solver_parameters.txt').read_text().splitlines();assert len(old)==len(new) and old[:2]==new[:2]
 a=old[2].split();b=new[2].split();assert all(a[i]==b[i] for i in range(len(a)) if i not in [7,8,9,len(a)-1]);assert float(b[7])==c['dt_s'] and float(b[9])==c['rho_kg_m3']
 pos=3
 for port in range(4):
  a=old[pos].split();b=new[pos].split();count=int(a[5]);assert all(a[i]==b[i] for i in range(len(a)) if i!=4);assert old[pos+1:pos+1+count]==new[pos+1:pos+1+count]
  assert float(b[4])==(1. if port==0 else c[f'outlet_0{port}_rho_lu']);pos+=count+1
 assert pos==len(old)
 for p in (N/'source').glob('*.hpp'):assert sha(p)==sha(B/'source'/p.name),'GPU core/header changed'
 save('verification/NUMERICS_STATIC_CHECK.json',dict(status='PASS',reference_reproduction='PASS',reference_rule='8ULP',tau=1.,pressure_unit_pa=c['pressure_unit_pa'],estimated_inlet_Mach=c['command_inlet_mach_estimate'],generated_dt_s=c['dt_s'],generated_outlet_rho_lu={k:c[k] for k in expected if k!='dt_s'},old_geometry_port_tokens_identical=True,source_generator_single_authority=True,scope='Static check only; runtime pending'))
 for p in (N/'contracts').iterdir():
  if p.is_file():p.chmod(p.stat().st_mode & ~0o222)
 tbb=B/'tbb/usr/lib/x86_64-linux-gnu';T=Path('/workspace/hemocell_restore/toolchains/nvhpc_26_5/Linux_x86_64/26.5')
 command('configure',['cmake','-S',str(N/'scripts/recipe'),'-B',str(N/'build'),'-DCMAKE_BUILD_TYPE=Release','-DCMAKE_C_COMPILER='+str(T/'compilers/bin/nvc'),'-DCMAKE_CXX_COMPILER='+str(T/'compilers/bin/nvc++'),'-DTBB_DIR='+str(tbb/'cmake/TBB'),'-DCMAKE_EXE_LINKER_FLAGS=-L'+str(tbb)])
 command('compile',['cmake','--build',str(N/'build'),'--parallel','1','--verbose'])
 exe=N/'build/vascularPoC';p=subprocess.run(['ldd',str(exe)],capture_output=True,text=True);assert p.returncode==0 and 'not found' not in p.stdout
 save('BUILD_PROVENANCE.json',dict(status='PASS',compiler='NVHPC26.5',cuda='13.2',GPU_arch='cc89',Palabos_commit='4127697e90169bbef982295f1d1c933cf6e90caa',native_library_sha256=sha(B/'build/official/libpalabos.a'),GPU_core_headers_sha256={p.name:sha(p) for p in (N/'source').glob('*.hpp')},new_driver_sha256=sha(N/'source/vascularPoC.cpp'),binary_sha256=sha(exe),binary=str(exe),ldd=p.stdout,flags=(N/'build/CMakeFiles/vascularPoC.dir/flags.make').read_text(),link_command=(N/'build/CMakeFiles/vascularPoC.dir/link.txt').read_text(),old_stage4_revalidated=False,old_binary_executed=False))
 state=dict(state='PASS',phase='BUILD_READY',unix=time.time(),solver_runs=0)
except BaseException as e:
 traceback.print_exc();state=dict(state='FAIL',phase=phase,error=repr(e),traceback=traceback.format_exc(),unix=time.time(),solver_runs=0)
save('BUILD_STATE.json',state);sys.exit(0 if state['state']=='PASS' else 2)
