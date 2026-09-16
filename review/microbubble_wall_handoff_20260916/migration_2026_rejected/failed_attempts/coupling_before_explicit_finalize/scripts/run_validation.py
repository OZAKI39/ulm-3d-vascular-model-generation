"""One-shot remote execution. Frozen inputs, serial case gates, no retries."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,time,signal,csv,shutil,traceback,platform
ROOT=Path(__file__).resolve().parents[1]
WORK=Path('/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019')
ENGINE=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
sys.path.insert(0,str(ROOT/'src'))
from finalize_palabos_lammps_coupling import integrity,unit_audit,case_audit,audit_all
ENV=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1')
os.environ.update(ENV)
def save(p,r):Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def command(args,timeout=120):return subprocess.check_output(args,text=True,stderr=subprocess.STDOUT,env=ENV,timeout=timeout)
def state(status,**kw):save(ROOT/'EXECUTION_STATE.json',{'status':status,'utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),**kw});print(status,kw,flush=True)
def tree_rss(pid):
 rows=[list(map(int,l.split())) for l in command(['ps','-eo','pid=,ppid=,rss=']).splitlines()];ids={pid}
 while True:
  new=ids|{p for p,parent,r in rows if parent in ids}
  if new==ids:break
  ids=new
 return sum(r for p,parent,r in rows if p in ids)
def baseline_audit(stage):
 before=ROOT/'provenance/BASELINE_SOURCE_BEFORE.json'
 if stage=='before':
  source={str(p.relative_to(ENGINE)):sha(p) for p in sorted((ENGINE/'upstream/lammps').rglob('*')) if p.is_file() and '.git' not in p.parts}
  for backend in ['cpu','gpu']:
   for file in ['lmp','liblammps.a']:source[f'build_{backend}/{file}']=sha(ENGINE/f'build_{backend}'/file)
  save(before,source)
 else:
  expected=json.loads(before.read_text());source={n:sha(ENGINE/n) for n in expected};assert source==expected,'LAMMPS original source/binary modified'
 checks={}
 for label,base in [('LAMMPS_RESULTS',Path('/workspace/microbubble_lammps/results/lammps_particle_engine_20260915_214759')),('PBS_BSA_RESULTS',Path('/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937'))]:
  fails=[];count=0
  for line in (base/'SHA256SUMS').read_text().splitlines():
   digest,name=line.split('  ',1);count+=1
   if sha(base/name)!=digest:fails.append(name)
  checks[label]={'count':count,'mismatches':fails};assert not fails,(label,fails)
 save(ROOT/'provenance'/('BASELINE_INTEGRITY_'+stage.upper()+'.json'),{'status':'PASS','lammps_source_binary_files':len(source),'checks':checks})
def run_case(sp):
 p=ROOT/'cases'/sp['name'];started=p/'RUN_STARTED.json';assert not started.exists(),'Refuse second solver launch'
 exe=WORK/'build'/('coupled_lmp_'+sp['engine']);args=['mpirun','--allow-run-as-root','--bind-to','none','-np',str(sp['mpi_ranks']),str(exe)]
 if sp['engine']=='gpu':args+=['-k','on','g','1','-sf','kk','-pk','kokkos','neigh','half','newton','off','gpu/aware','off']
 args+=['-in','in.lammps','-log','log.lammps']
 if sp['engine']=='gpu':
  args=['nsys','profile','--sample=none','--cpuctxsw=none','--trace=cuda,nvtx,osrt','--force-overwrite=false','--output='+str(p/'kokkos_host_fix_trace')]+args
 save(started,{'command':args,'env':{k:ENV[k] for k in ['HWLOC_COMPONENTS','OMP_NUM_THREADS','CUDA_VISIBLE_DEVICES']},'frozen_manifest_sha256':sha(ROOT/'FROZEN_INPUT_SHA256SUMS'),'epoch':time.time(),'no_automatic_retry':True})
 gpu=None;raw=None
 if sp['engine']=='gpu':
  raw=(p/'GPU_RESOURCE_RAW.csv').open('w');gpu=subprocess.Popen(['nvidia-smi','--query-gpu=timestamp,utilization.gpu,memory.used,power.draw','--format=csv,noheader,nounits','-lms','100'],stdout=raw,stderr=subprocess.STDOUT,env=ENV)
 peak=0;begin=time.perf_counter();timeout=False
 with (p/'stdout.log').open('w') as log,(p/'PROCESS_RESOURCE_TRACE.csv').open('w') as f:
  writer=csv.writer(f);writer.writerow(['elapsed_s','process_tree_rss_kib']);proc=subprocess.Popen(args,cwd=p,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll() is None:
   elapsed=time.perf_counter()-begin;rss=tree_rss(proc.pid);peak=max(peak,rss);writer.writerow([elapsed,rss]);f.flush()
   if elapsed>300:timeout=True;os.killpg(proc.pid,signal.SIGTERM);break
   time.sleep(.2)
  try:rc=proc.wait(timeout=10)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
 wall=time.perf_counter()-begin
 if gpu:
  gpu.terminate();gpu.wait(timeout=10);raw.close()
  with (p/'GPU_RESOURCE_TRACE.csv').open('w') as f:
   w=csv.writer(f);w.writerow(['timestamp','utilization_percent','memory_mib','power_w'])
   for row in csv.reader((p/'GPU_RESOURCE_RAW.csv').read_text().splitlines()):
    if len(row)==4:w.writerow([v.strip() for v in row])
 save(p/'RUN_METRICS.json',{'returncode':rc,'wall_seconds':wall,'process_tree_peak_rss_kib':peak,'timeout':timeout,'timing_role':'correctness run with diagnostics; GPU includes nsys overhead','command':args})
 if rc!=0:raise RuntimeError(f'{sp["name"]} process exit {rc}')
 r=case_audit(ROOT,sp['name']);save(p/'INDEPENDENT_CASE_AUDIT.json',r);print(json.dumps({'case':sp['name'],'status':r['status'],'seconds':wall,'failed':[k for k,v in r['checks'].items() if not v]}),flush=True)
 if r['status']!='PASS':raise RuntimeError(f'{sp["name"]} independent audit FAIL')
def main():
 marker=ROOT/'VALIDATION_STARTED.json'
 with marker.open('x') as f:json.dump({'epoch':time.time(),'no_retries':True},f)
 assert integrity(ROOT),'Frozen input mismatch';state('PREFLIGHT')
 free=shutil.disk_usage(ROOT).free/2**30;mem=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))/2**20
 gpu=command(['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu','--format=csv,noheader,nounits']);assert '4090' in gpu and free>5 and mem>4
 compute=command(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader']);assert not compute.strip(),'Existing GPU compute process'
 save(ROOT/'RESOURCE_PREFLIGHT.json',{'status':'PASS','available_ram_gib':mem,'free_disk_gib':free,'gpu':gpu,'existing_gpu_compute_processes':compute,'MPI_QUIRK':'HWLOC_COMPONENTS=-gl excludes GL display probing only'})
 baseline_audit('before')
 binaries={str(p):sha(p) for p in (WORK/'build').glob('*') if p.is_file() and (p.name.startswith('coupled_lmp_') or p.name in ['flow_audit_cli','libfrozen_flow_coupling.so'])}
 for filename in ['coupling_build_0.log','coupling_build_1.log','coupling_final_build.log']:
  source=WORK/'logs'/filename
  if source.exists():shutil.copy2(source,ROOT/'provenance'/filename)
 save(ROOT/'provenance/COUPLING_BUILD_PROVENANCE.json',{'binaries':binaries,'compiler':command(['g++','--version']),'MPI':command(['mpirun','--version']),'python':sys.version,'platform':platform.platform(),'LAMMPS_root':str(ENGINE),'adapter_type':'external embedding executable registers custom Fix; original static library unchanged','LAMMPS_core_modified':False,'baseline_sha256':{k:sha(ENGINE/k) for k in ['build_cpu/lmp','build_cpu/liblammps.a','build_gpu/lmp','build_gpu/liblammps.a']}})
 state('UNIT_AUDIT');out=ROOT/'validation/cpp';out.mkdir(exist_ok=False)
 jobs=[('uniform','UNIFORM_FIELD.h5'),('linear','LINEAR_FIELD.h5'),('real_nodes','FROZEN_FLOW_FIELD_V0.h5'),('real_interior','FROZEN_FLOW_FIELD_V0.h5'),('invalid_domain','UNIFORM_FIELD.h5'),('invalid_missing','MISSING_CORNER_FIELD.h5'),('invalid_solid','FROZEN_FLOW_FIELD_V0.h5'),('force','UNIFORM_FIELD.h5')]
 for label,field in jobs:
  command([str(WORK/'build/flow_audit_cli'),'force' if label=='force' else 'query',str(ROOT/'fields'/field),str(ROOT/'validation/queries'/(label+'.csv')),str(out/(label+'.csv'))])
 r=unit_audit(ROOT);save(ROOT/'validation/UNIT_AUDIT.json',r);assert r['status']=='PASS',r
 contract=json.loads((ROOT/'contracts/PALABOS_LAMMPS_COUPLING_V0_CONTRACT.json').read_text())
 for sp in contract['cases']:
  state('RUNNING_CASE',case=sp['name']);assert integrity(ROOT);run_case(sp)
 state('FINAL_AUDIT');r=audit_all(ROOT);save(ROOT/'REMOTE_INDEPENDENT_FINALIZER.json',r);assert r['status']=='PASS'
 baseline_audit('after');assert integrity(ROOT);state('PASS',cases=len(contract['cases']))
if __name__=='__main__':
 try:main()
 except Exception as exc:
  state('FAIL',exception=repr(exc),traceback=traceback.format_exc());raise
