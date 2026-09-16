"""One-shot gated execution. No physical/threshold edits and no automatic solver retry."""
from pathlib import Path
import sys,os,json,subprocess,hashlib,time,signal,csv,shutil,traceback
ROOT=Path(__file__).resolve().parents[1]
WORK=Path('/workspace/microbubble_lammps/work/passive_transport_v0_20260915_233516')
ENGINE=Path('/workspace/microbubble_lammps/work/lammps_particle_engine_20260915_214759')
sys.path.insert(0,str(ROOT/'src'))
from finalize_passive_transport_v0 import integrity,audit_case,audit_all,convergence,read,table
import numpy as np
ENV=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1');os.environ.update(ENV)
def save(p,r):Path(p).write_text(json.dumps(r,indent=2,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def cmd(args,timeout=120):return subprocess.check_output(args,text=True,stderr=subprocess.STDOUT,env=ENV,timeout=timeout)
def state(s,**kw):save(ROOT/'EXECUTION_STATE.json',{'status':s,'epoch':time.time(),**kw});print(s,kw,flush=True)
def baseline_audit(stage):
 roots={'COUPLING':Path('/workspace/microbubble_lammps/results/palabos_lammps_coupling_v0_20260915_224019'),'LAMMPS_RESULTS':Path('/workspace/microbubble_lammps/results/lammps_particle_engine_20260915_214759'),'PBS_BSA_RESULTS':Path('/workspace/hemocell_restore/results/new_medium_smoke_20260915_161937')};checks={}
 for label,p in roots.items():
  mismatches=[];count=0
  for line in (p/'SHA256SUMS').read_text().splitlines():
   h,n=line.split('  ',1);count+=1
   if sha(p/n)!=h:mismatches.append(n)
  checks[label]={'files':count,'mismatches':mismatches};assert not mismatches
 expected=read(roots['COUPLING']/'provenance/BASELINE_SOURCE_BEFORE.json');bad=[name for name,h in expected.items() if sha(ENGINE/name)!=h];assert not bad
 checks['LAMMPS_SOURCE_BINARY']={'files':len(expected),'mismatches':bad};save(ROOT/'provenance'/('BASELINE_'+stage+'.json'),{'status':'PASS','checks':checks})
def rss_tree(pid):
 rows=[list(map(int,l.split())) for l in cmd(['ps','-eo','pid=,ppid=,rss=']).splitlines()];ids={pid}
 while True:
  new=ids|{p for p,parent,r in rows if parent in ids}
  if new==ids:break
  ids=new
 return sum(r for p,parent,r in rows if p in ids)
def run_case(sp):
 p=ROOT/'cases'/sp['name'];assert not (p/'RUN_STARTED.json').exists();assert integrity(ROOT)
 args=['mpirun','--allow-run-as-root','--bind-to','none','-np',str(sp['mpi_ranks']),str(WORK/'build'/('passive_lmp_'+sp['engine']))]
 if sp['engine']=='gpu':args+=['-k','on','g','1','-sf','kk','-pk','kokkos','neigh','half','newton','off','gpu/aware','off']
 args+=['-in','in.lammps','-log','log.lammps'];save(p/'RUN_STARTED.json',{'command':args,'epoch':time.time(),'no_automatic_retry':True,'manifest_sha256':sha(ROOT/'FROZEN_INPUT_SHA256SUMS')})
 gpu=None;raw=None
 if sp['engine']=='gpu':
  raw=(p/'GPU_RESOURCE_RAW.csv').open('w');gpu=subprocess.Popen(['nvidia-smi','--query-gpu=timestamp,utilization.gpu,memory.used,power.draw','--format=csv,noheader,nounits','-lms','100'],stdout=raw,stderr=subprocess.STDOUT,env=ENV)
 peak=0;begin=time.perf_counter();timeout=False
 with (p/'stdout.log').open('w') as log,(p/'PROCESS_RESOURCE_TRACE.csv').open('w') as f:
  writer=csv.writer(f);writer.writerow(['elapsed_s','process_tree_rss_kib']);proc=subprocess.Popen(args,cwd=p,env=ENV,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  while proc.poll() is None:
   elapsed=time.perf_counter()-begin;rss=rss_tree(proc.pid);peak=max(peak,rss);writer.writerow([elapsed,rss]);f.flush()
   if elapsed>600:timeout=True;os.killpg(proc.pid,signal.SIGTERM);break
   time.sleep(.2)
  try:rc=proc.wait(timeout=10)
  except subprocess.TimeoutExpired:os.killpg(proc.pid,signal.SIGKILL);rc=proc.wait()
 elapsed=time.perf_counter()-begin
 if gpu:
  gpu.terminate();gpu.wait(timeout=10);raw.close()
  with (p/'GPU_RESOURCE_TRACE.csv').open('w') as f:
   w=csv.writer(f);w.writerow(['timestamp','utilization_percent','memory_mib','power_w'])
   for row in csv.reader((p/'GPU_RESOURCE_RAW.csv').read_text().splitlines()):
    if len(row)==4:w.writerow([v.strip() for v in row])
 save(p/'RUN_METRICS.json',{'returncode':rc,'wall_seconds':elapsed,'process_tree_peak_rss_kib':peak,'timeout':timeout,'performance_gate':'NONE','command':args})
 assert rc==0,sp['name']+' process failed'
 result=audit_case(ROOT,sp['name']);save(p/'INDEPENDENT_CASE_AUDIT.json',result);print(json.dumps({'name':sp['name'],'status':result['status'],'seconds':elapsed,'failed':[k for k,v in result['checks'].items() if not v],'metrics':result['metrics']}),flush=True)
 assert result['status']=='PASS',sp['name']+' independent audit failed'
def main():
 with (ROOT/'VALIDATION_STARTED.json').open('x') as f:json.dump({'epoch':time.time(),'no_automatic_retries':True},f)
 assert integrity(ROOT);state('PREFLIGHT');mem=int(next(l.split()[1] for l in Path('/proc/meminfo').read_text().splitlines() if l.startswith('MemAvailable:')))/2**20;disk=shutil.disk_usage(ROOT).free/2**30;gpu=cmd(['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu','--format=csv,noheader,nounits']);assert '4090' in gpu and mem>4 and disk>5;assert not cmd(['nvidia-smi','--query-compute-apps=pid,process_name','--format=csv,noheader']).strip();save(ROOT/'RESOURCE_PREFLIGHT.json',{'status':'PASS','ram_available_gib':mem,'disk_free_gib':disk,'gpu':gpu})
 baseline_audit('BEFORE');state('WALL_DISTANCE_UNIT_AUDIT')
 cmd([str(WORK/'build/wall_audit'),str(ROOT/'provenance/closed_geometry_m.stl'),str(ROOT/'validation/WALL_QUERY_REFERENCE.csv'),str(ROOT/'validation/WALL_QUERY_CPP.csv')]);a=table(ROOT/'validation/WALL_QUERY_REFERENCE.csv');b=table(ROOT/'validation/WALL_QUERY_CPP.csv');err=float(np.max(np.abs(a['wall_distance_m']-b['distance_m'])));assert err<=1e-12;save(ROOT/'validation/WALL_DISTANCE_AUDIT.json',{'status':'PASS','queries':len(a),'maximum_error_m':err,'reference':'independent VTK exact triangle distance'})
 cases=read(ROOT/'contracts/CASE_INDEX.json');lookup={s['name']:s for s in cases}
 for case in [0,1,2]:
  for suffix in ['c050','c025','c0125']:
   name=f'case{case}_{suffix}';state('RUNNING_CASE',name=name);run_case(lookup[name])
 selection=convergence(ROOT);save(ROOT/'contracts/PASSIVE_TRANSPORT_TIMESTEP_CONTRACT.json',selection);assert selection['status']=='PASS';label='c025' if selection['selected_c_adv']==.25 else 'c0125';state('PRODUCTION_TIMESTEP_FROZEN',selected=selection['selected_c_adv'])
 names=['case3_'+q+'_'+label for q in ['d10','d50','d90']]+['case4_mpi1_'+label,'case4_mpi4_'+label]
 if 'case5_'+label in lookup:names.append('case5_'+label)
 names.append('kokkos_'+label)
 for name in names:state('RUNNING_CASE',name=name);run_case(lookup[name])
 state('FINALIZING');result=audit_all(ROOT);save(ROOT/'REMOTE_INDEPENDENT_FINALIZER.json',result);assert result['status']=='PASS';baseline_audit('AFTER');assert integrity(ROOT)
 binary={str(p):sha(p) for p in (WORK/'build').glob('*') if p.is_file() and p.name in ['libfrozen_flow_coupling.so','wall_audit','passive_lmp_cpu','passive_lmp_gpu']};save(ROOT/'provenance/BUILD_PROVENANCE.json',{'status':'PASS','binaries':binary,'compiler':cmd(['g++','--version']),'MPI':cmd(['mpirun','--version']),'python':sys.version,'LAMMPS_source_commit':'9c5ab448c78a14fd534619622162ba418d6a1fb1','baseline_core_modified':False})
 for f in (WORK/'logs').glob('build_*.log'):shutil.copy2(f,ROOT/'provenance'/f.name)
 state('PASS',actual_runs=9+len(names),local_geometry_audit_pending=True)
if __name__=='__main__':
 try:main()
 except Exception as exc:state('FAIL',exception=repr(exc),traceback=traceback.format_exc());raise
