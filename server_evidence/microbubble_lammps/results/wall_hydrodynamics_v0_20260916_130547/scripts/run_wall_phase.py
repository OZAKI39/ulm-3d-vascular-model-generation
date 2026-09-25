from pathlib import Path
import os,sys,subprocess,time,json,hashlib
R=Path(__file__).resolve().parents[1];W=Path('/workspace/microbubble_lammps/work/wall_hydrodynamics_v0_20260916_130547')
env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',HWLOC_COMPONENTS='-gl',CUDA_VISIBLE_DEVICES='0')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
phase=sys.argv[1];state=R/('PHASE_'+phase.upper()+'.json');assert not state.exists()
cs=json.loads((R/'contracts/CASE_INDEX.json').read_text());cs=cs[:6] if phase=='contact' else cs[6:]
save(state,{'status':'RUNNING','pid':os.getpid(),'cases':[c['name'] for c in cs]})
try:
 for c in cs:
  d=R/'cases'/c['name'];receipt=d/'RUN_RECEIPT.json';assert not receipt.exists() and not (d/'RUN_STATE.json').exists(),'NO_AUTOMATIC_RERUN'
  exe=W/'build'/('wall_lmp_gpu' if c['kokkos'] else 'wall_lmp_cpu');cmd=['mpirun','--allow-run-as-root','--bind-to','none','-np',str(c['mpi_ranks']),str(exe),('case.cfg' if phase=='contact' else 'case.selected.cfg')]+(['kokkos'] if c['kokkos'] else [])
  rec={'status':'STARTING','command':cmd,'binary_sha256':sha(exe),'math_library_sha256':sha(W/'build/libwall_math.so'),'case_contract_sha256':sha(d/'CASE_CONTRACT.json'),'table_sha256':sha(R/'tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5')};assert rec['table_sha256']==c['table_sha256'];save(receipt,rec);print('START',c['name'],flush=True)
  with (d/'RUN.log').open('w') as log:
   t=time.monotonic();p=subprocess.Popen(cmd,cwd=d,stdout=log,stderr=subprocess.STDOUT,env=env);rec['pid']=p.pid;save(receipt,rec)
   try:rc=p.wait(timeout=1200)
   except subprocess.TimeoutExpired:p.terminate();p.wait();raise RuntimeError('BOUNDED_RUNTIME_TIMEOUT')
  rec.update(returncode=rc,wall_seconds=time.monotonic()-t,status='FINISHED' if rc==0 else 'FAILED');save(receipt,rec)
  assert rc==0,(c['name'],(d/'RUN.log').read_text()[-2000:])
  s=json.loads((d/'RUN_STATE.json').read_text());print('TERMINAL',c['name'],json.dumps(s),flush=True)
  assert s['reason'] in ['MAX_PHYSICAL_TIME','BLOCKED_LOCAL_PLANE_VALIDITY'],s
  if s['reason']=='MAX_PHYSICAL_TIME':assert abs(s['time_s']-s['target_time_s'])<1e-14
 save(state,{'status':'TERMINAL','phase':phase,'all_processes_exited':True})
except Exception as e:save(state,{'status':'FAILED','phase':phase,'error':repr(e)});raise
