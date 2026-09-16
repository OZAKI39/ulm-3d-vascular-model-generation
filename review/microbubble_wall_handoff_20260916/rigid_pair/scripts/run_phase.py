from pathlib import Path
import os,sys,subprocess,json,time,hashlib,shutil
R=Path(__file__).resolve().parents[1];W=Path('/workspace/microbubble_lammps/work/stable_rigid_sphere_near_field_v0_20260916_101617');PY='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0')
def run_logged(cmd,log,cwd=R):
 with open(log,'w') as f:
  t=time.time();p=subprocess.Popen(cmd,cwd=cwd,env=env,stdout=f,stderr=subprocess.STDOUT);code=p.wait()
 return {'command':cmd,'cwd':str(cwd),'pid':p.pid,'returncode':code,'wall_seconds':time.time()-t}
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def run_case(c):
 D=R/'cases'/c['name']
 if (D/'RUN_RECEIPT.json').exists():
  old=json.loads((D/'RUN_RECEIPT.json').read_text());assert old['status']=='FINISHED' and (D/'RUN_STATE.json').exists(),'Failed solver must not restart implicitly'
  print('REUSE_COMPLETED_SOLVER_VALIDATE_ONLY',c['name'],flush=True)
  r=run_logged([PY,'-B',str(R/'scripts/validate_cases.py'),c['name']],D/'INDEPENDENT_VALIDATION.log');assert r['returncode']==0,(c['name'],(D/'INDEPENDENT_VALIDATION.log').read_text()[-2000:]);return
 exe=W/'build'/('rigid_lmp_gpu' if c['gpu'] else 'rigid_lmp_cpu');cmd=['mpirun','--allow-run-as-root','--bind-to','none','-np',str(c['mpi_ranks']),str(exe),'case.cfg']+(['kokkos'] if c['gpu'] else [])
 receipt={'status':'STARTING','command':cmd,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),'library_sha256':hashlib.sha256((W/'build/librigid_math.so').read_bytes()).hexdigest(),'case_contract_sha256':hashlib.sha256((D/'CASE_CONTRACT.json').read_bytes()).hexdigest()};save(D/'RUN_RECEIPT.json',receipt)
 print('START',c['name'],flush=True);r=run_logged(cmd,D/'RUN.log',D);receipt.update(r);receipt['status']='FINISHED' if r['returncode']==0 else 'FAILED';save(D/'RUN_RECEIPT.json',receipt)
 if r['returncode']:print((D/'RUN.log').read_text()[-3000:],flush=True);raise RuntimeError('Solver failed: '+c['name'])
 r=run_logged([PY,'-B',str(R/'scripts/validate_cases.py'),c['name']],D/'INDEPENDENT_VALIDATION.log');assert r['returncode']==0,(c['name'],(D/'INDEPENDENT_VALIDATION.log').read_text()[-2000:]);print('PASS',c['name'],json.loads((D/'RUN_STATE.json').read_text()),flush=True)
phase=sys.argv[1];statusfile=R/('PHASE_'+phase.upper()+'_STATE.json');save(statusfile,{'status':'RUNNING','pid':os.getpid()})
try:
 cases=json.loads((R/'contracts/CASE_INDEX.json').read_text())
 if phase in ('initial','initial_resume'):
  if phase=='initial':
   r=run_logged(['cmake','--build',str(W/'build'),'-j4'],R/'logs/FINAL_INTEGRATION_BUILD.log');assert r['returncode']==0
   shutil.copy2(W/'build/librigid_math.so',R/'build_math/librigid_math.so')
   r=run_logged([PY,'-B',str(R/'scripts/audit_math.py')],R/'logs/REMOTE_MATH_AUDIT.log');assert r['returncode']==0,(R/'logs/REMOTE_MATH_AUDIT.log').read_text()[-2500:]
  else:assert json.loads((R/'validation/MATHEMATICAL_RELEASE_GATE.json').read_text())['status']=='PASS'
  for c in cases[:6]:run_case(c)
  r=run_logged([PY,'-B',str(R/'scripts/validate_cases.py'),'cgap'],R/'logs/CGAP.log');assert r['returncode']==0
  g=json.loads((R/'validation/CGAP_CONVERGENCE.json').read_text());assert g['status']=='PASS';save(R/'contracts/SELECTED_C_GAP.json',g)
 elif phase=='remaining':
  g=json.loads((R/'contracts/SELECTED_C_GAP.json').read_text());assert g['selected_C_gap']==.4,'Prepare selected-case configs before launch if selection differs'
  assert json.loads((R/'validation/MATHEMATICAL_RELEASE_GATE.json').read_text())['status']=='PASS'
  for c in cases[6:]:run_case(c)
 elif phase=='passive_pair':
  run_case(next(c for c in cases if c['name']=='J2_passive_separated_pair'))
 else:raise ValueError(phase)
 save(statusfile,{'status':'PASS','phase':phase});print('PHASE PASS',phase,flush=True)
except Exception as e:
 save(statusfile,{'status':'FAIL','phase':phase,'error':repr(e)});raise
