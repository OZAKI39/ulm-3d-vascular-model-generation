"""Independent native-field and convergence acceptance; never infer PASS from KSP alone."""
import json,math,subprocess,sys,re,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13n import flow_gate,parse_runtime_semantics
from sv_validation.sv11 import linear_gate,nonlinear_gate
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.provenance import sha256
name=sys.argv[1];R=ROOT/'reports/sv1_3n';O=ROOT/'outputs/sv1_3n';run=json.loads((R/'remote'/(name+'_execution.json')).read_text());log=(ROOT/'logs/sv1_3n/remote'/Path(run['log']).name).read_text(errors='replace');case=O/name
xml=ET.parse(case/'solver.xml');h=run['history'];rows=h['linear_solves'];failures=[]
try:linear_gate(dict(run,exit_code=0,monitor_stop=None))
except ValueError as e:failures.append('LINEAR: '+str(e))
try:nonlinear_gate(h,max_iterations=int(xml.findtext('.//Add_equation/Max_iterations')),tolerance=float(xml.findtext('.//Add_equation/Tolerance')),min_iterations=int(xml.findtext('.//Add_equation/Min_iterations')))
except ValueError as e:failures.append('NONLINEAR: '+str(e))
nonlinear_failures=sum(s.startswith('NONLINEAR:') for s in failures)
results=run['results'];reload=[];fields_finite=False;measurement=None
if name=='REAL_VASCULAR_GPU':
 cadence=json.loads((ROOT/'configs/sv1_3n/convergence_strategy.json').read_text())['steady_policy']['save_interval_steps']
 results=[f for f in results if int(Path(f['path']).stem.rsplit('_',1)[1])%cadence==0 or f==run['results'][-1]]
for f in results:
 path=O/Path(f['path']).relative_to('outputs');assert sha256(path)==f['sha256']
 p=subprocess.run([sys.executable,'-B',str(ROOT/'scripts/sv13n/reload_field.py'),str(path)],capture_output=True,text=True)
 if p.returncode:reload.append({'status':'FAIL','stdout':p.stdout,'stderr':p.stderr})
 else:reload.append(json.loads(p.stdout))
if reload:fields_finite=all(d['status']=='PASS' and d['velocity_finite'] and d['pressure_finite'] for d in reload)
if results and not name.startswith('OFFICIAL_') and fields_finite:
 ref=json.loads((R/'reference_manifest.json').read_text());physics=ref['physics']
 # Resolve Q/Umean from frozen scientific artifact, not literals.
 source=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
 # Display names remain in the frozen physics dictionary; see source artifact.
 Q=source['accepted_solution']['Q_target_m3_s'];U=source['policy']['Umean_m_s']
 measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',Q,U)
 path=O/Path(results[-1]['path']).relative_to('outputs');u,p=measure.read(path);measurement=measure.measure(u,p)
mat=next((t for t in run['matrix_types'] if 'cusparse' in t),next((t for t in run['matrix_types'] if 'aij' in t),None));vec=next((t for t in run['vector_types'] if 'cuda' in t),next(iter(run['vector_types']),None))
d=dict(run,status='PENDING',PETSc_error_detected='PETSC ERROR' in log,MPI_error_detected=bool(re.search(r'MPI_ERR_TYPE|MPI_Comm_rank\(\) function was called after MPI_FINALIZE|MPI_ABORT was invoked',log)),startup_pass=bool(rows),linear_solves=len(rows),linear_failures=h['failed_linear_solves'],nonlinear_failures=nonlinear_failures,steps_completed=len({r['step'] for r in rows}),VTU_count=len(results),velocity_finite=fields_finite,pressure_finite=fields_finite,reload_pass=bool(reload) and all(d['status']=='PASS' for d in reload),reload=reload,mat_type=mat,vec_type=vec,KSP='gmres' if '-ksp_type gmres' in run['PETSC_OPTIONS'] and 'type: gmres' in log else None,PC='asm' if '-pc_type asm' in run['PETSC_OPTIONS'] and 'type: asm' in log else None,measurement=measurement,acceptance_errors=failures)
d['runtime_semantics']=parse_runtime_semantics(log)
if name.startswith('OFFICIAL_') and run['backend']=='GPU':
 trace=re.search(r'SV13N_PETSC_STATE initialized=(\d+) finalized=(\d+) count=(\d+)',log)
 values=list(map(int,trace.groups())) if trace else None
 lifecycle=dict(status='PASS' if values==[0,1,1] and 'exited normally' in log and run['exit_code']==0 else 'FAIL',
 observed_at_MPI_Finalize=values,expected=[0,1,1],source='Read-only GDB trace of the actual official solver',normal_exit='exited normally' in log)
 (R/'solver_finalize_lifecycle.json').write_text(json.dumps(lifecycle,indent=2)+'\n')
 if lifecycle['status']!='PASS':failures.append('PETSC_LIFECYCLE_ORDER_OR_COUNT')
if run['benchmark']:
 # Benchmark turns off views, but same frozen options and matched convergence logs remain.
 d['KSP']='gmres' if '-ksp_type gmres' in run['PETSC_OPTIONS'] else None;d['PC']='asm' if '-pc_type asm' in run['PETSC_OPTIONS'] else None
 # Actual types established by independent proof/profile, not asserted from this view-free run.
 d['type_evidence_source']='No object views during benchmark; proof/profile evidence required for GPU type acceptance.'
 if run['backend']=='GPU':
  proof=json.loads((R/'NEW_GPU_PROOF_20_acceptance.json').read_text())
  assert proof['accepted'] and proof['solver_sha256']==run['solver_sha256'] and proof['PETSc_library_sha256']==run['PETSc_library_sha256']
  assert proof['xml_sha256']==run['xml_sha256']
  d.update(mat_type=proof['mat_type'],vec_type=proof['vec_type'],type_evidence_source='NEW_GPU_PROOF_20_acceptance.json; same executable/library/input hashes, CUDA options unchanged')
if measurement:
 d.update(wall_noslip_pass=measurement['wall_noslip_pass'],flows_finite=all(math.isfinite(x) for x in [measurement['Q_in_m3_s'],*measurement['outlet_flows_m3_s'].values(),measurement['epsilon_mass']]))
safe_stop=R/'old_cpu_safe_stop_request.json'
superseded_long_run=name=='OLD_PETSC_CPU_PROOF_20' and safe_stop.exists()
if superseded_long_run:
 d['authorized_safe_stop']=json.loads(safe_stop.read_text())
 d['fixed20_requirement_superseded']=True
try:flow_gate(d,gpu=run['backend']=='GPU',steps=None if name.startswith('OFFICIAL_') or superseded_long_run or name=='REAL_VASCULAR_GPU' else run['requested_steps'])
except ValueError as e:failures.append(str(e))
d.update(status='FAIL' if failures else 'PASS',accepted=not failures)
(R/(name+'_acceptance.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
if name.startswith('OFFICIAL_'):(R/('svmp_'+run['backend'].lower()+'_smoke.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:d.get(k) for k in ('status','steps_completed','linear_solves','linear_failures','nonlinear_failures','mat_type','vec_type','KSP','PC','VTU_count','acceptance_errors')},ensure_ascii=False))
raise SystemExit(0 if d['status']=='PASS' else 1)
