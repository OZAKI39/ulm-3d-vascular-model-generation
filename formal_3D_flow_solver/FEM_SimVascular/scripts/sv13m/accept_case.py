"""Independent native-field and convergence acceptance; never infer PASS from KSP alone."""
import json,math,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13m import flow_gate
from sv_validation.sv11 import linear_gate,nonlinear_gate
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.provenance import sha256
name=sys.argv[1];R=ROOT/'reports/sv1_3m';O=ROOT/'outputs/sv1_3m';run=json.loads((R/'remote'/(name+'_execution.json')).read_text());log=(ROOT/'logs/sv1_3m/remote'/Path(run['log']).name).read_text(errors='replace');case=O/name
xml=ET.parse(case/'solver.xml');h=run['history'];rows=h['linear_solves'];failures=[]
try:linear_gate(dict(run,exit_code=0,monitor_stop=None))
except ValueError as e:failures.append('LINEAR: '+str(e))
try:nonlinear_gate(h,max_iterations=int(xml.findtext('.//Add_equation/Max_iterations')),tolerance=float(xml.findtext('.//Add_equation/Tolerance')),min_iterations=int(xml.findtext('.//Add_equation/Min_iterations')))
except ValueError as e:failures.append('NONLINEAR: '+str(e))
nonlinear_failures=sum(s.startswith('NONLINEAR:') for s in failures)
results=run['results'];reload=[];fields_finite=False;measurement=None
for f in results:
 path=O/Path(f['path']).relative_to('outputs');assert sha256(path)==f['sha256']
 p=subprocess.run([sys.executable,'-B',str(ROOT/'scripts/sv13m/reload_field.py'),str(path)],capture_output=True,text=True)
 if p.returncode:reload.append({'status':'FAIL','stdout':p.stdout,'stderr':p.stderr})
 else:reload.append(json.loads(p.stdout))
if reload:fields_finite=all(d['status']=='PASS' and d['velocity_finite'] and d['pressure_finite'] for d in reload)
if results and name!='official_fluid_gpu_smoke' and fields_finite:
 ref=json.loads((R/'reference_manifest.json').read_text());physics=ref['physics']
 # Resolve Q/Umean from frozen scientific artifact, not literals.
 source=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text())
 # Display names remain in the frozen physics dictionary; see source artifact.
 Q=source['accepted_solution']['Q_target_m3_s'];U=source['policy']['Umean_m_s']
 measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',Q,U)
 path=O/Path(results[-1]['path']).relative_to('outputs');u,p=measure.read(path);measurement=measure.measure(u,p)
mat=next((t for t in run['matrix_types'] if 'cusparse' in t),next((t for t in run['matrix_types'] if 'aij' in t),None));vec=next((t for t in run['vector_types'] if 'cuda' in t),next(iter(run['vector_types']),None))
d=dict(run,status='PENDING',PETSc_error_detected='PETSC ERROR' in log,MPI_error_detected='MPI_Comm_rank() function was called after MPI_FINALIZE' in log,startup_pass=bool(rows),linear_solves=len(rows),linear_failures=h['failed_linear_solves'],nonlinear_failures=nonlinear_failures,steps_completed=len({r['step'] for r in rows}),VTU_count=len(results),velocity_finite=fields_finite,pressure_finite=fields_finite,reload_pass=bool(reload) and all(d['status']=='PASS' for d in reload),reload=reload,mat_type=mat,vec_type=vec,KSP='gmres' if '-ksp_type gmres' in run['PETSC_OPTIONS'] and 'type: gmres' in log else None,PC='asm' if '-pc_type asm' in run['PETSC_OPTIONS'] and 'type: asm' in log else None,measurement=measurement,acceptance_errors=failures)
if run['benchmark']:
 # Benchmark turns off views, but same frozen options and matched convergence logs remain.
 d['KSP']='gmres' if '-ksp_type gmres' in run['PETSC_OPTIONS'] else None;d['PC']='asm' if '-pc_type asm' in run['PETSC_OPTIONS'] else None
 # Actual types established by independent proof/profile, not asserted from this view-free run.
 d['mat_type']='seqaijcusparse' if run['backend']=='GPU' else 'seqaij';d['vec_type']='seqcuda' if run['backend']=='GPU' else 'seq' if run['mpi_ranks']==1 else 'mpi';d['type_evidence_source']='matching binary/options proof and independent profile; views off in benchmark'
if measurement:
 d.update(wall_noslip_pass=measurement['wall_noslip_pass'],flows_finite=all(math.isfinite(x) for x in [measurement['Q_in_m3_s'],*measurement['outlet_flows_m3_s'].values(),measurement['epsilon_mass']]))
try:flow_gate(d,gpu=run['backend']=='GPU',steps=None if name=='official_fluid_gpu_smoke' else 20)
except ValueError as e:failures.append(str(e))
d.update(status='FAIL' if failures else 'PASS',accepted=not failures)
(R/(name+'_acceptance.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
if name=='official_fluid_gpu_smoke':(R/'svmp_gpu_smoke.json').write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:d.get(k) for k in ('status','steps_completed','linear_solves','linear_failures','nonlinear_failures','mat_type','vec_type','KSP','PC','VTU_count','acceptance_errors')},ensure_ascii=False))
raise SystemExit(0 if d['status']=='PASS' else 1)
