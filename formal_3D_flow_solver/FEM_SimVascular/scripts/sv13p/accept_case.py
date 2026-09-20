"""Read actual solver evidence and independently reload the final native files."""
import json,math,re,subprocess,sys,xml.etree.ElementTree as ET
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13p import *
from sv_validation.postprocess import SolutionMeasurements
from sv_validation.provenance import sha256
name=sys.argv[1];R=ROOT/'reports/sv1_3p';case=ROOT/'outputs/sv1_3p'/name
record=json.loads((R/'remote'/(name+'_execution.json')).read_text())
log=(ROOT/'logs/sv1_3p/remote'/Path(record['log']).name).read_text(errors='replace')
header=ROOT/'external/petsc325/petsc-3.25.5/include/petscksp.h'
reasons={k:int(v) for k,v in re.findall(r'KSP_((?:CONVERGED|DIVERGED)_\w+)\s*=\s*(-?\d+)',header.read_text())}
errors=[];semantics=parse_runtime_semantics(log);reload=None;checkpoint=None;measurement=None
try:
 for semantic in semantics:semantics_gate(semantic,record['PETSC_OPTIONS'])
except (ValueError,KeyError) as exc:errors.append('FROZEN_RUNTIME: '+str(exc))
try:solver_health_gate(record,log,case/'solver.xml',reasons)
except (ValueError,KeyError,IndexError) as exc:errors.append('SOLVER_HEALTH: '+str(exc))
rows=record['history']['linear_solves'];final=max((r['step'] for r in rows),default=0)
baseline=ROOT/'outputs/sv1_3o/REAL_VASCULAR_GPU_PERF/solver.xml'
try:scientific_xml_gate(baseline,case/'solver.xml')
except ValueError as exc:errors.append(str(exc))
try:
 steps=[int(Path(f['path']).stem.rsplit('_',1)[1]) for f in record['results']]
 output_policy_gate(steps,final,record['start_step'])
 if record['start_step']:
  policy=json.loads((ROOT/'configs/sv1_3p/policy.json').read_text());require(record['initial_checkpoint_sha256']==policy['checkpoint_sha256'],'CHECKPOINT_CHANGED')
 for f in record['results']+record['checkpoints']:
  p=ROOT/'outputs/sv1_3p'/Path(f['path']).relative_to('outputs');require(sha256(p)==f['sha256'],'NATIVE_FILE_HASH')
 vtu=case/f'1-procs/result_{final:03d}.vtu';cp=case/f'1-procs/stFile_{final:03d}.bin'
 reload=json.loads(subprocess.check_output([sys.executable,'-B',ROOT/'scripts/sv13p/reload_field.py',vtu],text=True))
 require(reload['velocity_finite'] and reload['pressure_finite'],'NONFINITE_FIELD')
 checkpoint=checkpoint_one_rank(cp,final,record['dt_s']);require(checkpoint['nodes']==reload['points'],'CHECKPOINT_FIELD_NODE_MISMATCH')
 require(math.isclose(reload['time_s'],checkpoint['time_s'],rel_tol=1e-12),'VTU_CHECKPOINT_TIME_MISMATCH')
 if True:
  ref=json.loads((ROOT/'reports/sv1_3/reference_freeze.json').read_text());policy=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
  measure=SolutionMeasurements(ROOT/'outputs/sv1/SV_MESH/mesh_arrays.npz',ref['accepted_solution']['Q_target_m3_s'],policy['Umean_m_s'])
  u,p=measure.read(vtu);measurement=measure.measure(u,p)
  require(measurement['wall_noslip_pass'],'WALL_CONDITION')
  require(all(math.isfinite(v) for v in [measurement['Q_in_m3_s'],*measurement['outlet_flows_m3_s'].values(),measurement['epsilon_Q'],measurement['epsilon_mass']]),'NONFINITE_FLOW_DIAGNOSTICS')
except (ValueError,KeyError,FileNotFoundError,subprocess.CalledProcessError) as exc:errors.append('FIELD_OUTPUT: '+str(exc))
restart=int(re.search(r'-ksp_gmres_restart\s+(\d+)',record['PETSC_OPTIONS'])[1])
stats=iteration_statistics(record['history'],restart) if rows else {}
memory=[float(s['raw'].split(',')[1]) for s in record['GPU_samples'] if s['exit_code']==0 and s['raw']]
d=dict(status='FAIL' if errors else 'PASS',name=name,candidate=record['candidate'],mode=record['mode'],wall_time_s=record['wall_time_s'],start_step=record['start_step'],stop_step=final,steps=len({r['step'] for r in rows}),linear_failures=record['history']['failed_linear_solves'],nonlinear_failures=0 if not any(x.startswith('SOLVER_HEALTH') for x in errors) else None,errors=errors,statistics=stats,reload=reload,checkpoint=checkpoint,measurement=measurement,VTU_count=len(record['results']),VTU_steps=[int(Path(f['path']).stem.rsplit('_',1)[1]) for f in record['results']],initial_checkpoint_sha256=record['initial_checkpoint_sha256'],config_sha256=record['config_sha256'],solver_sha256=record['solver_sha256'],PETSc_library_sha256=record['PETSc_library_sha256'],PETSC_OPTIONS=record['PETSC_OPTIONS'],runtime_semantics=semantics,Mat='seqaijcusparse' if 'type: seqaijcusparse' in log else None,Vec='seqcuda' if 'type: seqcuda' in log else None,factor_packages=sorted(set(re.findall(r'package used to perform factorization:\s*(\S+)',log))),sampled_device_memory_max_MiB=max(memory) if memory else None,scientific_equivalence='DEFERRED',KSP_reason_values={k:reasons[k] for k in set(r['petsc_reason']['reason'] for r in rows if r.get('petsc_reason'))},profile_collected_in_timed_run=record['profiling'],CPU_fallback_warning_lines=[l for l in log.splitlines() if re.search(r'fallback|fall back',l,re.I)],residency='NOT_MEASURED',normal_exit=record['exit_code']==0)
profile_path=case/'petsc_profile.txt'
d['profile']=profile_events(profile_path.read_text()) if profile_path.exists() else None
d['observed_KSP_reasons']=record['history']['petsc_reasons']
d['observed_reason_iterations']=sum(r['iterations'] or 0 for r in d['observed_KSP_reasons'])
d['unassigned_last_iteration']=record['history']['unassigned_petsc_monitor'][-1]['iteration'] if record['history']['unassigned_petsc_monitor'] else None
d['reuse']=None
if record['candidate']=='P1' and d['status']=='PASS':
 try:d['reuse']=reuse_gate(log,d['profile'],record['history'])
 except (ValueError,KeyError) as exc:d['errors'].append('REUSE: '+str(exc));d['status']='FAIL'
d['gain_classification']=gain_classification(json.loads((R/'reference_freeze.json').read_text())['window_baseline']['wall_time_s'],d) if record['mode']=='window' else 'SMOKE_ONLY'
(R/(name+'_acceptance.json')).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
print(json.dumps({k:d[k] for k in ['name','status','wall_time_s','steps','statistics','errors']},ensure_ascii=False),flush=True)
raise SystemExit(0 if d['status']=='PASS' else 1)
