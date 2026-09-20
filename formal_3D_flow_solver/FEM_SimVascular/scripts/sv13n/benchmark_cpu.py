"""Serial same-host, t=0 to 20 timing: two samples plus third on >10% spread."""
raise SystemExit('CANCELLED by USER_STRATEGY_UPDATE.txt: no repeated CFD benchmarks')
import json,statistics,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13n import timing_gate
R=ROOT/'reports/sv1_3n';S=ROOT/'scripts/sv13n'
assert json.loads((R/'old_new_cpu_science.json').read_text())['status']=='PASS'
groups={};suffix=sys.argv[1] if len(sys.argv)>1 else 'INITIAL'
for version,stack in [('OLD','old'),('NEW','cpu')]:
 for ranks in (1,4):
  key=f'{version}_CPU{ranks}';runs=[]
  if ranks==1 and suffix=='INITIAL':
   name=f'{version}_PETSC_CPU_PROOF_20'
   d=json.loads((R/(name+'_acceptance.json')).read_text());assert d['accepted'] and not d['profiling']
   d['role']='warm-up (complete non-profiled science proof; identical I/O)';runs.append(d)
  for i in range(len(runs)+1,4):
   if i==3 and abs(runs[1]['wall_time_s']-runs[0]['wall_time_s'])/runs[0]['wall_time_s']<=.10:break
   name=f'BENCH_{key}_{suffix}_{i}'
   subprocess.run([sys.executable,'-B',S/'run_case.py',name,str(ranks),'CPU',stack],check=True)
   d=json.loads((R/(name+'_acceptance.json')).read_text());d['role']='warm-up' if i==1 else 'measurement'
   runs.append(d)
  median=timing_gate(runs);times=[r['wall_time_s'] for r in runs]
  count=[sum(x['linear_iterations'] for x in r['history']['linear_solves']) for r in runs]
  solves=[r['linear_solves'] for r in runs]
  groups[key]=dict(median_s=median,raw_s=times,spread=(max(times)-min(times))/statistics.median(times),samples=len(runs),sample_roles=[r['role'] for r in runs],first_two_difference_relative=abs(times[1]-times[0])/times[0],names=[r['name'] for r in runs],iterations=count,linear_solve_counts=solves,wall_per_step_s=[t/20 for t in times],wall_per_KSP_iteration_s=[t/n for t,n in zip(times,count)],iterations_per_solve=[n/s for n,s in zip(count,solves)],KSP_time_per_solve='NOT_MEASURED: no heavy profiling in benchmark',xml_sha256=runs[0]['xml_sha256'],solver_sha256=runs[0]['solver_sha256'],PETSc_library_sha256=runs[0]['PETSc_library_sha256'],options=runs[0]['PETSC_OPTIONS'],ranks=ranks)
  (R/('cpu_benchmark_'+suffix+'.json')).write_text(json.dumps(dict(status='RUNNING',groups=groups),indent=2)+'\n')
speedups={str(n):groups[f'OLD_CPU{n}']['median_s']/groups[f'NEW_CPU{n}']['median_s'] for n in (1,4)}
classify=lambda s:'NEW_PETSC_CPU_SIMILAR' if abs(s-1)<=.05 else 'NEW_PETSC_CPU_FASTER' if s>1 else 'NEW_PETSC_CPU_SLOWER'
d=dict(status='PASS',groups=groups,speedups=speedups,classifications={k:classify(v) for k,v in speedups.items()},median_policy='All two/three samples, including warm-up',same_server=True,same_window=True)
for n in (1,4):
 a,b=groups[f'OLD_CPU{n}'],groups[f'NEW_CPU{n}'];assert a['options']==b['options'] and a['xml_sha256']==b['xml_sha256']
 d.setdefault('iteration_comparison',{})[str(n)]=dict(old=a['iterations'],new=b['iterations'],options_equal=True,comment='Printed histories and exact frozen runtime options retained; do not infer kernel speed from total alone.')
(R/('cpu_benchmark_'+suffix+'.json')).write_text(json.dumps(d,indent=2)+'\n');(R/'cpu_version_benchmark.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d,indent=2))
