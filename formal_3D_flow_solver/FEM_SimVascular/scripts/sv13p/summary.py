"""Close candidate decisions using preserved runs, never projected or failed timings."""
import json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13p import choose_winner,gain_classification
R=ROOT/'reports/sv1_3p';C=ROOT/'configs/sv1_3p'
def read(p):return json.loads(p.read_text())
def maybe(p):return read(p) if p.exists() else None
ref=read(R/'reference_freeze.json');allrows=[];rows={}
for key,label in [('P1','同一步内复用 ILU'),('P2','hypre GPU ILU0'),('P3','hypre BoomerAMG'),('P4','PETSc GAMG'),('P5','NVIDIA AmgX')]:
 smoke=maybe(R/(key+'_SMOKE_acceptance.json'));window=maybe(R/(key+'_WINDOW_acceptance.json'));config=maybe(C/'candidates'/(key+'.json'))
 if window:result=gain_classification(ref['window_baseline']['wall_time_s'],window);allrows.append(window)
 elif smoke:result='GAMG_NOT_SUITABLE_CURRENT_MONOLITHIC_LAYOUT' if key=='P4' and smoke['status']=='FAIL' else 'REJECT_UNHEALTHY' if smoke['status']=='FAIL' else 'SMOKE_PASS_WINDOW_MISSING'
 elif key=='P5':result=read(R/'amgx_feasibility.json')['result']
 elif key in ('P2','P3'):result='BUILD_FAILED' if read(R/'remote/petsc_hypre_build.json')['status']=='FAIL' else 'SETUP_NOT_COMPLETED'
 else:raise RuntimeError('Candidate evidence missing: '+key)
 d=window or smoke
 rows[key]=dict(label=label,result=result,config=config,smoke=smoke,window=window,scientific_parameters_changed=False,scientific_equivalence='DEFERRED',timing_scope=d['mode'] if d else None,PC_setup_count=d['profile']['events'].get('PCSetUp',{}).get('count') if d and d.get('profile') else None,PC_rebuild_count=d.get('reuse',{}).get('PC_rebuild_count') if d and d.get('reuse') else None)
rows['P1']['strategy_result']='P1_PROMISING' if rows['P1']['result'] in ('USEFUL_GAIN','STRONG_GPU_CANDIDATE') else 'P1_REJECT'
w=choose_winner(ref['window_baseline']['wall_time_s'],allrows)
if w['fastest']:
 k=w['fastest'];w.update(candidate=k,**{a:rows[k]['config'][a] for a in ('PETSC_OPTIONS','build_report')})
w.update(baseline='Stage O PERF_A, reused unchanged',threshold=.10,rank_basis='Actual complete healthy 10-step process wall only',single_run=True,profiling_limit=read(R/'profiling_adjustment.json')['comparison_limit'])
(R/'winner.json').write_text(json.dumps(w,indent=2)+'\n')
(R/'candidate_summary.json').write_text(json.dumps(dict(status='COMPLETE',candidates=rows,baseline=ref['window_baseline'],winner=w),indent=2,ensure_ascii=False)+'\n')
if not w['full_required']:(R/'full_steady_decision.json').write_text(json.dumps(dict(status='NOT_REQUIRED',reason='No healthy >=10% fixed-window gain',full_run_count=0),indent=2)+'\n')
else:(R/'full_steady_decision.json').write_text(json.dumps(dict(status='REQUIRED',candidate=w['candidate'],reason='Fastest healthy >=10% fixed-window gain',maximum_runs=1),indent=2)+'\n')
print(json.dumps(w,indent=2))
