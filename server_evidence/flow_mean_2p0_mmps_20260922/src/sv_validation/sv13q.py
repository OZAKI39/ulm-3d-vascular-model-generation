"""Read-only Stage Q validation of native policy traces and scientific evidence."""
import re
from .sv13p import semantics_gate,profile_events
from .sv13o import scientific_xml_gate,output_policy_gate,iteration_statistics,parse_runtime_semantics
from .sv13m import require
from .sv11 import linear_gate,nonlinear_gate
from .sv13n import checkpoint_one_rank

def trace_rows(log):
 rows=[]
 for line in log.splitlines():
  if not line.startswith(('SV13Q_BEGIN ','SV13Q_END ','SV13Q_RECOVERY ')):continue
  d=dict(re.findall(r'(\w+)=([^ ]+)',line));d['kind']=line.split()[0]
  for k in ('logical','attempt','step','equation','reuse','age','ref','builds','attempts','reason','iterations','healthy','pending','old_reason','old_iterations','original_rhs_restored','initial_guess_zero'):
   if k in d:d[k]=int(d[k])
  if 'ratio' in d:d['ratio']=float(d['ratio'])
  rows.append(d)
 return rows

def attempts(log):
 trace=trace_rows(log);ends={(r['logical'],r['attempt']):r for r in trace if r['kind']=='SV13Q_END'}
 return [dict(r,outcome=ends.get((r['logical'],r['attempt']))) for r in trace if r['kind']=='SV13Q_BEGIN']

def validate_reuse(log,profile,history,candidate):
 rows=attempts(log);require(bool(rows),'NO_POLICY_TRACE')
 last=-1;ref=0;pending=False;builds=0;recoveries=0;seen=set();previous=None
 for n,r in enumerate(rows,1):
  e=r['outcome'];require(e is not None,'INCOMPLETE_ATTEMPT')
  require(r['policy']==candidate,'POLICY_CHANGED')
  expected_age=0 if last<0 else r['step']-last
  require(r['age']==expected_age,'BAD_REUSE_AGE')
  require(r['ref']==ref,'BAD_REFERENCE')
  if r['attempt']==1:
   require(candidate=='RA' and previous and previous['reuse']==1 and previous['logical']==r['logical'] and not previous['outcome']['healthy'],'ILLEGAL_RECOVERY')
   expected_reason='STALE_FAILURE';recoveries+=bool(e['healthy'])
  else:
   require(r['attempt']==0,'MULTIPLE_RETRIES')
   require(r['logical'] not in seen,'DUPLICATE_LOGICAL_SOLVE');seen.add(r['logical'])
   limit=5 if candidate=='RA' else int(candidate[1:])
   expected_reason='FIRST_SOLVE' if last<0 else 'ITERATION_GROWTH' if candidate=='RA' and pending else ('MAX_AGE' if candidate=='RA' else 'FIXED_INTERVAL') if expected_age>=limit else 'REUSE'
  reuse=expected_reason=='REUSE'
  require(r['rebuild_reason']==expected_reason and r['reuse']==int(reuse),'POLICY_DECISION_MISMATCH')
  if not reuse:last=r['step'];builds+=1;pending=False
  require(r['builds']==builds and r['attempts']==n,'BAD_COUNTERS')
  if e['healthy'] and not reuse:ref=e['iterations']
  if e['healthy'] and candidate=='RA' and reuse and e['iterations']>1.5*ref:pending=True
  require(e['ref']==ref and e['pending']==int(pending),'BAD_ADAPTIVE_UPDATE')
  require(abs(e['ratio']-(e['iterations']/ref if ref>0 else 0.))<1e-12,'BAD_TRIGGER_RATIO')
  previous=r
 require(len(seen)==len(history['linear_solves']),'LOGICAL_SOLVE_COVERAGE')
 require(len({r['ksp'] for r in rows})==len({r['pc'] for r in rows})==1,'KSP_PC_RECREATED')
 if profile:
  require(profile['events']['MatLUFactorNum']['count']==builds,'ACTUAL_FACTORIZATION_COUNT')
  require(profile['events']['KSPSolve']['count']==len(rows),'ACTUAL_KSP_ATTEMPT_COUNT')
 recovery=[r for r in trace_rows(log) if r['kind']=='SV13Q_RECOVERY']
 require(len(recovery)==sum(r['attempt']==1 for r in rows),'RECOVERY_TRACE_MISSING')
 for r in recovery:require(r['original_rhs_restored']==r['initial_guess_zero']==1,'RECOVERY_CHANGED_RHS_OR_GUESS')
 return dict(status='PASS',trace=rows,ILU_rebuild_count=builds,ILU_reuse_count=sum(r['reuse'] for r in rows),KSP_attempts=len(rows),logical_solves=len(seen),recovery_count=recoveries,failed_attempts=[r for r in rows if not r['outcome']['healthy']],recovery_trace=recovery)

def solver_health_gate(record,log,xml,reason_values):
 import xml.etree.ElementTree as ET
 require(record['MPI_ranks']==record['GPUs']==record['OMP_NUM_THREADS']==1,'SINGLE_GPU_SCOPE')
 require(not re.search(r'PETSC ERROR|MPI_ERR_TYPE|MPI_ABORT|SIGSEGV|Resetting restart flag|out of memory',log,re.I),'HARD_ERROR')
 linear_gate(record);eq=ET.parse(xml).find('Add_equation')
 nonlinear_gate(record['history'],int(eq.findtext('Max_iterations')),float(eq.findtext('Tolerance')),int(eq.findtext('Min_iterations')))
 rows=record['history']['linear_solves'];steps=sorted({r['step'] for r in rows});require(bool(steps),'NO_STEPS')
 require(steps==list(range(record['start_step']+1,steps[-1]+1)),'TIMESTEP_COVERAGE')
 if record['mode']!='full':require(steps==list(range(record['start_step']+1,record['end_step']+1)),'WINDOW_NOT_COMPLETED')
 for r in rows:require(reason_values.get(r['petsc_reason']['reason'],0)>0,'KSP_REASON')
 a=attempts(log);semantics=parse_runtime_semantics(log)
 require(len(semantics)==len(a),'KSP_VIEW_COVERAGE')
 for s in semantics:semantics_gate(s,record['PETSC_OPTIONS'])
 for r in a:
  e=r['outcome'];require(e is not None,'INCOMPLETE_ATTEMPT')
  if not e['healthy']:
   require(record['candidate']=='RA' and r['reuse'] and r['attempt']==0,'UNRECOVERABLE_FAILURE')
   retry=[x for x in a if x['logical']==r['logical'] and x['attempt']==1]
   require(len(retry)==1 and retry[0]['outcome']['healthy'],'FRESH_RETRY_FAILED')
 if record['mode']=='smoke':require(sum(not r['outcome']['healthy'] for r in a)<=1,'SMOKE_MULTIPLE_STALE_FAILURES')
 require('type: seqaijcusparse' in log and 'type: seqcuda' in log,'GPU_BACKEND')
 return semantics

def gain_classification(baseline_s,candidate):
 if candidate['status']!='PASS':return 'REJECT_UNHEALTHY'
 ratio=candidate['wall_time_s']/baseline_s
 return 'STRONG_GAIN' if ratio<=.8 else 'USEFUL_GAIN' if ratio<=.9 else 'SMALL_GAIN' if ratio<=.95 else 'NO_MEANINGFUL_GAIN'
