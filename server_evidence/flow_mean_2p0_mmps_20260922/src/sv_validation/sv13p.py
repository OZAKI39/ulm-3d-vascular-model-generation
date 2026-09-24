"""Bounded single-GPU preconditioner evaluation, with science frozen."""
import math,re
from .sv13o import scientific_xml_gate,output_policy_gate,iteration_statistics,parse_runtime_semantics
from .sv13m import require,GateError
from .sv11 import linear_gate,nonlinear_gate
from .sv13n import checkpoint_one_rank

def semantics_gate(row,options):
    expected=dict(KSP='gmres',side='right',norm='UNPRECONDITIONED',rtol=1e-10,atol=1e-24,max_iterations=2000,divergence=10000.,diagonal_scale=True,restart=100)
    for k,v in expected.items():require(row.get(k)==v,'FROZEN_KSP_'+k)
    pc=re.search(r'(?:^|\s)-pc_type\s+(\S+)',options)[1]
    require(pc in ('asm','hypre','gamg','amgx') and row['PC']==pc,'PC_STRATEGY')
    if pc=='asm':
        for k,v in dict(overlap=2,sub_KSP='preonly',sub_PC='ilu',fill_level=2,ordering='natural',zero_pivot=2.22045e-14).items():require(row.get(k)==v,'BASELINE_PC_'+k)

def solver_health_gate(record,log,xml,reason_values):
    import xml.etree.ElementTree as ET
    require(record['MPI_ranks']==record['GPUs']==record['OMP_NUM_THREADS']==1,'SINGLE_GPU_SCOPE')
    require(not re.search(r'PETSC ERROR|MPI_ERR_TYPE|MPI_ABORT|SIGSEGV|Resetting restart flag|out of memory',log,re.I),'HARD_ERROR')
    linear_gate(record);eq=ET.parse(xml).find('Add_equation')
    nonlinear_gate(record['history'],int(eq.findtext('Max_iterations')),float(eq.findtext('Tolerance')),int(eq.findtext('Min_iterations')))
    rows=record['history']['linear_solves'];steps=sorted({r['step'] for r in rows})
    require(steps==list(range(record['start_step']+1,steps[-1]+1)),'TIMESTEP_COVERAGE')
    if record['mode']!='full':require(steps==list(range(61,record['end_step']+1)),'WINDOW_NOT_COMPLETED')
    for r in rows:require(reason_values.get(r['petsc_reason']['reason'],0)>0,'KSP_REASON')
    semantics=parse_runtime_semantics(log);require(len(semantics)==len(rows),'KSP_VIEW_COVERAGE')
    for s in semantics:semantics_gate(s,record['PETSC_OPTIONS'])
    require('type: seqaijcusparse' in log and 'type: seqcuda' in log,'GPU_BACKEND')
    return semantics

def profile_events(text):
    events={};by_stage={};stages={};stage=None
    for line in text.splitlines():
        s=re.match(r'^\s*\d+:\s+(.+?):\s+([0-9.eE+-]+)\s+([0-9.]+)%',line)
        if s:stages[s[1]]=dict(time_s=float(s[2]),percent_total=float(s[3]))
        s=re.match(r'^--- Event Stage \d+: (.+)$',line)
        if s:stage=s[1];by_stage.setdefault(stage,{});continue
        m=re.match(r'^([A-Za-z][A-Za-z0-9_: ]*?)\s{2,}(\d+)\s+([0-9.]+)\s+([0-9.eE+-]+|n/a)\s+([0-9.]+|n/a)\s+',line)
        if m and stage is not None:
            event=m[1].rstrip();row=dict(count=int(m[2]),time_s=None if m[4]=='n/a' else float(m[4]),raw=line.strip());by_stage[stage][event]=row
            total=events.setdefault(event,dict(count=0,time_s=0.,stage_rows=[]));total['count']+=row['count'];total['time_s']=total['time_s']+row['time_s'] if total['time_s'] is not None and row['time_s'] is not None else None;total['stage_rows'].append(dict(stage=stage,**row))
    return dict(events=events,events_by_stage=by_stage,stages=stages,nested_not_additive=True)

def reuse_gate(log,profile,history):
    pattern=r'SV13P_REUSE step=(\d+) equation=(\d+) reuse=(\d+) rebuild_requests=(\d+) solve_requests=(\d+) ksp=(\S+) pc=(\S+)'
    rows=[dict(zip(('step','equation','reuse','rebuild_requests','solve_requests','ksp','pc'),[int(v) if i<5 else v for i,v in enumerate(m)])) for m in re.findall(pattern,log)]
    require(len(rows)==len(history['linear_solves']),'REUSE_TRACE_COVERAGE')
    seen=set()
    for r in rows:
        key=(r['equation'],r['step']);require(r['reuse']==int(key in seen),'REUSE_OUTSIDE_TIMESTEP');seen.add(key)
    require(len({r['ksp'] for r in rows})==len({r['pc'] for r in rows})==1,'OBJECT_RECREATED')
    factors=profile['events']['MatLUFactorNum']['count']
    require(factors==len(seen)<len(rows),'ACTUAL_REUSE_NOT_OBSERVED')
    return dict(status='PASS',scope='REUSE_WITHIN_TIMESTEP',trace=rows,PC_rebuild_count=factors,KSP_solve_count=len(rows))

def gain_classification(baseline_s,candidate):
    if candidate['status']!='PASS':return 'REJECT_UNHEALTHY'
    ratio=candidate['wall_time_s']/baseline_s
    return 'STRONG_GPU_CANDIDATE' if ratio<=.8 else 'USEFUL_GAIN' if ratio<=.9 else 'SMALL_GAIN' if ratio<=.95 else 'NO_MEANINGFUL_GAIN'

def choose_winner(baseline_s,candidates):
    eligible=[d for d in candidates if d['status']=='PASS' and d['mode']=='window' and d['steps']==10]
    fastest=min(eligible,key=lambda d:d['wall_time_s']) if eligible else None
    return dict(fastest=fastest['candidate'] if fastest else None,improvement=1-fastest['wall_time_s']/baseline_s if fastest else None,full_required=bool(fastest and fastest['wall_time_s']<=.9*baseline_s))
