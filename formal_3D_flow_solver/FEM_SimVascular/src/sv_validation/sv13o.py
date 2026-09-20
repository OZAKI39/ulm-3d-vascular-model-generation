"""Single-rank performance gates; no CPU/GPU scientific-equivalence claim."""
import math,re,statistics,xml.etree.ElementTree as ET
from .sv13m import GateError,require
from .sv11 import linear_gate,nonlinear_gate
from .sv13n import parse_runtime_semantics as _parse_stage_N_semantics,checkpoint_one_rank

def parse_runtime_semantics(log):
    # PETSc PCView uses singular "1 level of fill" for the Stage O ILU(1)
    # candidate. Keep the frozen Stage N parser unchanged.
    return _parse_stage_N_semantics(re.sub(r'\b(\d+) level of fill\b',r'\1 levels of fill',log))

def scientific_xml_gate(baseline,candidate):
    allowed={'Continue_previous_simulation','Number_of_time_steps','Increment_in_saving_VTK_files','Increment_in_saving_restart_files'}
    def normalize(p):
        root=ET.parse(p).getroot()
        def entry(e,parent=''):
            text=(e.text or '').strip()
            if parent=='GeneralSimulationParameters' and e.tag in allowed:text='<runtime control>'
            return e.tag,sorted(e.attrib.items()),text,[entry(c,e.tag) for c in e]
        return entry(root)
    require(normalize(baseline)==normalize(candidate),'SCIENTIFIC_XML_CHANGED')
    return True

def semantics_gate(row,options):
    def opt(k):
        m=re.search(r'(?:^|\s)'+re.escape(k)+r'\s+(\S+)',options)
        return m[1] if m else None
    expected=dict(KSP='gmres',side='right',norm='UNPRECONDITIONED',rtol=1e-10,atol=1e-24,max_iterations=2000,divergence=10000.,diagonal_scale=True,restart=int(opt('-ksp_gmres_restart')),PC=opt('-pc_type'),ordering='natural',zero_pivot=2.22045e-14)
    if expected['PC']=='asm':expected.update(overlap=2,sub_KSP='preonly',sub_PC='ilu',fill_level=2)
    else:
        require(expected['PC']=='ilu','UNAUTHORIZED_PC')
        expected['fill_level']=int(opt('-pc_factor_levels'))
        require(row['sub_PC'] is None and row['sub_KSP'] is None,'ASM_WRAPPER_STILL_ACTIVE')
    require(expected['restart'] in (100,200) and expected['fill_level'] in (1,2),'UNAUTHORIZED_TUNING')
    for k,v in expected.items():require(row.get(k)==v,'RUNTIME_SEMANTICS_'+k)
    return True

def solver_health_gate(record,log,xml,reason_values):
    require(record['MPI_ranks']==record['GPUs']==record['OMP_NUM_THREADS']==1,'SINGLE_GPU_SCOPE')
    require(not re.search(r'PETSC ERROR|MPI_ERR_TYPE|MPI_ABORT|SIGSEGV|Resetting restart flag',log),'HARD_ERROR')
    linear_gate(record)
    eq=ET.parse(xml).find('Add_equation')
    nonlinear_gate(record['history'],int(eq.findtext('Max_iterations')),float(eq.findtext('Tolerance')),int(eq.findtext('Min_iterations')))
    rows=record['history']['linear_solves'];steps=sorted({r['step'] for r in rows})
    require(steps==list(range(record['start_step']+1,steps[-1]+1)),'INCOMPLETE_TIMESTEP_COVERAGE')
    if record['mode'] in ('window','profile'):require(steps==list(range(61,71)),'FIXED_WINDOW_CHANGED')
    for r in rows:
        reason=r['petsc_reason']['reason'];require(reason_values.get(reason,0)>0,'NONPOSITIVE_KSP_REASON')
    semantics=parse_runtime_semantics(log)
    require(len(semantics)==len(rows),'KSP_VIEW_COVERAGE')
    for s in semantics:semantics_gate(s,record['PETSC_OPTIONS'])
    mats=re.findall(r'(?:Mat|Matrix) Object:[^\n]*\n\s*type:\s*(\S+)',log)
    vecs=re.findall(r'(?:Vec|Vector) Object:[^\n]*\n\s*type:\s*(\S+)',log)
    require('seqaijcusparse' in mats and 'seqcuda' in vecs,'CUDA_BACKEND_NOT_OBSERVED')
    return semantics

def output_policy_gate(steps,final,start=0,cadence=10):
    expected=list(range((start//cadence+1)*cadence,final+1,cadence))
    if final not in expected:expected.append(final)
    require(sorted(steps)==expected,'OUTPUT_CADENCE_OR_FINAL_VTU')
    return True

def selection_action(incumbent_s,candidate):
    require(math.isfinite(incumbent_s) and incumbent_s>0,'INVALID_INCUMBENT_TIME')
    if candidate['status']!='PASS':return dict(action='REJECT_UNHEALTHY',reduction=None)
    t=candidate['wall_time_s'];require(math.isfinite(t) and t>0,'INVALID_CANDIDATE_TIME')
    reduction=1-t/incumbent_s
    action='ACCEPT' if t<=.90*incumbent_s else 'CONFIRM_ONCE' if t<=.95*incumbent_s else 'KEEP_INCUMBENT'
    return dict(action=action,reduction=reduction)

def confirm_selection(incumbent_s,first,second):
    require(selection_action(incumbent_s,first)['action']=='CONFIRM_ONCE','UNAUTHORIZED_REPEAT')
    if second['status']!='PASS':return dict(accepted=False,reason='CONFIRMATION_UNHEALTHY')
    times=[first['wall_time_s'],second['wall_time_s']];mean=statistics.mean(times)
    return dict(accepted=all(t<incumbent_s for t in times) and 1-mean/incumbent_s>=.05,arithmetic_mean_s=mean,reduction=1-mean/incumbent_s,repetitions=2,not_a_median=True)

def iteration_statistics(history,restart):
    values=[r['linear_iterations'] for r in history['linear_solves']]
    return dict(KSP_solves=len(values),total_iterations=sum(values),mean_iterations=statistics.mean(values),max_iterations=max(values),GMRES_cycles=sum(math.ceil(i/restart) for i in values),GMRES_restarts=sum(max(math.ceil(i/restart)-1,0) for i in values),cycle_count_source='Derived from actual iterations and fixed restart length; no extra trace run')
