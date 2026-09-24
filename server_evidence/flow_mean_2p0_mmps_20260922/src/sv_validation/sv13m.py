"""Stage M acceptance gates: numerical evidence overrides labels and exit status."""
import math, statistics
class GateError(ValueError):pass
def require(condition,reason):
    if not condition:raise GateError(reason)
def finite(x):return isinstance(x,(float,int)) and math.isfinite(x)
def ghost_run_gate(r):
    require(r['exit_code']==0 and not r.get('timeout') and not r.get('hard_error'),'GHOST_API_ERROR')
    n=r['ranks'];require(len(r['types'])==len(r['locals'])==len(r['checks'])==n,'GHOST_RANK_COVERAGE')
    require({int(t[0]) for t in r['types']}==set(range(n)),'GHOST_RANK_COVERAGE')
    expected_type=('seq' if n==1 else 'mpi')+('cuda' if r['variant'].startswith('CUDA') else '')
    require(all(t[1]==expected_type for t in r['types']),'GHOST_BACKEND_TYPE')
    for rank,t,size,expected,present in r['locals']:
        require(present=='1' and int(size)==int(expected)==(4 if n==1 else 6),'GHOST_LOCAL_SIZE')
    for rank,f,b,size,expected,errors in r['checks']:
        require(finite(float(f)) and float(f)<=1e-13,'GHOST_FORWARD_VALUES')
        require(finite(float(b)) and float(b)<=1e-13,'GHOST_REVERSE_VALUES')
        require(size==expected and int(errors)==0,'GHOST_LOCAL_OR_VALUES')
    values=r['values'];require(len(values)==n*(4 if n==1 else 6)*2,'GHOST_VALUE_COVERAGE')
    keys={(ph,int(rank),int(i)) for ph,rank,i,a,e in values}
    require(keys=={(ph,rank,i) for ph in ('forward','reverse') for rank in range(n) for i in range(4 if n==1 else 6)},'GHOST_VALUE_COVERAGE')
    # Independently reconstruct each expected value; do not trust a probe's printed expected value.
    for ph,rank,i,actual,expected in values:
        rank,i=int(rank),int(i);a=float(actual);e=float(expected)
        want=rank*4+i+1 if i<4 else (1-rank)*4+(i-4)+1
        if ph=='reverse':
            if i>=4:want=100+10*rank+i-4
            elif n==2 and i<2:want+=100+10*(1-rank)+i
        require(finite(a) and finite(e) and abs(a-want)<=1e-13 and abs(e-want)<=1e-13,'GHOST_ELEMENT_MISMATCH')
    return True
def ghost_gate(d):
    require(len(d['runs'])==12,'GHOST_REPEATS')
    for name in ('CPU_seq','CUDA_seq','CPU_MPI','CUDA_MPI'):
        rows=[r for r in d['runs'] if r['variant']==name]
        require({r['repetition'] for r in rows}=={1,2,3} and len(rows)==3,'GHOST_VARIANT_COVERAGE')
        for r in rows:ghost_run_gate(r)
    return True
def patch_scope_gate(modified,expected,solver_changes=()):
    require(set(modified)==set(expected),'UNEXPECTED_PETSC_FILE_CHANGE')
    require(not solver_changes,'SVMULTIPHYSICS_SOURCE_CHANGED')
    require(all(p.startswith(('src/vec/','include/petsc/private/')) for p in modified),'PATCH_OUTSIDE_VEC_LAYER')
def flow_gate(r,gpu=True,steps=None):
    require(not r.get('PETSc_error_detected') and not r.get('MPI_error_detected') and not r.get('monitor_stop'),'PETSC_HARD_ERROR_OVERRIDES_CONVERGENCE')
    require(r['exit_code']==0,'FLOW_EXIT_ERROR')
    require(r['linear_solves']>0 and r['linear_failures']==r['nonlinear_failures']==0,'FLOW_CONVERGENCE')
    require(r['VTU_count']>0 and r['velocity_finite'] and r['pressure_finite'] and r['reload_pass'],'FLOW_FIELD_OR_RELOAD')
    if gpu:require(r['mat_type']=='seqaijcusparse' and r['vec_type']=='seqcuda','GPU_BACKEND_INACTIVE')
    require(r['KSP']=='gmres' and r['PC']=='asm','SOLVER_POLICY_CHANGED')
    if steps is not None:
        require(r['steps_completed']==steps and r['initial_state']=='t=0','PROOF_STEPS')
        require(r['wall_noslip_pass'] and r['flows_finite'],'PROOF_FIELDS')
    return True
SCIENCE_LIMITS={'velocity_relative_volume_L2':1e-5,'pressure_relative_volume_L2':1e-5,'relative_Qin':1e-6,'max_relative_Qout':1e-5,'mass_error_difference':1e-6,'relative_max_velocity':1e-5}
def science_gate(d):
    require(d['pressure_shift_applied'] is False,'PRESSURE_SHIFT_FORBIDDEN')
    for k,limit in SCIENCE_LIMITS.items():require(finite(d['errors'][k]) and 0<=d['errors'][k]<=limit,'SCIENCE_'+k)
    return True
def timing_gate(runs):
    require(len(runs) in (2,3),'BENCHMARK_REPEATABILITY')
    for r in runs:
        require(not r.get('profiling') and not r.get('heavy_profiling'),'PROFILE_IS_NOT_BENCHMARK')
        require(r.get('accepted') and r.get('steps_completed')==20 and r.get('initial_state')=='t=0','BENCHMARK_INCOMPLETE')
        require(finite(r['wall_time_s']) and r['wall_time_s']>0,'BENCHMARK_TIME')
    times=[r['wall_time_s'] for r in runs]
    require(abs(times[1]-times[0])/times[0]<=.10 or len(times)==3,'THIRD_RUN_REQUIRED')
    return statistics.median(times)
