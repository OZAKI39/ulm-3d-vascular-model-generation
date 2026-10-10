"""Stage N gates: source identity, compatibility scope and matched numerical evidence."""
import hashlib,re
from .parallel_checks import (GateError,require,finite,ghost_run_gate,ghost_gate,science_gate,SCIENCE_LIMITS,
                    flow_gate as base_flow_gate,timing_gate as base_timing_gate)

def production_monitor_policy_gate(policy,frozen):
    for key in ('dt_s','save_interval_steps','steady_last_intervals','velocity_change_limit','flow_change_limit','mass_limit'):
        require(policy.get(key)==frozen.get(key),'STEADY_POLICY_CHANGED_'+key)
    return True

def lifecycle_gate(d):
    require(d.get('normal_exit') and d.get('observed_at_MPI_Finalize')==[0,1,1],'FINALIZE_ORDER_OR_COUNT')
    return True

def checkpoint_one_rank(path,step,dt):
    import struct,math,numpy as np
    from pathlib import Path
    p=Path(path);data=p.read_bytes()
    require(p.suffix=='.bin' and len(data)>=56,'INCOMPLETE_CHECKPOINT')
    h=struct.unpack_from('<8i3d',data)
    require(h[:3]==(1,1,1) and h[4:8]==(0,4,0,step),'CHECKPOINT_LAYOUT_OR_STEP')
    nodes=h[3];require(nodes>0 and len(data)==56+64*nodes,'INCOMPLETE_CHECKPOINT')
    require(math.isclose(h[8],step*dt,rel_tol=1e-12),'CHECKPOINT_TIME')
    values=np.frombuffer(data,dtype='<f8',offset=56)
    require(np.isfinite(values).all() and np.isfinite(h[8:]).all(),'NONFINITE_CHECKPOINT')
    return dict(status='PASS',path=str(p),step=step,time_s=h[8],nodes=nodes,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),complete_integration_history=True,fields=['Y_n velocity/pressure','A_n time derivative','step/time','equation initial norm'])

def source_gate(d):
    require(bool(re.fullmatch(r'[0-9a-f]{40}',d.get('commit',''))),'UNPINNED_PETSC_SOURCE')
    require(bool(re.fullmatch(r'[0-9a-f]{64}',d.get('source_archive_sha256',''))),'SOURCE_ARCHIVE_SHA')
    require(d.get('version')=='3.25.5' and d.get('source_patch') is False,'UNAPPROVED_PETSC_SOURCE')
    if d.get('candidate')=='official tag':require(d.get('tag')=='v3.25.5','WRONG_TAG')
    else:
        require(d.get('candidate')=='Candidate R' and d.get('branch')=='release','MAIN_FORBIDDEN')
        require(all(d.get(k) for k in ('tag_failure_evidence','specific_fix_commit','official_diff','root_cause')),'RELEASE_FALLBACK_UNJUSTIFIED')
        require(bool(re.fullmatch(r'[0-9a-f]{40}',d['specific_fix_commit'])),'FIX_NOT_PINNED')
    return True

def linkage_gate(d,prefix,version='3.25.5'):
    require(d.get('version')==version,'WRONG_PETSC_VERSION')
    require(d.get('resolved_PETSc')==prefix+'/lib/libpetsc.so.3.25.5' or d.get('resolved_PETSc')==d.get('expected_PETSc') and d.get('resolved_PETSc','').startswith(prefix+'/lib/'),'WRONG_PETSC_LINK')
    require(d.get('resolved_MPI')==d.get('expected_MPI') and bool(d.get('resolved_MPI')),'WRONG_MPI_LINK')
    return True

def compatibility_gate(patch,manifest):
    require(manifest.get('repair_iterations',0)<=3,'REPAIR_BUDGET_EXHAUSTED')
    names=re.findall(r'^\+\+\+ b/(.+)$',patch,re.M)
    require(all(n in {'Code/Source/solver/PetscLinearAlgebra.cpp','Code/Source/solver/petsc_impl.cpp','Code/Source/solver/petsc_impl.h'} for n in names),'SCIENCE_SOURCE_MODIFIED')
    additions='\n'.join(l[1:] for l in patch.splitlines() if l.startswith('+') and not l.startswith('+++'))
    require(not re.search(r'\(void\)\s*Petsc|catch\s*\(\.\.\.\)|ierr\s*=\s*0\s*;',additions),'PETSC_ERROR_SWALLOWED')
    require(not re.search(r'\b(?:KSPSetTolerances|KSPSolve|MatSetValues|petsc_set_values|petsc_solve)\s*\(',additions),'SCIENTIFIC_OPERATION_CHANGED')
    require(hashlib.sha256(patch.encode()).hexdigest()==manifest.get('patch_sha256'),'PATCH_NOT_REVIEWED')
    require(set(names)==set(manifest.get('modified_files',[])),'PATCH_FILE_SCOPE')
    return True

def flow_gate(d,gpu=True,steps=None):
    base_flow_gate(d,gpu=gpu,steps=steps)
    require(len(d.get('runtime_semantics',[]))==d['linear_solves'],'RUNTIME_SOLVER_VIEW_COVERAGE')
    for observed in d['runtime_semantics']:runtime_semantics_gate(observed)
    if steps is not None:
        require(d.get('requested_steps')==steps,'REQUESTED_TIMESTEPS_CHANGED')
        require({r['step'] for r in d['history']['linear_solves']}==set(range(1,steps+1)),'TIMESTEP_COVERAGE')
    return True

SOLVER_SEMANTICS=dict(KSP='gmres',restart=100,max_iterations=2000,rtol=1e-10,atol=1e-24,divergence=10000.,side='right',norm='UNPRECONDITIONED',diagonal_scale=True,PC='asm',overlap=2,sub_KSP='preonly',sub_PC='ilu',fill_level=2)

def parse_runtime_semantics(log):
    blocks=re.findall(r'^KSP Object:.*?(?=^KSP Object:|^\s*NS\s+\d+-|\Z)',log,re.M|re.S)
    rows=[]
    for b in blocks:
        def get(pattern,cast=str):
            m=re.search(pattern,b,re.S)
            return cast(m[1]) if m else None
        rows.append(dict(KSP=get(r'\n\s*type:\s*(\S+)'),restart=get(r'restart=(\d+)',int),max_iterations=get(r'maximum iterations=(\d+)',int),rtol=get(r'tolerances:\s+relative=([^,]+)',float),atol=get(r'tolerances:\s+relative=[^,]+, absolute=([^,]+)',float),divergence=get(r'tolerances:\s+relative=[^,]+, absolute=[^,]+, divergence=([0-9.eE+-]+)',float),side=get(r'\b(left|right) preconditioning'),norm=get(r'using (\w+) norm type for convergence'),diagonal_scale='diagonally scaled system' in b,PC=get(r'\nPC Object:.*?\n\s*type:\s*(\S+)'),overlap=get(r'amount of overlap\s*=\s*(\d+)',int),sub_KSP=get(r'KSP Object: \(sub_\).*?\n\s*type:\s*(\S+)'),sub_PC=get(r'PC Object: \(sub_\).*?\n\s*type:\s*(\S+)'),fill_level=get(r'(\d+) levels of fill',int),ordering=get(r'matrix ordering:\s*(\S+)'),zero_pivot=get(r'tolerance for zero pivot\s+(\S+)',float)))
    return rows

def runtime_semantics_gate(d):
    for k,v in SOLVER_SEMANTICS.items():require(d.get(k)==v,'RUNTIME_SOLVER_SEMANTICS_'+k)
    require(d.get('ordering')=='natural' and d.get('zero_pivot')==2.22045e-14,'RUNTIME_ILU_DEFAULTS_CHANGED')
    return True

def timing_gate(runs):
    median=base_timing_gate(runs)
    require(len({r.get('xml_sha256') for r in runs})==1 and bool(runs[0].get('xml_sha256')),'BENCHMARK_INPUT_MISMATCH')
    require(all(r.get('requested_steps')==20 for r in runs),'BENCHMARK_TIMESTEPS_CHANGED')
    for k in ('solver_sha256','PETSc_library_sha256','PETSC_OPTIONS','mpi_ranks','OMP_NUM_THREADS'):
        require(len({r.get(k) for r in runs})==1 and runs[0].get(k) is not None,'BENCHMARK_IDENTITY_MISMATCH_'+k)
    return median
