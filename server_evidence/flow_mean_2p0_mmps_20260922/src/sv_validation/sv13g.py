"""Fail-closed acceptance gates for native single-GPU validation.

These gates distinguish unexecuted stages from observed passes. GPU utilization,
successful configuration, and exit zero alone are never scientific acceptance.
"""
from pathlib import Path
import math

class GateError(ValueError):pass

def require(condition,reason):
    if not condition:raise GateError(reason)

def mpi_process_gate(record,ranks):
    require(not record['timeout'] and record['exit_code']==0,'MPI_RUNTIME_UNUSABLE')
    require(set(record['stdout'].splitlines())=={f'rank={i} size={ranks}' for i in range(ranks)},'MPI_RANK_OUTPUT_INVALID')

def mpi_stack_gate(compiler_family,launcher_family,compiler_prefix,launcher_prefix):
    require(compiler_family==launcher_family and Path(compiler_prefix)==Path(launcher_prefix),'MPI_STACK_MISMATCH')

def explicit_root_gate(family,uid,command,environment):
    if uid==0 and family=='Open MPI':
        require('--allow-run-as-root' in command or (environment.get('OMPI_ALLOW_RUN_AS_ROOT')=='1' and environment.get('OMPI_ALLOW_RUN_AS_ROOT_CONFIRM')=='1'),'OPENMPI_ROOT_OVERRIDE_NOT_EXPLICIT')

def cuda_types_gate(mat_type,vec_type):
    require(mat_type.lower() in ('aijcusparse','seqaijcusparse','mpiaijcusparse'),'PETSC_GPU_TYPES_NOT_ACTIVE: Mat')
    require(vec_type.lower() in ('cuda','seqcuda','mpicuda'),'PETSC_GPU_TYPES_NOT_ACTIVE: Vec')

def linkage_gate(resolved_petsc,expected_petsc):
    require(Path(resolved_petsc)==Path(expected_petsc),'SVMULTIPHYSICS_LINKED_WRONG_PETSC')

def restart_gate(initial_state,checkpoint_ranks,run_ranks):
    require(initial_state=='t=0' or (initial_state=='native' and checkpoint_ranks==run_ranks),'INCOMPATIBLE_NATIVE_RESTART')

def proof_gate(record):
    require(record.get('executed') is True and record.get('steps')==20,'GPU_PROOF_NOT_COMPLETE')
    require(record['linear_failures']==0 and record['nonlinear_failures']==0,'GPU_PROOF_CONVERGENCE_FAIL')
    require(record['velocity_finite'] and record['pressure_finite'],'GPU_FIELD_NONFINITE')
    require(all(math.isfinite(x) for x in [record['Qin'],*record['Qout'].values(),record['mass_error']]),'GPU_FIELD_NONFINITE')
    require(record['wall_noslip'] and record['reload_pass'],'GPU_PROOF_QC_FAIL')

def equivalence_gate(metrics,policy):
    limits={'velocity_relative_L2':policy['velocity_relative_L2'],'pressure_relative_L2':policy['pressure_relative_L2'],'Qin_relative':policy['Qin_relative'],'mass_error_difference':policy['mass_error_difference']}
    for key,limit in limits.items():require(math.isfinite(metrics[key]) and 0<=metrics[key]<=limit,'GPU_NUMERICAL_DIFFERENCE: '+key)
    require(bool(metrics['Qout_relative']),'GPU_NUMERICAL_DIFFERENCE: missing outlets')
    for value in metrics['Qout_relative'].values():require(math.isfinite(value) and 0<=value<=policy['Qout_relative'],'GPU_NUMERICAL_DIFFERENCE: Qout')

def residency_gate(evidence):
    require(evidence.get('measured') is True,'GPU_RESIDENCY_NOT_MEASURED')
    cuda_types_gate(evidence['mat_type'],evidence['vec_type'])
    require(evidence['matrix_resident'] and evidence['vectors_resident'],'GPU_RESIDENCY_FAIL')
    require(not evidence['large_transfer_per_iteration'],'GPU_WORKS_BUT_EXCESSIVE_TRANSFER')
    require(evidence['pc_apply_device']=='GPU','GPU_PC_NOT_RESIDENT')

def benchmark_time(times,threshold=.1):
    require(len(times)>=2 and all(math.isfinite(t) and t>0 for t in times),'BENCHMARK_REPEATABILITY_MISSING')
    if abs(times[1]-times[0])/times[0]>threshold:
        require(len(times)>=3,'BENCHMARK_THIRD_REPEAT_REQUIRED')
        return sorted(times[:3])[1]
    return times[1]

def gpu_performance(cpu1,cpu4,gpu,science_pass,residency_pass):
    require(science_pass,'GPU_NUMERICAL_DIFFERENCE')
    require(residency_pass,'GPU_TRANSFER_BOUND')
    speedup=min(benchmark_time(cpu1),benchmark_time(cpu4))/benchmark_time(gpu)
    return {'speedup':speedup,'classification':'GPU_A_PROMISING' if speedup>=1.25 else 'GPU_WORKS_BUT_NOT_WORTH_COMPLEXITY' if speedup>1 else 'GPU_SLOWER_THAN_CPU'}
