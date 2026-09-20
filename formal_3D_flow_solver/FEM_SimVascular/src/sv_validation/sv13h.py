"""Acceptance gates for the isolated CUDA12 bring-up; no scientific solver edits."""
import math,statistics
from pathlib import Path
from .sv13g import GateError,require,cuda_types_gate,mpi_process_gate,linkage_gate,restart_gate,proof_gate,equivalence_gate

def cuda_environment_gate(record,expected_prefix):
    prefix=Path(expected_prefix)
    require(record.get('wrapper_enabled') is True,'CUDA12_WRAPPER_NOT_ACTIVE')
    require('release 12.6' in record['nvcc_version'],'WRONG_CUDA_TOOLKIT')
    require(Path(record['nvcc_path']).is_relative_to(prefix),'WRONG_CUDA_TOOLKIT')
    for key in ('libcudart_path','cusparse_path','curand_path','include_path'):
        require(Path(record[key]).is_relative_to(prefix),'CUDA_LIBRARY_PREFIX_MISMATCH: '+key)

def unchanged_driver_and_cuda13(before,after):
    require(before['driver_files']==after['driver_files'],'DRIVER_CHANGED')
    require(before['cuda_default_link']==after['cuda_default_link'],'SYSTEM_CUDA_LINK_CHANGED')
    require(before['cuda13_prefix']==after['cuda13_prefix'] and before['cuda13_files']==after['cuda13_files'],'CUDA13_CHANGED')
    require(before['ld_configuration']==after['ld_configuration'],'SYSTEM_LINKER_CONFIGURATION_CHANGED')

def petsc_cuda_build_gate(record):
    require(record['configure_exit']==record['make_exit']==0,'PETSC_CUDA12_BUILD_FAIL')
    require(record.get('cuda_enabled') is True,'PETSC_CUDA_NOT_ENABLED')
    require(record.get('library_sha256') and record.get('source_unmodified') is True,'PETSC_BUILD_PROVENANCE_MISSING')

def gpu_smoke_gate(record):
    require(record['exit_code']==0,'GPU_SMOKE_PROCESS_FAILED')
    cuda_types_gate(record['mat_type'],record['vec_type'])
    require(record['converged_reason']>0 and record['linear_failures']==0,'GPU_SMOKE_KSP_DIVERGED')
    require(math.isfinite(record['true_relative_residual']) and record['true_relative_residual']<=record['tolerance'],'GPU_SMOKE_RESIDUAL_FAILED')

def science_gate(metrics,policy):
    equivalence_gate(metrics,policy)
    v=metrics['max_velocity_relative']
    require(math.isfinite(v) and 0<=v<=policy['max_velocity_relative'],'GPU_NUMERICAL_DIFFERENCE: max_velocity')

def median_benchmark(runs):
    require(len(runs)>=2,'BENCHMARK_REPEATABILITY_MISSING')
    require(all(x.get('profiling') is False for x in runs),'PROFILING_NOT_A_BENCHMARK')
    require(all(x['steps']==20 and x['initial_state']=='t=0' and x['science_pass'] for x in runs),'BENCHMARK_WINDOW_OR_SCIENCE_INVALID')
    times=[x['wall_time_s'] for x in runs]
    require(all(math.isfinite(t) and t>0 for t in times),'BENCHMARK_TIME_INVALID')
    if abs(times[1]-times[0])/times[0]>.1:require(len(times)>=3,'BENCHMARK_THIRD_REPEAT_REQUIRED')
    return statistics.median(times)

def performance_classification(cpu1,cpu4,gpu,science_pass,residency):
    require(science_pass,'GPU_SCIENCE_FAIL')
    speedup=min(median_benchmark(cpu1),median_benchmark(cpu4))/median_benchmark(gpu)
    if residency=='GPU_TRANSFER_BOUND':classification='GPU_TRANSFER_BOUND'
    elif residency!='GPU_RESIDENCY_PASS':classification='BLOCKED'
    elif speedup>=1.25:classification='GPU_A_PROMISING'
    elif speedup>1:classification='GPU_WORKS_BUT_NOT_WORTH_COMPLEXITY'
    else:classification='GPU_SLOWER_THAN_CPU'
    return {'speedup':speedup,'classification':classification,'production':'CPU_EARLY_STOP_PRODUCTION'}
