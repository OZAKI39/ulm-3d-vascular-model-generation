#!/usr/bin/env python3
"""Record blocked GPU prerequisites and select the verified CPU configuration."""
import json,sys
from pathlib import Path
import yaml
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from sv_validation.sv13 import *
from sv_validation.provenance import write_json,sha256,now
check_reference();cpu=load('cpu_validation');build=load('cuda_petsc_build');ref=load('reference_freeze')
assert cpu['status']=='PASS' and build['status']=='BLOCKED'
assert load('solution_reload')['status']=='PASS'
assert load('cpu_performance')['measured_continuation_speedup']>1
reason='远端 MPI 最小启动测试超时，PETSc 3.19.6 CUDA 配置未完成'
decision={'status':'GPU_BLOCKED','reason_code':'CUDA_PETSC_BUILD','reason':reason,
    'evidence':['outputs/sv1_3/remote_build_attempt2/configure_detail.log','reports/sv1_3/remote_mpi_probe.json'],
    'environment':'RTX 4090 / CUDA toolkit 13.2 / OpenMPI 4.1.6',
    'cuda_build_attempts':2,'first_attempt':'default library discovery requires unavailable nvToolsExt',
    'second_attempt':'explicit supported CUDA library list passes library checks; MPI launcher test times out',
    'PETSc_upgraded':False,'solver_scientific_source_changed':False,
    'CUDA_13_PETSc_source_incompatibility_proven':False,
    'gpu_A_status':'NOT_RUN_PREREQUISITE_BLOCKED','gpu_A_scientific_pass':False,
    'gpu_A_PC_or_transfer_dominated':None,'gpu_B_attempt_count':0,'gpu_B_status':'NOT_RUN',
    'gpu_C_or_later_attempt_count':0,'gpu_full_run':'NOT_RUN',
    'gpu_speedup':None,'new_production':'CPU_EARLY_STOP_PRODUCTION'}
write_json(REPORT/'gpu_decision.json',decision)
unavailable={'status':'NOT_RUN','reason':reason,'reason_code':'CUDA_PETSC_BUILD','measured':False}
for name in ('cuda_smoke','cuda_types','gpu_ksp_view','gpu_transfer_audit','gpu_memory',
             'fixed_window_equivalence','benchmark_repeatability','gpu_full_run'):
    write_json(REPORT/(name+'.json'),dict(unavailable))
write_json(REPORT/'gpu_residency_audit.json',dict(unavailable,matrix_resident=None,vectors_resident=None,
    transfer_scales_with=None,excessive_host_device_transfer=None,
    interpretation='No CUDA PETSc KSP executed. Hardware availability and utilization are not residency evidence.'))
write_json(ROOT/'benchmarks/sv1_3/status.json',dict(unavailable,
    CPU_1R=None,CPU_4R=None,GPU_1R=None,window_steps=frozen_policy()['benchmark_window_steps'],
    same_checkpoint_comparison_performed=False,speedup=None,
    checkpoint_note='Native checkpoint headers bind MPI rank count and local partition. No 4-rank checkpoint was silently reused as a 1-rank restart or replaced with VTU.'))
accepted=load('cpu_accepted_solution');selected=select_production(cpu)
write_json(REPORT/'production_selection.json',{'status':'PASS','selected':selected,
    'accepted_field':accepted['path'],'accepted_field_sha256':accepted['sha256'],
    'selection_basis':'Only candidate with completed numerical, steady, mass, equivalence, restart, reload and measured speed gates',
    'candidate_measurement':'same-WSL continuation from legal native step10',
    'GPU_candidates_excluded_before_performance_ranking':True,'timestamp':now()})
config={'designation':'FAST_PRODUCTION_CONFIGURATION','mode':'CPU_EARLY_STOP',
    'MPI_ranks':ref['PETSc']['mpi_ranks'],'GPU_count':0,'OMP_NUM_THREADS':ref['PETSc']['OMP_NUM_THREADS'],
    'PETSc_options':ref['PETSc']['PETSC_OPTIONS'],'PETSc_version':ref['build']['PETSc_version'],
    'svMultiPhysics_commit':ref['build']['commit'],'binary_sha256':ref['build']['binary_sha256'],
    'scientific_xml':'configs/sv1_2/sv_flow.xml',
    'validated_initial_state':{'kind':'full native checkpoint','step':load('cpu_execution')['initial_step'],
        'source':'outputs/sv1_2/vascular_flow/4-procs/stFile_010.bin','sha256':load('cpu_initial_state')['sha256'],
        'required_MPI_ranks':ref['PETSc']['mpi_ranks'],'preserves':['Y_n','A_n','physical time','equation initial residual']},
    'steady_stop_policy':{k:frozen_policy()[k] for k in ('velocity_change_limit','flow_change_limit','steady_last_intervals','save_interval_steps','mass_limit','stop_policy')},
    'expected_termination_mechanism':'native STOP_SIM at a saved timestep',
    'terminal_validation':['finite velocity and pressure','zero solver failures','wall no slip','mass','five steady intervals','complete native checkpoint','reference equivalence','fresh-process reload'],
    'monitor_implementation':'src/sv_validation/sv13.py:SteadyStopMonitor',
    'validation_runner':'scripts/sv13/run_cpu.py (one protected evidence run; do not overwrite completed output)',
    'termination_step':'computed dynamically; 80 is the observed validation result, not a configured target',
    'validation_reference':{'stage':'SV1.2','designation':'VALIDATION_REFERENCE','step':ref['accepted_solution']['step'],
        'path':ref['accepted_solution']['path'],'sha256':ref['accepted_solution']['sha256']},
    'GPU_status':decision['status'],'performance':load('cpu_performance')}
(CONFIG/'production_performance.yaml').write_text(yaml.safe_dump(config,sort_keys=False,allow_unicode=True))
print(selected,flush=True)
