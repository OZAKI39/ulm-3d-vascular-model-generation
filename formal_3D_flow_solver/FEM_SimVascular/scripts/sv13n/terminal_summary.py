"""Requested terminal record; no scientific equivalence or performance claim."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];R=ROOT/'reports/sv1_3n'
def read(n):return json.loads((R/(n+'.json')).read_text())
s=read('petsc325_source');b=read('petsc_gpu13_build');a=read('compatibility_adapter');smoke=read('svmp_gpu_smoke');c=read('gpu_steady_candidate');run=read('REAL_VASCULAR_GPU_acceptance');state=read('stage_result');tests=read('test_summary')
assert state['completed']
text=f"""Stage SV1.3N GPU development completed.

Modern PETSc:
    version = {s['version']}
    tag = {s['tag']}
    commit = {s['commit']}
    source SHA = {s['source_archive_sha256']}

GPU stack:
    CUDA = {b['CUDA_version']}
    MPI = Open MPI 4.1.6, full Fortran datatype support
    compiler = GCC/G++ 13.3; accepted by actual CUDA/PETSc configure
    GPU = RTX4090

compatibility repairs:
    repair_01 = PETSc object cleanup and exactly-once finalize before MPI; PASS
    repair_02 = CUDAFLAGS configure option / case-insensitive success-log parser; PASS
    repair_03 = NOT NEEDED

PETSc GPU:
    build = {b['status']}
    Mat = {read('petsc_cuda_types')['mat_type']}
    Vec = {read('petsc_cuda_types')['vec_type']}
    ghost forward = PASS (adopted 1 rank)
    ghost reverse = PASS (adopted 1 rank; no remote ghosts)
    optional CUDA MPI2 = FAIL, reverse/local-form stale owned data; NOT SUPPORTED
    standalone smoke = {read('petsc_gpu_smoke')['status']}

svMultiPhysics GPU:
    build = {read('svmp_gpu_build')['status']}
    linkage = {read('svmp_gpu_link')['status']}
    official smoke = {smoke['status']}
    runtime Mat = {run['mat_type']}
    runtime Vec = {run['vec_type']}
    KSP = GMRES, right, restart100, unchanged tolerance
    PC = ASM overlap2 / ILU(2)

real vascular GPU:
    started = YES, t=0
    linear failures = {run['linear_failures']}
    nonlinear failures = {run['nonlinear_failures']}
    velocity finite = {c['reload']['velocity_finite']}
    pressure finite = {c['reload']['pressure_finite']}

steady:
    first joint steady step = {c['first_full_steady_step']}
    required consecutive intervals = {c['required_consecutive_intervals']}
    safe stop step = {c['stop_step']}
    physical time = {c['physical_time_s']:.12g}
    E_u = {c['last_interval']['E_u']:.12g}
    E_Q = {c['last_interval']['E_Q']:.12g}
    mass = {c['measurement']['epsilon_mass']:.12g}
    reload = {c['reload']['status']}

GPU candidate:
    VTU = {c['VTU']}
    VTU SHA = {c['VTU_sha256']}
    checkpoint = {c['checkpoint']['path']}
    config SHA = {c['config_sha256']}
    solver binary SHA = {c['solver_sha256']}

observed runtime:
    wall time = {c['wall_time_s']:.6f} s
    steps = {c['steps']}
    KSP solves = {c['linear_solves']}
    total KSP iterations = {c['total_KSP_iterations']}
    sampled device memory max = {read('observed_gpu_memory')['sampled_device_memory_max_MiB']:g} MiB (device-wide, 15 s samples; not exact process peak)
    note = OBSERVATIONAL ONLY

deferred by user decision:
    OLD vs NEW PETSc science = DEFERRED
    CPU vs GPU science equivalence = DEFERRED
    formal CPU/GPU benchmark = DEFERRED
    transfer profiling = DEFERRED
    GPU-native PC tuning = DEFERRED

production:
    current accepted production =
    CPU_EARLY_STOP_PRODUCTION

new result classification:
    GPU_STEADY_CANDIDATE
    NOT YET SCIENTIFICALLY VALIDATED

tests:
    Stage N = {tests['stage']}
    Full suite = {tests['full']}
    Historical = {tests['historical']}

report:
    reports/sv1_3n/REPORT.md

STAGE SV1.3N STATUS:
    {state['status']}

All expensive computation stopped. No further CFD, comparison or benchmark.
"""
(R/'TERMINAL_SUMMARY.txt').write_text(text);print(text)
