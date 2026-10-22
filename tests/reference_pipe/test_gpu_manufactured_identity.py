"""Guards against silent fallback and accidental scientific-option changes."""
from pathlib import Path

from human_vessel_flow.reference_pipe import gpu_manufactured as MODULE


def test_requested_gpu_does_not_certify_native_backend():
    assert not MODULE.backend_evidence('-mat_type aijcusparse -vec_type cuda')['verified']
    assert not MODULE.backend_evidence('  type: mpiaij\n  type: mpi\n')['verified']


def test_cpu_local_factorization_is_allowed_with_cuda_objects():
    text = '  type: mpiaijcusparse\n  type: mpicuda\n  type: seqaij\n  type: lu\n'
    assert MODULE.backend_evidence(text)['verified']


def test_gpu_switch_keeps_scientific_options_and_detects_tolerance_change():
    common = '-ksp_rtol 1e-10 -ksp_atol 1e-24 -pc_type asm -sv_pc_adaptive_rebuild true'
    cpu = common + ' -mat_type aij -vec_type standard'
    gpu = common + ' -mat_type aijcusparse -vec_type cuda -log_view_gpu_time'
    assert MODULE.science_options(cpu) == MODULE.science_options(gpu)
    assert MODULE.science_options(cpu) != MODULE.science_options(gpu.replace('1e-10', '1e-8'))


def test_explicit_failures_stop_but_slow_convergence_does_not():
    assert MODULE.fatal_reason('Linear solve converged due to CONVERGED_RTOL iterations 1200') is None
    assert MODULE.fatal_reason('Linear solve did not converge due to DIVERGED_ITS')
    assert MODULE.fatal_reason('KSP unpreconditioned resid norm nan')
    assert MODULE.fatal_reason('CUDA error 2 out of memory')
