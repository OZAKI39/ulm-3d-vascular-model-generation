from sv_validation.sv13 import gpu_gate
from sv13_support import gpu_evidence
def good():
    return dict(mat_type='seqaijcusparse',vec_type='seqcuda',peak_memory_fraction=.2,
        scientific_equivalence=True,mass_pass=True,linear_convergence=True,matrix_resident=True,
        vectors_resident=True,transfer_scales_with='nonlinear_solves',speedup=2)
def test_transfer_per_iteration_rejected():
    d=good();d['transfer_scales_with']='KSP_iterations'
    assert gpu_gate(d)=='EXCESSIVE_HOST_DEVICE_TRANSFER'
def test_actual_event_audit_has_evidence():
    d=gpu_evidence('gpu_transfer_audit')
    assert d['actual_petsc_version']=='3.19.6' and d['event_names'] and d['separate_profiling_run']
