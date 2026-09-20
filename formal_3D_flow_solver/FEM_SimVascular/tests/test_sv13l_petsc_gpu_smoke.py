from sv13l_support import *
def test_actual_gpu_smoke_three_times():
    d=actual('petsc_gpu_smoke');assert len(d['runs'])==3
    for r in d['runs']:
        assert r['accepted'] and r['exit_code']==0 and r['converged_reason']>0
        cuda_types_gate(r['mat_type'],r['vec_type'])
        assert r['KSP']=='gmres' and r['true_relative_residual']<=1e-10 and r['solution_error_inf']<=1e-10
