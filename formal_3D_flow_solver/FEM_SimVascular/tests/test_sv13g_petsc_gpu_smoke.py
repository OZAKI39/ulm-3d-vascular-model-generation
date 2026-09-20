from sv13g_support import *
def test_three_actual_petsc_gpu_smokes():
    d=actual('petsc_gpu_smoke');assert len(d['runs'])==3
    for r in d['runs']:assert r['exit_code']==0 and r['converged_reason']>0 and r['relative_residual']<=1e-10
