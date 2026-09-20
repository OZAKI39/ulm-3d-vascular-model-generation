from sv13j_support import *
def test_actual_cpu_mpi_cuda_self_tests():
    d=actual('petsc_self_test')
    assert all(d['observed'].values()) and d['check']['exit_code']==0
def test_three_actual_gpu_sparse_solves():
    d=actual('petsc_gpu_smoke');assert len(d['runs'])==3
    for r in d['runs']:
        assert r['exit_code']==0 and r['accepted'] and r['KSP']=='gmres' and r['converged_reason']>0
        assert r['iterations']>0 and r['true_relative_residual']<=1e-10 and r['solution_error_inf']<=1e-10
        assert sha256(ROOT/'logs/sv1_3j/remote'/Path(r['log']).name)==r['sha256']

