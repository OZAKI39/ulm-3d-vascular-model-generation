from sv13n_support import *
def test_one_required_gpu_sparse_solve():
 d=accepted('petsc_gpu_smoke');assert len(d['runs'])==1
 for r in d['runs']:
  assert r['accepted'] and 'cusparse' in r['mat_type'] and 'cuda' in r['vec_type']
  assert r['true_relative_residual']<=1e-10 and r['solution_error_inf']<=1e-10
