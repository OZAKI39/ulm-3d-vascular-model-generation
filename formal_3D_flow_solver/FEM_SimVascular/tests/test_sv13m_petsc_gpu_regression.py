from sv13m_support import *
from sv_validation.sv13m import *
def test_three_same_gpu_sparse_solves():
 d=accepted('petsc_gpu_smoke');old=read('reference_manifest')['stack']['petsc_gpu_smoke']['runs'][0]
 assert len(d['runs'])==3
 for r in d['runs']:
  assert r['accepted'] and r['mat_type']=='seqaijcusparse' and r['vec_type']=='seqcuda'
  assert r['iterations']==old['iterations']
  assert abs(r['true_relative_residual']-old['true_relative_residual'])<1e-14
  assert abs(r['solution_error_inf']-old['solution_error_inf'])<1e-14
