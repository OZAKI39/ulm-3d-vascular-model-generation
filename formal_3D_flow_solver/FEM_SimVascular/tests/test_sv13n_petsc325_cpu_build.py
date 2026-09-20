from sv13n_support import *
def test_cpu_build_and_rank_tests():
 d=accepted('petsc_cpu_build');assert d['source_unmodified'] and d['CUDA_version']=='disabled'
 assert d['selftest_observed']=={'CPU':True,'MPI2':True}
 assert all(s['exit_code']==0 and not s['timeout'] for s in d['steps'])
