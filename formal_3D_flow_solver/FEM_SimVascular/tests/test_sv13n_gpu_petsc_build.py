from sv13n_support import *
def test_modern_cuda_install_and_checks():
 d=accepted('petsc_gpu13_build');assert d['version']=='3.25.5' and d['CUDA_version']=='13.2'
 assert all(d['selftest_observed'][k] for k in ('CPU','MPI2','CUDA'))
 assert d['source_unmodified'] and all(s['exit_code']==0 for s in d['steps'])
