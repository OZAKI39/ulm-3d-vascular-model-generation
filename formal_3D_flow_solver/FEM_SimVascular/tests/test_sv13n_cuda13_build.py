from sv13n_support import *
def test_cuda13_actual_build():
 d=accepted('petsc_gpu13_build');assert d['CUDA_version']=='13.2' and d['source_unmodified'] and d['selftest_observed']['CUDA']
