from sv13n_support import *
def test_correct_modern_gpu_library():
 d=accepted('svmp_gpu_build');linkage_gate(accepted('svmp_gpu_link'),d['PETSc_prefix'])
