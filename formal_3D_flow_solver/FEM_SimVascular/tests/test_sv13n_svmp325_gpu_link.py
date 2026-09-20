from sv13n_support import *
def test_modern_gpu_linkage():
 d=accepted('svmp_gpu_build');l=accepted('svmp_gpu_link');linkage_gate(l,d['PETSc_prefix'])
 assert 'sv1_3n/external/petsc325/install_gpu' in l['resolved_PETSc']
