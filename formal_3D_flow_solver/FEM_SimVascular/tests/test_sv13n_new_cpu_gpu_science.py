from sv13n_support import *
def test_new_cpu_gpu_native_fields():
 d=accepted('new_cpu_gpu_science');science_gate(d)
 a,b=[read(n+'_acceptance') for n in d['cases']]
 assert a['solver_sha256']==b['solver_sha256'] and a['PETSc_library_sha256']==b['PETSc_library_sha256']
