from sv13m_support import *
from sv_validation.sv13m import *
def test_correct_link_and_clean_solver_source():
 d=accepted('svmp_gpu_link');assert d['resolved_PETSc']==d['expected_PETSc'] and 'ghostfix' in d['resolved_PETSc']
 assert d['resolved_MPI']==d['expected_MPI'] and d['source_unmodified']
