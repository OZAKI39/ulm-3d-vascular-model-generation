from sv13m_support import *
from sv_validation.sv13m import *
def test_actual_patch_scope():
 d=accepted('patch_integrity');patch_scope_gate(d['modified_files'],d['expected_files'],d['solver_changes'])
@pytest.mark.parametrize('extra,solver',[(['src/ksp/ksp/interface/itfunc.c'],[]),([],['petsc_impl.cpp'])])
def test_unexpected_change_rejected(extra,solver):
 expected=['src/vec/vec/impls/mpi/commonmpvec.c']
 with pytest.raises(GateError):patch_scope_gate(expected+extra,expected,solver)
