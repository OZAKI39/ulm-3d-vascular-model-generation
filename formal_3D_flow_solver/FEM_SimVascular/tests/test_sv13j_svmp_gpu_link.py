from sv13j_support import *
def test_actual_same_commit_solver_and_winner_link():
    d=actual('svmp_gpu_link');linkage_gate(d['resolved_PETSc'],d['expected_PETSc'])
    b=actual('svmp_gpu_build')
    assert b['source_unmodified'] and b['source_modifications']==[]
    assert b['commit']=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
    assert 'not found' not in b['linked_libraries']
def test_cpu_petsc_link_rejected():
    with pytest.raises(GateError):linkage_gate('/cpu/libpetsc.so','/gpu/libpetsc.so')

