from sv13j_support import *
def test_clean_source_architecture_and_no_old_objects():
    d=actual('cuda123_source_integrity_before');source_integrity_gate(d)
    b=actual('petsc_cuda123_build')
    assert '/sv1_3j/external/compat_cuda/' in b['source']
    assert b['PETSC_ARCH']=='arch-sv13j-cuda123'
    isolation_gate([b])
def test_reused_objects_rejected():
    d=load('cuda123_source_integrity_before');d['previous_objects_reused']=True
    with pytest.raises(GateError,match='OBJECT_REUSE'):source_integrity_gate(d)
def test_shared_architecture_rejected():
    with pytest.raises(GateError,match='ISOLATION'):
        isolation_gate([dict(source='/A',prefix='/pA',PETSC_ARCH='same'),dict(source='/B',prefix='/pB',PETSC_ARCH='same')])

