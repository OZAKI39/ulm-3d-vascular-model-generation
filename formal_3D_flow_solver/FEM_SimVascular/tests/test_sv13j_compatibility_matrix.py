from sv13j_support import *
def test_actual_matrix_winner_is_cuda123():
    d=load('compatibility_matrix')
    assert d['winner']=='cuda123' and matrix_order_gate(d['candidates'])=='12.3.2'
    compatibility_gate(candidate())
    b=actual('petsc_cuda123_build')
    assert b['configure']['exit_code']==b['make']['exit_code']==0 and b['source_unmodified']
    assert len(b['library_sha256'])==64
def test_configure_pass_make_fail_is_not_compatible():
    d=candidate();d['make']='FAIL'
    with pytest.raises(GateError,match='NOT_COMPATIBLE'):compatibility_gate(d)

