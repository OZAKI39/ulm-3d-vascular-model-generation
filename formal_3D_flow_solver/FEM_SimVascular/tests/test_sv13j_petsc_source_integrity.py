from sv13j_support import *
def test_original_sources_unchanged_before_after_and_final():
    before=actual('cuda123_source_integrity_before');after=actual('cuda123_source_integrity_after');final=actual('final_source_integrity')
    assert before['source_files']==after['source_files']==final['PETSc_original_files']==11035
    assert before['modifications']==after['modifications']==final['PETSc_changes']==[]
    assert not final['svMultiPhysics_changes'] and not final['svMultiPhysics_git_status']
    assert not final['PETSc_upgraded'] and not final['source_patches_applied']
def test_patched_source_rejected():
    d=load('cuda123_source_integrity_before');d['modifications']=['vecseqcupm.hpp']
    with pytest.raises(GateError,match='SOURCE_MODIFIED'):source_integrity_gate(d)

