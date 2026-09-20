from sv13l_support import *
def test_three_actual_compiler_smokes():
    d=actual('compiler_environment')
    assert set(d['compilers'])=={'gcc-12','g++-12','gfortran-12'}
    assert all(c['number'].split('.')[0]=='12' for c in d['compilers'].values())
    assert len(d['smokes'])==3 and all(s['run']['exit_code']==0 and 'result=42' in s['run']['stdout'] for s in d['smokes'])
    assert d['default_compilers_unchanged']
def test_actual_nvcc_host_lookup_correction():
    d=actual('compiler_path_correction')
    assert d['host_gcc'].split('.')[0]=='12' and d['compile_check']['exit_code']==0
    assert d['fresh_source_required'] and load('petsc_configure_attempt1')['exit_code']!=0
