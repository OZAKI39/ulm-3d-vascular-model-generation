from sv13h_support import *
def test_actual_petsc_cuda12_build_acceptance():
    # An observed failed hard gate remains a real failure, never skipped or xfailed.
    petsc_cuda_build_gate(load('petsc_cuda12_build'))
def test_source_unchanged_and_failure_diagnosis_is_preserved():
    a=actual('petsc_failure_audit')
    assert a['source_files_verified']==11035 and not a['changed_original_files']
    assert a['source_unmodified'] and not a['source_patch_applied'] and not a['PETSc_upgraded']
    assert not a['usable_library_exists']
    d=load('failure_classification')
    assert d['reason']=='PETSC_CUDA12_BUILD_FAIL'
    assert d['failure_class']=='CUDA12_THRUST_TUPLE_API_INCOMPATIBILITY'
    assert not d['previous_API_errors_recurred']
    assert any('has no member "get"' in e['text'] for e in d['compiler_errors'])
    assert sha256(ROOT/d['error_log'])==d['error_log_sha256']
def test_successful_compile_without_cuda_is_rejected():
    d=dict(configure_exit=0,make_exit=0,cuda_enabled=False,library_sha256='a'*64,source_unmodified=True)
    with pytest.raises(GateError,match='PETSC_CUDA_NOT_ENABLED'):petsc_cuda_build_gate(d)

