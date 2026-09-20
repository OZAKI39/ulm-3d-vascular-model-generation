from sv13g_support import *
def test_actual_cuda_configure():
    d=actual('cuda_configure');assert d['exit_code']==0 and d['MPI_verified_before_configure']
    assert d['version']==load('reference_manifest')['build']['PETSc_version']
def test_build_failure_is_specifically_classified():
    d=load('cuda_build')
    if d['status'] in ('FAIL','BLOCKED'):assert d['reason'] and d['failure_layer'] and d['diagnostic_evidence']
def test_observed_api_errors_support_toolkit_classification():
    d=load('cuda_build')
    if d['reason']!='CUDA_TOOLKIT_COMPATIBILITY':pytest.skip('No toolkit compatibility blocker')
    txt=(ROOT/d['diagnostic_evidence']['path']).read_text()
    assert "has no member named 'clockRate'" in txt
    assert "has no member named 'memoryClockRate'" in txt
    assert 'namespace "thrust" has no member class "unary_function"' in txt
    assert d['configure']=='PASS' and d['make']=='FAIL'
    assert load('cuda_toolkit_inventory')['usable_cuda12_found'] is False
def test_native_preconfigure_smoke_actually_passed():
    d=actual('preconfigure_smoke');assert len(d['tests'])==6
    assert all(r['exit_code']==0 and not r['timeout'] for r in d['tests'])
