from sv13j_support import *
def test_official_side_by_side_gcc12_and_smokes():
    d=actual('gcc12_environment')
    assert d['default_compilers_unchanged'] and not d['unsupported_compiler_override']
    for n in ('gcc-12','g++-12'):assert d['compilers'][n]['number'].split('.')[0]=='12'
    for s in d['smokes']:assert s['compile']['exit_code']==s['run']['exit_code']==0 and 'result=42' in s['run']['stdout']
    assert len(d['libstdcpp_sha256'])==64
@pytest.mark.parametrize('major,flags',[(13,[]),(12,['--allow-unsupported-compiler']),(12,['-allow-unsupported-compiler'])])
def test_unsupported_host_compiler_and_override_rejected(major,flags):
    with pytest.raises(GateError):
        toolchain_gate('12.3.2',major,'release 12.3','/cuda123/bin/nvcc','/cuda123',flags)

