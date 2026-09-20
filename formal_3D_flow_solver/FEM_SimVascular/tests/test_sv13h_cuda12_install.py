from sv13h_support import *
def test_verified_official_source_used_for_toolkit_only():
    s=actual('cuda12_source_manifest');i=actual('cuda12_install')
    assert s['url'].startswith('https://developer.download.nvidia.com/compute/cuda/12.6.3/')
    assert s['actual_md5']==s['official_md5'] and len(s['sha256'])==64
    assert s['frozen_before_installer_execution'] and s['final_12_6_patch_verified']
    assert s['sha256']==i['source_sha256']
    assert i['exit_code']==0 and i['toolkit_only'] and not i['driver_install']
    assert '--toolkit' in i['command'] and '--driver' not in i['command']
    assert i['installer_uid']==65534
def test_driver_cuda13_and_linker_preservation():
    d=actual('final_cuda13_preservation')
    assert all(d['checks'].values())
    unchanged_driver_and_cuda13(load('pre_install_environment'),load('final_environment'))
def test_modified_driver_rejected():
    before=load('pre_install_environment');after=copy.deepcopy(before)
    key=next(iter(after['driver_files']));after['driver_files'][key]['sha256']='0'*64
    with pytest.raises(GateError,match='DRIVER_CHANGED'):unchanged_driver_and_cuda13(before,after)

