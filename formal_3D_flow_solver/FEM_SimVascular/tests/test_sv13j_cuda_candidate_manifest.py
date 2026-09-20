from sv13j_support import *
def test_official_candidate_checksum_and_toolkit_only():
    d=actual('cuda123_source_manifest')
    assert d['url'].startswith('https://developer.download.nvidia.com/compute/cuda/12.3.2/')
    assert d['official_md5']==d['actual_md5'] and len(d['sha256'])==64
    assert d['frozen_before_installer_execution']
    i=actual('cuda123_installation')
    assert i['source_sha256']==d['sha256'] and i['toolkit_only'] and not i['driver_install']
    assert '--toolkit' in i['command'] and '--driver' not in i['command']
    assert not (R/'cuda122_source_manifest.json').exists() and not (R/'cuda121_source_manifest.json').exists()

