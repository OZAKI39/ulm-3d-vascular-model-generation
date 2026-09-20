from sv_validation.sv11 import load
from sv_validation.provenance import sha256

def test_real_petsc_binary_and_dependency():
    data=load('petsc_build_manifest')
    assert data['status']=='PASS'
    assert data['commit']=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
    assert sha256(data['executable'])==data['executable_sha256']
    assert 'libpetsc.so' in data['linked_libraries'] and 'not found' not in data['linked_libraries']
    assert not data['additional_algebra_packages']
    assert data['source_clean'] and not data['mesh_library_recompiled']
