from sv13j_support import *
def test_frozen_references_match_actual_files():
    d=actual('reference_manifest')
    assert len(d['files'])==53 and d['CPU_production_read_only']
    for f in d['files']:assert sha256(ROOT/f['path'])==f['sha256'],f['path']
    assert d['build']['commit']=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
def test_protected_remote_environment_unchanged():
    d=actual('environment_preservation');assert all(d['checks'].values())
    assert not load('stage_result')['production_changed']

