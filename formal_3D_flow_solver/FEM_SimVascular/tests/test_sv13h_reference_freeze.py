from sv13h_support import *
def test_actual_frozen_reference_checksums():
    d=actual('reference_manifest')
    assert len(d['files'])>=45 and d['CPU_production_read_only']
    for f in d['files']: assert sha256(ROOT/f['path'])==f['sha256'], f['path']
def test_versions_and_production_are_frozen():
    d=load('reference_manifest')
    assert d['build']['PETSc_version']=='3.19.6'
    assert d['build']['commit']=='c3f0bb892b765b718f61069ecd9726dbc6d177fd'
    s=load('stage_result')
    assert not s['production_changed'] and s['production']=='CPU_EARLY_STOP_PRODUCTION'
    assert not s['mesh_convergence_started'] and not s['GPU_full_production_started']

