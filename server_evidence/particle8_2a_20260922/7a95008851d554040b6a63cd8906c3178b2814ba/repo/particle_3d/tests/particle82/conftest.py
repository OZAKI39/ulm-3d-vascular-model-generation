from pathlib import Path
import sys,pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'src'))

@pytest.fixture
def remote(monkeypatch):
    import particle_3d.particle82_provenance as p
    monkeypatch.setattr(p.socket,'gethostname',lambda:'test-remote')
    monkeypatch.setattr(p.platform,'release',lambda:'6.8-linux')
    return dict(hostname='test-remote',role='HEAVY_COMPUTE_PRIMARY',ssh_endpoint='root@50.115.148.16:4159',
                source_git_commit='test-commit',frozen_input_sha256={'input':'abc'})

@pytest.fixture(scope='session')
def evidence_root():
    import os
    value=os.environ.get('PARTICLE82_RESULTS_ROOT')
    if not value:pytest.skip('Set PARTICLE82_RESULTS_ROOT to the actual completed remote evidence tree')
    path=Path(value)
    assert path.is_dir()
    return path

@pytest.fixture(scope='session')
def point_root(evidence_root):
    import os
    return Path(os.environ.get('PARTICLE82_POINT_ROOT',str(evidence_root/'point_basin_100000')))

@pytest.fixture(scope='session')
def saved_scene(evidence_root):
    from particle_3d.particle82_results import Scene
    return Scene(evidence_root)
