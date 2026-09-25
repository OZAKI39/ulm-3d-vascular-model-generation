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
