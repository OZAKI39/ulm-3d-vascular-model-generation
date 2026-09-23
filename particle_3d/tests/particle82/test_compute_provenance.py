import pytest
from particle_3d.particle82_provenance import require_remote,atomic_json

def test_authenticated_remote_contract(remote):
    assert require_remote(remote,'test-remote')

@pytest.mark.parametrize('field,value',[('role','DEVELOPMENT_ORCHESTRATION_ONLY'),('hostname','localhost'),('ssh_endpoint','unknown')])
def test_invalid_host_receipt_rejected(remote,field,value):
    remote[field]=value
    with pytest.raises(ValueError):require_remote(remote,'test-remote')

def test_wsl_cannot_be_relabelled_remote(remote,monkeypatch):
    monkeypatch.setattr('particle_3d.particle82_provenance.platform.release',lambda:'6.6-microsoft-WSL2')
    with pytest.raises(ValueError,match='forbidden on WSL'):require_remote(remote,'test-remote')

def test_atomic_json_preserves_old_file_on_invalid_payload(tmp_path):
    p=tmp_path/'receipt.json';atomic_json(p,{'old':True})
    with pytest.raises(ValueError):atomic_json(p,{'bad':float('nan')})
    assert p.read_text().find('old')>=0
