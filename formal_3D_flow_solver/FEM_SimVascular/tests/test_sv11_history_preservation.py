import gzip,json
from pathlib import Path
import pytest
from sv_validation.sv11 import load,REPORT,ROOT,require_file_hash,require_frozen_inputs
from sv_validation.validation import ValidationError
from sv_validation.provenance import sha256

def test_history_and_old_fem_audit():
    data=load('preservation_audit');assert data['status']=='PASS'
    assert data['old_fem_git_unchanged']
    original=json.loads(gzip.decompress((REPORT/'history_baseline.json.gz').read_bytes()))
    for path,item in original['files'].items():
        if 'sha256' in item:assert sha256(ROOT/path)==item['sha256']

def test_content_changed_mesh_rejected(tmp_path):
    frozen=ROOT/'outputs/sv1/SV_MESH/mesh-complete.mesh.vtu'
    data=bytearray(frozen.read_bytes());data[-1]^=1
    changed=tmp_path/'changed_mesh.vtu';changed.write_bytes(data)
    with pytest.raises(ValidationError,match='Frozen input'):require_file_hash(changed,sha256(frozen))

def test_actual_all_frozen_inputs_unchanged():
    assert require_frozen_inputs()
