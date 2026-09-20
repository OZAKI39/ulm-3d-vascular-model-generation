import pytest
from sv_validation.sv12 import ROOT,load,checkpoint_audit
from sv_validation.provenance import sha256
from sv_validation.validation import ValidationError

def test_actual_native_restart_has_full_state():
    decision=load('restart_decision');dt=load('frozen_input_manifest')
    import json
    policy=json.loads((ROOT/'configs/time_policy.json').read_text())
    audit=checkpoint_audit(decision['path'],10,policy['dt_s'])
    assert audit['sha256']==decision['source_checkpoint_sha256']
    assert all(r['Y_values']==r['A_values']>0 and r['finite'] for r in audit['records'])
    assert len(audit['records'])==4

def test_vtu_not_restart():
    with pytest.raises(ValidationError,match='VTU'):checkpoint_audit(ROOT/'outputs/sv1_1/vascular_short/4-procs/result_010.vtu',10,1.)

def test_truncated_checkpoint_rejected(tmp_path):
    source=load('restart_decision')['path'];p=tmp_path/'bad.bin'
    from pathlib import Path
    p.write_bytes(Path(source).read_bytes()[:-100])
    import json
    dt=json.loads((ROOT/'configs/time_policy.json').read_text())['dt_s']
    with pytest.raises(ValidationError):checkpoint_audit(p,10,dt)
