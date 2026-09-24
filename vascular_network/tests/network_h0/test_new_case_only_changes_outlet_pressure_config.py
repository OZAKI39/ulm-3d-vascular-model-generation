import json
from pathlib import Path
import pytest
from conftest import DATA
from network_1d0d.audit import sha256
from network_1d0d.fem_h0_case import audit_xml_change

def test_production_inputs_and_only_outlet_values_changed(bc):
    diff=json.loads((DATA/'3d_case_config_diff.json').read_text())
    old,new=Path(diff['original_case']),Path(diff['new_case'])
    values=[r['pressure_cap_shifted_Pa'] for r in bc['ports']]
    fields=audit_xml_change((old/'run/solver.xml').read_bytes(),(new/'run/solver.xml').read_bytes(),values)
    assert sum(r['numerically_changed'] for r in fields)==2 # O3 remains gauge zero.
    for name,digest in diff['unchanged_input_hashes'].items():
        assert sha256(old/name)==sha256(new/name)==digest
    assert sha256(old/'run/solver.xml')==diff['original_solver_sha256']
    tampered=(new/'run/solver.xml').read_bytes().replace(b'<Density>1056.0',b'<Density>1000.0')
    with pytest.raises(ValueError,match='non-outlet-pressure'):
        audit_xml_change((old/'run/solver.xml').read_bytes(),tampered,values)
