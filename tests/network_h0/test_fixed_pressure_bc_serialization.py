import json
import numpy as np
import pytest
from mouse_vessel_flow.paths import project_path
DATA=project_path('inputs/network/reports/a_network_1d0d_boundary_v2_idealized/data')
from mouse_vessel_flow.network_models.fem_h0_case import serialize_fixed_pressure_xml,audit_xml_change

def test_saved_json_to_solver_xml_roundtrip(bc):
    diff=json.loads((DATA/'3d_case_config_diff.json').read_text())
    from pathlib import Path
    old=(Path(diff['original_case'])/'run/solver.xml').read_bytes()
    values=[r['pressure_cap_shifted_Pa'] for r in bc['ports']]
    encoded=serialize_fixed_pressure_xml(old,values)
    fields=audit_xml_change(old,encoded,values)
    assert [r['new'] for r in fields]==values
    for bad in ([0,np.nan,0],[0,np.inf,0],[0,1]):
        with pytest.raises(ValueError): serialize_fixed_pressure_xml(old,bad)
