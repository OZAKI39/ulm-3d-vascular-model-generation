from copy import deepcopy
import pytest
from stage03_helpers import *
from fem3d.vascular import validate_reference

def test_inlet_has_global_flow_constraint_without_profile_or_imposed_split():
    bc=preflight()['boundary_conditions']
    assert bc['inlet_velocity_profile'] is None and bc['outlet_velocity_profile'] is None
    for key in ('inlet_velocity_profile','outlet_flow_split'):
        c=deepcopy(config());c['boundary_model'][key]='imposed'
        with pytest.raises(ValueError):validate_reference(c)
