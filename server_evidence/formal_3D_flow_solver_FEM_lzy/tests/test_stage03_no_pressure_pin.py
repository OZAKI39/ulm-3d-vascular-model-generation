from copy import deepcopy
import pytest
from stage03_helpers import *
from fem3d.vascular import validate_reference

def test_pressure_has_only_natural_traction_reference():
    bc=preflight()['boundary_conditions']
    assert bc['pressure_dirichlet_dofs']==0 and not bc['pressure_nullspace_registered']
    for key in ('pressure_pin','pressure_dirichlet'):
        c=deepcopy(config());c['boundary_model'][key]=True
        with pytest.raises(ValueError):validate_reference(c)
