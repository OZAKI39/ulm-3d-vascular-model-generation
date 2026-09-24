import numpy as np
import pytest
from particle_3d.microbubble import MicrobubbleState,FieldUnavailableError
from particle_3d.integrator import advance_single_microbubble
from particle_3d.particle1_cases import AffineValidationField,single_step_case
from particle_3d.field import FrozenFEMField


def test_euler_uses_old_position_field_and_refreshes_endpoint_state():
    g=np.diag([1.,2.,3.]); field=AffineValidationField(np.array([2e-4,-1e-4,1e-4]),g)
    state=MicrobubbleState(0,[1e-5,2e-5,3e-5],1e-6,[100.,200.,300.],[0.,0.,0.])
    new=advance_single_microbubble(state,field,.001)
    expected=state.position_m+.001*field.sample(state.position_m).velocity_m_s
    np.testing.assert_array_equal(new.position_m,expected)
    np.testing.assert_array_equal(new.velocity_m_s,field.sample(expected).velocity_m_s)
    assert new is not state and new.particle_id==state.particle_id and new.radius_m==state.radius_m


def test_outside_endpoint_rejected_without_returning_invalid_active_state():
    vertices=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    field=FrozenFEMField(vertices,np.array([[0,1,2,3]]),np.tile([1.,0.,0.],(4,1)),np.zeros(4))
    state=MicrobubbleState(0,[.1,.1,.1],1e-6,[0.,0.,0.],[0.,0.,0.])
    with pytest.raises(FieldUnavailableError): advance_single_microbubble(state,field,2.)
