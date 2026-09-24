import numpy as np
from particle_3d.microbubble import MicrobubbleState
from particle_3d.integrator import advance_single_microbubble
from particle_3d.particle1_cases import AffineValidationField


def test_step_leaves_every_input_field_unchanged():
    state=MicrobubbleState(7,[1.,2.,3.],1e-6,[4.,5.,6.],[7.,8.,9.])
    before=[state.position_m.tobytes(),state.velocity_m_s.tobytes(),state.angular_velocity_s_inv.tobytes()]
    new=advance_single_microbubble(state,AffineValidationField(np.array([.1,.2,.3]),np.zeros((3,3))),.01)
    assert before==[state.position_m.tobytes(),state.velocity_m_s.tobytes(),state.angular_velocity_s_inv.tobytes()]
    assert state.radius_m==1e-6 and state.particle_id==7
    assert not np.shares_memory(state.position_m,new.position_m)
