from types import SimpleNamespace
import numpy as np
import pytest
from particle_3d.particle82_checkpoint import pending_intervals,restore_stepper_recording,preceding_nominal_state


def test_pending_right_children_retain_original_depth_and_nominal_grid():
    # Last accepted leaf [1/4,3/8] of the original [0,1] step.
    pending,k=pending_intervals(np.array([[.25],[.375]]),1.)
    assert pending==[(.375,.5,3),(.5,1.,1)]
    assert k==2
    assert pending_intervals(np.array([[.5],[1.]]),1.)==([],2)
    with pytest.raises(ValueError):pending_intervals(np.array([[.2],[.375]]),1.)


def test_restore_keeps_every_saved_float_and_restores_endpoint():
    samples=np.zeros((3,20));samples[:,0]=[0.,.5,1.];samples[:,1]=[1.,2.,2.]
    meta=dict(accepted_steps=2,rejected_trials=3,maximum_refinement_depth=4,minimum_original_wall_gap_m=1e-9,
              maximum_position_identity_error_m=1e-20,nearfield_states={'FAR':3})
    step=SimpleNamespace(read=lambda:[SimpleNamespace(position=samples[-1,1:4])],time_s=1.)
    restore_stepper_recording(step,meta,samples)
    assert np.array(step.samples).tobytes()==samples.tobytes()
    assert step.accepted_count==2 and step.rejected_count==3 and step.maximum_depth==4
    position,unchanged=preceding_nominal_state(samples,.5)
    assert np.array_equal(position,samples[-1,1:4]) and unchanged==1
    step.time_s=.5
    with pytest.raises(ValueError,match='exact saved endpoint'):restore_stepper_recording(step,meta,samples)
