from dataclasses import fields
import inspect
import numpy as np
import pytest
from particle_3d.rbc import RBCState
from particle_3d.rbc_integrator import advance_single_rbc
from particle_3d.particle1_cases import AffineValidationField


def test_state_shapes_immutable_arrays_no_inertia(geometries):
    x=np.array([1e-5,2e-5,3e-5]);q=np.array([2.,0.,0.,0.])
    state=RBCState(1,x,q,[0,0,0],[0,0,0],geometries[0])
    x[:]=1;q[:]=1
    np.testing.assert_array_equal(state.position_m,[1e-5,2e-5,3e-5])
    np.testing.assert_array_equal(state.quaternion_wxyz,[1,0,0,0])
    for name in ["position_m","quaternion_wxyz","velocity_m_s","angular_velocity_s_inv"]:
        a=getattr(state,name)
        assert a.dtype==np.float64
        with pytest.raises(ValueError):a.setflags(write=True)
    assert {f.name for f in fields(RBCState)}=={"particle_id","position_m","quaternion_wxyz","velocity_m_s","angular_velocity_s_inv","geometry"}
    assert inspect.signature(advance_single_rbc).parameters["dt_s"].default is inspect.Parameter.empty


@pytest.mark.parametrize("q",[[0,0,0,0],[1,2,3],[np.nan,0,0,1]])
def test_invalid_quaternion_inputs(q,geometries):
    with pytest.raises(ValueError):RBCState(0,[0,0,0],q,[0,0,0],[0,0,0],geometries[0])


@pytest.mark.parametrize("dt",[0,-1,np.nan,np.inf,None,True])
def test_invalid_dt_rejected_without_mutating_state(dt,geometries):
    state=RBCState(0,[0,0,0],[1,0,0,0],[0,0,0],[0,0,0],geometries[0])
    field=AffineValidationField(np.ones(3),np.zeros((3,3)))
    with pytest.raises((ValueError,TypeError)):advance_single_rbc(state,field,dt)
    np.testing.assert_array_equal(state.position_m,[0,0,0])
