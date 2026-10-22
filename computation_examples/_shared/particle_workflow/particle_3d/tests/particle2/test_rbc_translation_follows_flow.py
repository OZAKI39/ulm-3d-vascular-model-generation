import numpy as np
from particle_3d.rbc import RBCState
from particle_3d.rbc_integrator import advance_single_rbc
from particle_3d.particle1_cases import AffineValidationField
from particle_3d.rbc_orientation import quaternion_from_short_axis


def test_translation_uses_old_velocity_and_does_not_mutate(geometries):
    field=AffineValidationField(np.array([2e-4,-1e-4,3e-4]),np.array([[2,3,0],[-1,0,2],[1,1,-2]],float))
    for geometry in geometries:
        old=RBCState(0,[1e-5,2e-5,3e-5],quaternion_from_short_axis([1,2,3]),[10,20,30],[0,0,0],geometry)
        saved=old.position_m.copy();q=old.quaternion_wxyz.copy()
        new=advance_single_rbc(old,field,.001)
        np.testing.assert_array_equal(new.position_m,saved+.001*field.sample(saved).velocity_m_s)
        np.testing.assert_array_equal(new.velocity_m_s,field.sample(new.position_m).velocity_m_s)
        np.testing.assert_array_equal(old.position_m,saved);np.testing.assert_array_equal(old.quaternion_wxyz,q)
        assert not np.shares_memory(old.position_m,new.position_m)


def test_uniform_translation_roundoff_budget(static_case_result):
    rows,_=static_case_result
    bound=16*81*np.finfo(float).eps*1e-4
    assert all(r["position_error_m"]<=bound for r in rows)
