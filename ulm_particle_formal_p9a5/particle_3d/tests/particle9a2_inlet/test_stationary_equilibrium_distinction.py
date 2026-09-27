import numpy as np
from particle_3d.particle9a2_equilibrium_audit import independent_contact_equilibrium


def test_feasible_fem_direction_does_not_disprove_driven_stationarity():
    rng=np.random.default_rng(842);a=rng.normal(size=(6,6));R=a.T@a+np.eye(6)
    J=np.column_stack([np.eye(3),np.zeros((3,3))])
    U=np.array([0.,0.,0.,.2,-.3,.1]);lam=np.array([1.,2.,3.])
    b=R@U-J.T@lam
    assert np.all(J@np.array([1.,1.,1.,0.,0.,0.])>0)
    c=independent_contact_equilibrium(R,b,J,U)
    assert c['recorded_solution_verified'] and c['independent_translation_zero_within_roundoff']


def test_spurious_zero_translation_is_rejected():
    J=np.column_stack([np.eye(3),np.zeros((3,3))])
    c=independent_contact_equilibrium(np.eye(6),np.array([1.,1.,1.,0.,0.,0.]),J,np.zeros(6))
    assert not c['recorded_solution_verified']
    assert not c['independent_translation_zero_within_roundoff']
