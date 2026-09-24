import numpy as np
from particle_3d.particle9a_diagnostics import closure_metrics, decomposition

def test_linear_shear_has_velocity_units_and_scale_covariance():
    a,h,gamma=1e-6,2e-8,500.; g=np.zeros((3,3));g[0,2]=gamma
    u=np.array([gamma*(a+h),0,0.]); n=np.array([0.,0.,1.])
    d=closure_metrics(u,g,n,a,h)
    np.testing.assert_array_equal(d['tangential_bulk_velocity_xyz'],d['Hs_velocity_scale_xyz'])
    assert d['closure_relative_mismatch']==0 and d['closure_speed_ratio']==1
    e=closure_metrics(1000*u,g,n,1000*a,1000*h)
    np.testing.assert_allclose(e['linear_shear_speed'],1000*d['linear_shear_speed'])
    assert closure_metrics(np.zeros(3),np.zeros((3,3)),n,a,h)['closure_angle_deg'] is None
    c=decomposition(u,np.array([0.,gamma/2,0.]),g,n,a,h,.00345312)
    np.testing.assert_array_equal(c['q_target'],c['q_wall'])  # w == 1
