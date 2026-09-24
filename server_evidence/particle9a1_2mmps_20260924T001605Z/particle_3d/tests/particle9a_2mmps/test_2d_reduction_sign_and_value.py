import numpy as np
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block,legacy_planar_coefficients_si

def test_2d_reduction_sign_and_value(legacy):
    n=np.array([.6,0,.8]);t=np.array([-.8,0,.6]);a=1e-6;mu=.00345312;gamma=500.
    gradient=gamma*np.outer(t,n)
    for xi in [.003,.01,.1,.5,1.]:
        velocity=gamma*a*(1+xi)*t;omega=np.array([0,-gamma/2,0.])
        _,_,d=planar_wall_affine_block(a,mu,xi*a,n,gradient,np.r_[velocity,omega])
        vx,vz,wy,w,s=legacy.background_hydrodynamic_velocity_scalar(velocity[0]*1e6,velocity[2]*1e6,
            gradient[0,0],gradient[0,2],gradient[2,0],gradient[2,2],n[0],n[2],a*1e6,xi*a*1e6,mu)
        vt=np.array(d['target_tangential_velocity_xyz'])@t;wr=np.array(d['target_tangential_omega_xyz'])[1]
        oldvt=(vx*t[0]+vz*t[2])*1e-6
        c=legacy_planar_coefficients_si(xi);err=c['projection_entry_error'];base=1/(6*np.pi*mu*a)
        # Analytic projection error from the published 4-decimal constants,
        # not a fitted tolerance. Roundoff allowance is < 1e-12 relative.
        fs,ty,_=legacy.wall_shear_load_scalar(mu,a*1e6,xi*a*1e6,gamma)
        assert abs(vt-oldvt)<=w*base*err*abs(ty*1e-18/a)+1e-15
        assert abs(wr-wy)<=w*base*err*abs(fs*1e-12)/a+1e-9
        if xi<=.1:assert wr>0 and vt>0
    for xi in [.001,.003,.01,.03,.1,.3,1.]:
        c=legacy_planar_coefficients_si(xi,freeze_asymptote=False)
        ref=legacy.wall_dimensionless_mobility_entries_scalar(xi)
        np.testing.assert_allclose([c['m_tt'],c['m_tR'],c['m_Rt'],c['m_rr']],[ref[0],ref[2],ref[3],ref[4]],rtol=3e-15)
        assert abs(c['reciprocity_numerator'])<=c['rounding_numerator_bound']
