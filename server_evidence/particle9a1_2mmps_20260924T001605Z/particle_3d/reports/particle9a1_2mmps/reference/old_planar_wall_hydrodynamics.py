"""P9-A: one nearest planar wall, spherical MB, SI, reciprocal tangent block.

The legacy four-decimal constants permit an antisymmetric numerator of 8.75e-5;
the observed 5e-5 is projected to its reciprocal mean. No normal mobility is
imported. Geometric gaps are never changed by the logarithm evaluation floor.
"""
import numpy as np

XI_MIN=.001
XI_NEAR=.1
XI_FAR=1.
SHEAR_FORCE_COEFFICIENT=1.7005
SHEAR_TORQUE_COEFFICIENT=.9440
ROUNDING_NUMERATOR_BOUND=.5e-4*(.75+1.)


def stable_wall_frame(normal):
    n=np.asarray(normal,dtype=float)
    if n.shape!=(3,) or not np.isfinite(n).all() or np.linalg.norm(n)==0:
        raise ValueError('Finite nonzero wall normal required')
    n=n/np.linalg.norm(n)
    axis=np.eye(3)[np.argmin(np.abs(n))]
    t1=axis-n*np.dot(axis,n);t1/=np.linalg.norm(t1)
    return np.column_stack((t1,np.cross(n,t1),n))


def local_shear_vector(gradient,normal):
    g=np.asarray(gradient,dtype=float);n=stable_wall_frame(normal)[:,2]
    if g.shape!=(3,3) or not np.isfinite(g).all():raise ValueError('Finite 3x3 FEM gradient required')
    return (np.eye(3)-np.outer(n,n))@g@n


def wall_weight(xi):
    if not np.isfinite(xi):raise ValueError('Finite authoritative gap ratio required')
    if xi<=XI_NEAR:return 1.
    if xi>=XI_FAR:return 0.
    s=(xi-XI_NEAR)/(XI_FAR-XI_NEAR)
    return 1.-3*s*s+2*s*s*s


def legacy_planar_coefficients_si(gap_ratio,*,freeze_asymptote=True):
    """Dimensionless entries; dimensionalization is done using SI mu and radius."""
    xi=max(float(gap_ratio),XI_MIN)
    if freeze_asymptote:xi=min(xi,XI_NEAR)
    if not np.isfinite(xi):raise ValueError('Finite coefficient argument required')
    log=np.log(xi)
    c_ft=(8/15)*log-.9588;c_tt=.1*log-.1895
    c_fr=(2/15)*log-.2526;c_tr=.4*log-.3187
    det=c_tt*c_fr-c_ft*c_tr
    residual=.75*c_fr-c_tt
    if abs(residual)>ROUNDING_NUMERATOR_BOUND+32*np.finfo(float).eps:
        raise ValueError('Legacy coupling differs beyond published decimal rounding')
    tr=.75*c_fr/det;rt=c_tt/det
    return dict(coefficient_gap_ratio=xi,m_tt=c_tr/det,m_tR=tr,m_Rt=rt,
        m_cross=.5*(tr+rt),m_rr=.75*c_ft/det,determinant=det,
        reciprocity_numerator=residual,rounding_numerator_bound=ROUNDING_NUMERATOR_BOUND,
        projection_entry_error=abs(tr-rt)/2)


def local_mobility_4d(radius_m,mu,gap_ratio):
    if not np.isfinite([radius_m,mu]).all() or radius_m<=0 or mu<=0:
        raise ValueError('Positive finite SI radius and viscosity required')
    c=legacy_planar_coefficients_si(gap_ratio);w=wall_weight(gap_ratio)
    # C v = v cross n; the X-Z tangent convention t cross n = +Y fixes signs.
    cross_n=np.array([[0.,1.],[-1.,0.]])
    near=np.block([[c['m_tt']*np.eye(2),c['m_cross']*cross_n.T],
                   [c['m_cross']*cross_n,c['m_rr']*np.eye(2)]])
    free=np.diag([1.,1.,.75,.75]);base=1/(6*np.pi*mu*radius_m)
    return base*((1-w)*free+w*near),base*free,base*near,c


def planar_wall_affine_block(radius_m,mu,gap_m,normal,gradient,free_velocity):
    a=float(radius_m);h=float(gap_m);u=np.asarray(free_velocity,dtype=float)
    if u.shape!=(6,) or not np.isfinite(u).all():raise ValueError('Finite six-component free velocity required')
    frame=stable_wall_frame(normal);t=frame[:,:2];n=frame[:,2]
    xi=h/a;w=wall_weight(xi)
    effective,free,near,c=local_mobility_4d(a,mu,xi)
    L=np.zeros((4,6));L[:2,:3]=t.T;L[2:,3:]=a*t.T
    shear=local_shear_vector(gradient,n)
    force=SHEAR_FORCE_COEFFICIENT*6*np.pi*mu*a*(a+max(h,0.))*shear
    torque=SHEAR_TORQUE_COEFFICIENT*4*np.pi*mu*a**3*np.cross(shear,n)
    load=np.r_[t.T@force,t.T@torque/a]
    bulk=L@u;wall=near@load;target=(1-w)*bulk+w*wall
    # Only the local 4x4 block is inverted. No inverse of a 6N matrix.
    r_eff=np.linalg.solve(effective,np.eye(4));r_free=np.linalg.solve(free,np.eye(4))
    dr4=r_eff-r_free;db4=r_eff@target-r_free@bulk
    if w==0:dr4[:]=0.;db4[:]=0.
    dr=L.T@dr4@L;db=L.T@db4
    diag=dict(gap_m=h,gap_ratio=xi,wall_weight=w,wall_normal_xyz=n.tolist(),
        shear_vector_xyz=shear.tolist(),shear_force_N=force.tolist(),shear_torque_N_m=torque.tolist(),
        wall_normal_shear_force_N=float(force@n),m_tt=c['m_tt'],m_cross=c['m_cross'],m_rr=c['m_rr'],
        effective_m_tt=float(effective[0,0]/free[0,0]),effective_m_rr=float(effective[2,2]/free[0,0]),
        effective_m_cross=w*c['m_cross'],bulk_velocity_xyz=u[:3].tolist(),bulk_omega_xyz=u[3:].tolist(),
        target_tangential_velocity_xyz=(t@target[:2]).tolist(),
        target_tangential_omega_xyz=(t@target[2:]/a).tolist(),symmetry_error=float(np.max(abs(dr-dr.T))),
        coefficient_gap_ratio=c['coefficient_gap_ratio'],projection_entry_error=c['projection_entry_error'])
    return dr,db,diag
