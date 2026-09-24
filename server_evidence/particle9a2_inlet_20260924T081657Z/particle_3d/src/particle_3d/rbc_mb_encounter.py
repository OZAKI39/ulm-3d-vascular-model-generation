"""Exact P4 translational contact on one fixed capsule end cap.

This adapter changes time integration, not P4's kinematic contact metric.
It is intentionally restricted to a certified hemispherical contact region,
with no wall contact or capsule rotation. Unsupported scenes fail explicitly.
"""
import numpy as np
from .particle_shapes import Sphere,Capsule,roundoff_length
from .particle4_motion import isolated_sphere_contact_path
from .pair_geometry import ParticlePairGapResult


def capsule_sphere_gap(cap,bubble,i=101,j=203):
    d=bubble.center_m-cap.center_m
    q=np.clip(d@cap.axis_world,-cap.cylindrical_length_m/2,cap.cylindrical_length_m/2)
    tip=cap.center_m+q*cap.axis_world;delta=tip-bubble.center_m;distance=np.linalg.norm(delta)
    if distance==0:raise ValueError('Overlapping centerline is outside encounter domain')
    n=delta/distance;gap=distance-cap.radius_m-bubble.radius_m
    pad=roundoff_length(cap.center_m,bubble.center_m,cap.bounding_radius_m,bubble.radius_m)
    state='SEPARATED' if gap>pad else 'PENETRATING' if gap< -pad else 'TOUCHING'
    return ParticlePairGapResult(i,j,cap.mode,bubble.mode,float(gap),state,
        tip-cap.radius_m*n,bubble.center_m+bubble.radius_m*n,n,abs(gap),0,pad,(i,j))


def advance_encounter(cap,bubble,free_cap,free_bubble,dt):
    """Frozen-free-velocity substep: free flight, exact hit, curved contact, release.

Output motion can be sampled at any partial dt by calling this same function.
No center correction, artificial shell, effective RBC sphere, or force law.
The sphere at the end cap is only the exact local capsule surface geometry.
"""
    if dt<0:raise ValueError('Nonnegative duration required')
    if not isinstance(cap,Capsule) or not isinstance(bubble,Sphere):raise TypeError('Capsule and sphere required')
    vc=np.asarray(free_cap);vb=np.asarray(free_bubble);axis=cap.axis_world
    relative=bubble.center_m-cap.center_m;half=cap.cylindrical_length_m/2
    side=1. if relative@axis>=0 else -1.
    # Relative speed under the P4 orthogonal velocity projection cannot grow.
    # This certifies that no surface point crosses into the cylindrical regime.
    margin=side*(relative@axis)-half
    if margin<=np.linalg.norm(vb-vc)*dt:raise ValueError('CAP_HEMISPHERE_CERTIFICATE_FAILED')
    gap=capsule_sphere_gap(cap,bubble)
    if gap.state=='PENETRATING':raise ValueError('Initial overlap')
    tip=cap.center_m+side*half*axis;d=tip-bubble.center_m;w=vc-vb
    radius=cap.radius_m+bubble.radius_m
    hit=None
    if gap.state=='TOUCHING' and d@w<0:hit=0.
    elif d@w<0 and w@w>0:
        c=d@d-radius**2;dw=d@w;disc=dw*dw-(w@w)*c
        if disc>=0:
            root=c/(-dw+np.sqrt(disc))
            if 0<=root<=dt:hit=float(root)
    contact_duration=0.;vel={101:vc.copy(),203:vb.copy()}
    if hit is None:
        end_cap=cap.moved(cap.center_m+dt*vc);end_mb=bubble.moved(bubble.center_m+dt*vb)
    else:
        at_cap=cap.moved(cap.center_m+hit*vc);at_mb=bubble.moved(bubble.center_m+hit*vb)
        virtual={101:Sphere(at_cap.center_m+side*half*axis,cap.radius_m),203:at_mb}
        exact=isolated_sphere_contact_path(virtual,{101:vc,203:vb},dt-hit)
        if exact is None:
            end_cap=at_cap.moved(at_cap.center_m+(dt-hit)*vc);end_mb=at_mb.moved(at_mb.center_m+(dt-hit)*vb)
        else:
            end,vel,record=exact;end_cap=at_cap.moved(end[101].center_m-side*half*axis);end_mb=end[203]
            contact_duration=float(record['contact_duration_s'])
    # P4 translation retains the exact axis bits.
    object.__setattr__(end_cap,'axis_world',cap.axis_world)
    end_gap=capsule_sphere_gap(end_cap,end_mb)
    if end_gap.state=='PENETRATING':raise ValueError('Analytic encounter overlap')
    record=dict(method='EXACT_P4_FIXED_CAPSULE_HEMISPHERE_CONTACT',dt_s=dt,hit_time_s=hit,
        contact_duration_s=contact_duration,hemisphere_side=side,hemisphere_margin_m=margin,
        relative_motion_bound_m=float(np.linalg.norm(w)*dt),position_projection=False,
        cap_free_velocity_m_s=vc.tolist(),mb_free_velocity_m_s=vb.tolist(),
        pair_gap_m=end_gap.gap_m,roundoff_budget_m=end_gap.roundoff_budget_m,
        contact_active_at_endpoint=bool(contact_duration>0 and end_gap.state=='TOUCHING'))
    return end_cap,end_mb,vel,record
