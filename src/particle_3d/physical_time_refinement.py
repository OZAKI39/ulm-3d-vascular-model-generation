"""Conservative interval certificates and real-time binary subdivision.

Every accepted child advances its own time; the right child starts from the
completed left state. Failed trials never consume time or alter a center.
"""
from dataclasses import asdict,dataclass
import numpy as np
from numpy.polynomial import Polynomial as P
from .particle_shapes import Sphere,Ellipsoid,Capsule,EPS,roundoff_length
from .convex_triangle import triangle_gap,capsule_triangle_many

MAX_REFINEMENT_DEPTH=48  # integer implementation guard; not a production dt


class PhysicalTimeRefinementError(RuntimeError):
    def __init__(self,reason,interval,gap=None,triangle_id=None,shape_mode=None,state=None):
        self.record=dict(reason=reason,interval=list(interval),gap_m=gap,triangle_id=triangle_id,shape_mode=shape_mode,state=state)
        super().__init__(str(self.record))


class TrialNeedsSubdivision(RuntimeError):
    def __init__(self,reason,gap=None,triangle_id=None):
        self.reason=reason;self.gap=gap;self.triangle_id=triangle_id
        super().__init__(reason)


def refine_interval(state,end_time,trial,*,max_depth=MAX_REFINEMENT_DEPTH,depth=0,ledger=None,event_count=0):
    """trial(state,dt) returns a certified new state, or requests a true split."""
    ledger=[] if ledger is None else ledger
    begin=float(state.time_s);end=float(end_time)
    if not np.isfinite([begin,end]).all() or end<=begin:raise ValueError('strictly increasing finite physical time required')
    try:
        new,proof=trial(state,end-begin)
        terminal=getattr(new,'boundary_event','ACTIVE').startswith('OUTLET_')
        handoff=getattr(new,'projection',{}).get('handoff_event')
        if handoff and not terminal and begin<new.time_s<end:
            ledger.append(dict(t0_s=begin,t1_s=float(new.time_s),dt_s=float(new.time_s)-begin,
                requested_t1_s=end,depth=depth,accepted=True,proof=proof,handoff_event=handoff))
            if event_count>=max_depth:
                raise PhysicalTimeRefinementError('HANDOFF_EVENT_COUNT_LIMIT',(new.time_s,end),state=new.to_dict())
            return refine_interval(new,end,trial,max_depth=max_depth,depth=depth,ledger=ledger,event_count=event_count+1)
        if new.time_s!=end and not terminal:
            if abs(new.time_s-end)>4*np.spacing(max(abs(end),np.finfo(float).tiny)):
                raise ValueError('trial did not cover requested physical time')
        actual_end=float(new.time_s) if terminal else end
        ledger.append(dict(t0_s=begin,t1_s=actual_end,dt_s=actual_end-begin,requested_t1_s=end,depth=depth,accepted=True,proof=proof))
        return new,ledger
    except TrialNeedsSubdivision as error:
        mid=begin+(end-begin)/2
        if depth>=max_depth or mid<=begin or mid>=end:
            snapshot=state.to_dict() if hasattr(state,'to_dict') else repr(state)
            raise PhysicalTimeRefinementError('DEPTH_LIMIT' if depth>=max_depth else 'TIME_ULP_FLOOR',
                    (begin,end),error.gap,error.triangle_id,getattr(state,'shape_mode',None),snapshot) from error
        left,_=refine_interval(state,mid,trial,max_depth=max_depth,depth=depth+1,ledger=ledger,event_count=event_count)
        if getattr(left,'boundary_event','ACTIVE').startswith('OUTLET_'):return left,ledger
        return refine_interval(left,end,trial,max_depth=max_depth,depth=depth+1,ledger=ledger,event_count=event_count)


def capsule_plane_minimum(start,end,normal,triangle):
    """Continuous support-plane minimum for a volume-preserving capsule path.

Center and radius interpolate linearly; L(R) is recomputed from fixed volume.
Axis is normalized linear interpolation of endpoint axes. Stationary points
come from a degree <=12 polynomial, with all real roots, endpoints and the
abs(dot(axis,n)) kink checked. This is not endpoint-only collision detection.
"""
    n=np.asarray(normal);scale=max(start.bounding_radius_m,end.bounding_radius_m)
    r=P([start.radius_m/scale,(end.radius_m-start.radius_m)/scale])
    a=start.axis_world;b=end.axis_world-a
    k=P([a@a,2*(a@b),b@b]);m=P([n@a,n@b])
    if k(.5)<=128*EPS:raise TrialNeedsSubdivision('Capsule axis interpolation approaches antipodal singularity')
    volume=start.volume_m3/scale**3
    if abs(end.volume_m3/start.volume_m3-1)>2048*EPS:raise ValueError('capsule interval must preserve RBC volume')
    A=volume/np.pi;B=P([A/2])-(2/3)*r**3;D=-A*r.coef[1]-(2/3)*r.coef[1]*r**3
    dc=float(n@(end.center_m-start.center_m)/scale);C=dc-r.coef[1]
    numerator=(D*m+B*r*m.deriv())*k-.5*B*r*m*k.deriv()
    polynomial=numerator if abs(C)<=64*EPS*max(1.,abs(dc),abs(r.coef[1])) else C*C*r**6*k**3-numerator**2
    coefficients=polynomial.coef
    if np.max(np.abs(coefficients))>0:
        polynomial=P(coefficients/np.max(np.abs(coefficients))).trim(64*EPS)
        roots=polynomial.roots()
    else:roots=[]
    times=[0.,1.]
    # Nearly real conjugate pairs occur for a squared derivative root. Their
    # real parts are also evaluated; value error is second order at stationarity.
    for root in roots:
        if abs(np.imag(root))<=np.sqrt(EPS)*max(1.,abs(root)) and 0<float(np.real(root))<1:
            times.append(float(np.real(root)))
    if len(m.coef)>1 and m.coef[1]!=0:
        kink=-m.coef[0]/m.coef[1]
        if 0<kink<1:times.append(float(kink))
    support_triangle=float(np.max((np.asarray(triangle)-start.center_m)@n)/scale)
    values=[]
    for t in times:
        radius=r(t);half_length=A/(2*radius**2)-2*radius/3
        values.append((dc*t-radius-half_length*abs(m(t))/np.sqrt(k(t))-support_triangle)*scale)
    return float(min(values))


def swept_clearance_certificate(start,end,wall,*,rotation_angle=0.):
    """Certify every possible swept WALL triangle; false means subdivide."""
    if type(start) is not type(end):return False,'MODE_CHANGE_REQUIRES_RECOMPUTED_PATH',None
    radius=max(start.bounding_radius_m,end.bounding_radius_m)
    ids=wall.swept_candidates(start.center_m,end.center_m,radius)
    pad=roundoff_length(start.center_m,end.center_m,radius,wall.triangles)
    translation=np.linalg.norm(end.center_m-start.center_m)
    if isinstance(start,Capsule):
        angle=np.arccos(np.clip(start.axis_world@end.axis_world,-1.,1.))
        movement=translation+abs(end.radius_m-start.radius_m)+abs(end.cylindrical_length_m-start.cylindrical_length_m)/2+radius*angle
        first=capsule_triangle_many(start,wall.triangles[ids]);last=capsule_triangle_many(end,wall.triangles[ids])
        fixed=bool(np.array_equal(start.axis_world,end.axis_world))
    elif isinstance(start,Sphere):
        movement=translation;fixed=True
        cs=Capsule(start.center_m,[0,0,1],start.radius_m,0.)
        ce=Capsule(end.center_m,[0,0,1],end.radius_m,0.)
        first=capsule_triangle_many(cs,wall.triangles[ids]);last=capsule_triangle_many(ce,wall.triangles[ids])
    else:
        movement=translation+radius*abs(rotation_angle);fixed=bool(np.array_equal(start.rotation,end.rotation))
        first=last=None
    for j,i in enumerate(ids):
        tri=wall.triangles[i]
        if first is None:
            a=triangle_gap(start,tri);b=triangle_gap(end,tri);g0,g1=a.gap_m,b.gap_m;normals=[a.normal_inward,b.normal_inward]
        else:g0,g1=first[0][j],last[0][j];normals=[first[3][j],last[3][j]]
        if min(g0,g1)>movement+pad:continue
        if min(g0,g1)<-pad:return False,'PENETRATING_ENDPOINT',int(i)
        safe=False
        for n in normals:
            if isinstance(start,Capsule) and not fixed:
                minimum=capsule_plane_minimum(start,end,n,tri)
            else:
                # For fixed orientation h_capsule(R) is convex in R, so both
                # endpoint separating-plane bounds also certify any volume-
                # preserving radius interpolation. Fixed rigid translation is affine.
                minimum=min(float(n@(s.support(-n)-tri[0])-np.max((tri-tri[0])@n)) for s in [start,end])
                if not fixed:minimum-=radius*abs(rotation_angle)
            if minimum>=-pad:
                safe=True;break
        if not safe:return False,'SWEEP_NOT_CERTIFIED',int(i)
    return True,'CONTINUOUS_SUPPORT_OR_SURFACE_MOTION_BOUND',None


def quasistatic_capsule_path_certificate(start,end,wall,area_radius_min):
    """Certify existence of the instantaneous maximum-feasible-radius path.

The model radius at time t is the largest feasible R, NOT a linear interpolation
of the two endpoint maximizers. A volume/area-valid auxiliary capsule path
which clears WALL at every t proves the feasible set is nonempty throughout;
the maximum-feasible choice is then clear by its defining constraint. This
proof-only witness never replaces the accepted particle and never changes V.
Its rotating axis has the same normalized-linear local-direction path.

The witness search is only a sufficient certificate, not an infeasibility
decision. Failure requests real time subdivision. No coarse radius grid is
used to claim the actual maximal radius (that uses the separate bounded solver).
"""
    clear,proof,triangle=swept_clearance_certificate(start,end,wall)
    if clear:return clear,'QUASI_STATIC_FEASIBLE_PATH; '+proof,triangle
    volume=start.volume_m3;upper=min(start.radius_m,end.radius_m)
    lower=float(area_radius_min)
    for count in range(16):  # proof search safety guard only
        radius=upper-(upper-lower)/(2**(count+1))
        length=volume/(np.pi*radius**2)-4*radius/3
        a=Capsule(start.center_m,start.axis_world,radius,length)
        b=Capsule(end.center_m,end.axis_world,radius,length)
        safe,_,failed=swept_clearance_certificate(a,b,wall)
        if safe:
            return True,f'QUASI_STATIC_MAX_FEASIBLE_RADIUS_EXISTENCE; PROOF_ONLY_WITNESS_R_M={radius:.17g}',None
        triangle=failed
    return False,'QUASI_STATIC_FEASIBILITY_NOT_CERTIFIED; SUBDIVIDE_PHYSICAL_TIME',triangle
