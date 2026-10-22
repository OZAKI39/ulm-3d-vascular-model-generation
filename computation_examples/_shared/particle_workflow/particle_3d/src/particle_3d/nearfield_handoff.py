"""Exact-geometry continuum handoff constraints and continuous certificates.

All shapes/triangles remain unchanged. A support-plane lower bound is compared
with h_lower, never with a physically inflated particle or displaced wall.
"""
import numpy as np
from .particle_shapes import Capsule,roundoff_length
from .convex_triangle import capsule_triangle_many,feature
from .pair_geometry import pair_gap
from .pair_broadphase import all_pairs
from .wall_gap import wall_gap,touching_contacts
from .kinematic_contact import ContactConstraint
from .hydrodynamic_resistance import PhysicalNearField,require_sphere


class InitialBelowContinuumHandoff(ValueError):
    def __init__(self,records):
        self.record=dict(status='INITIAL_STATE_BELOW_CONTINUUM_HANDOFF',interactions=records,
            action='REJECT_WITHOUT_MOVING_CENTERS_OR_CHANGING_GEOMETRY')
        super().__init__(str(self.record))


def exact_interactions(shapes,wall,policy,mu,pairs=None):
    specs=[];gaps={}
    if wall is not None:
        for i,s in sorted(shapes.items()):
            g=wall_gap(s,wall);gaps[('WALL',i)]=g
            specs.append(PhysicalNearField(i,None,g.gap_m,g.normal_inward,g.roundoff_m))
    for i,j in all_pairs(shapes) if pairs is None else pairs:
        g=pair_gap(shapes[i],shapes[j],i,j);gaps[('PAIR',i,j)]=g
        specs.append(PhysicalNearField(i,j,g.gap_m,g.normal_j_to_i,g.roundoff_budget_m))
    records=[policy.evaluate(s,shapes,mu)[2] for s in specs]
    return specs,records,gaps


def require_admissible_initial(shapes,wall,policy,mu):
    for s in shapes.values():require_sphere(s)
    _,rows,_=exact_interactions(shapes,wall,policy,mu)
    bad=[r for r in rows if not r['continuum_state_admissible']]
    if bad:raise InitialBelowContinuumHandoff(bad)
    return rows


def handoff_constraints(shapes,wall,policy):
    contacts=[];records=[]
    for i,j in all_pairs(shapes):
        g=pair_gap(shapes[i],shapes[j],i,j);lower=policy.lower_handoff_gap(policy.reference_length(shapes[i],shapes[j]))['h_lower_m']
        if g.state=='TOUCHING':
            contacts.append(ContactConstraint.pair(g));records.append(dict(contact_type='GEOMETRIC_HARD_CONTACT',pair=[i,j]))
        if abs(g.gap_m-lower)<=g.roundoff_budget_m:
            c=ContactConstraint(('CONTINUUM_HANDOFF_CONTACT','PAIR',i,j),i,j,g.normal_j_to_i,g.point_i_m,g.point_j_m)
            contacts.append(c);records.append(dict(contact_type='CONTINUUM_HANDOFF_CONTACT',canonical_id=c.canonical_id,
                h_geom_m=g.gap_m,h_lower_m=lower,g_nf_m=g.gap_m-lower,roundoff_m=g.roundoff_budget_m,normal=c.normal.tolist()))
    if wall is not None:
        for i,s in sorted(shapes.items()):
            for g in touching_contacts(s,wall):
                contacts.append(ContactConstraint.wall(i,g));records.append(dict(contact_type='GEOMETRIC_HARD_CONTACT',particle_id=i))
            lower=policy.lower_handoff_gap(policy.reference_length(s))['h_lower_m'];pad=roundoff_length(s.center_m,s.radius_m,wall.triangles)
            ids=wall.candidates(s.center_m,s.radius_m+lower)
            if not len(ids):continue
            cap=Capsule(s.center_m,[0,0,1],s.radius_m,0.)
            gs,wp,pp,ns,bary=capsule_triangle_many(cap,wall.triangles[ids])
            for k in np.flatnonzero(np.abs(gs-lower)<=pad):
                c=ContactConstraint(('CONTINUUM_HANDOFF_CONTACT','WALL',i,int(ids[k]),feature(bary[k])),i,None,ns[k],pp[k])
                contacts.append(c);records.append(dict(contact_type='CONTINUUM_HANDOFF_CONTACT',canonical_id=c.canonical_id,
                    h_geom_m=float(gs[k]),h_lower_m=lower,g_nf_m=float(gs[k]-lower),roundoff_m=pad,normal=c.normal.tolist(),
                    wall_point_m=wp[k].tolist(),particle_point_m=pp[k].tolist()))
    return contacts,records


def pair_handoff_certificate(a,b,end_a,end_b,lower):
    r=a.center_m-b.center_m;delta=(end_a.center_m-end_b.center_m)-r
    denom=float(delta@delta);t=0. if denom==0 else float(np.clip(-(r@delta)/denom,0.,1.))
    h=float(np.linalg.norm(r+t*delta)-a.radius_m-b.radius_m)
    pad=roundoff_length(a.center_m,b.center_m,end_a.center_m,end_b.center_m,a.radius_m,b.radius_m)
    return h-lower>=-pad,dict(proof='EXACT_LINEAR_SPHERE_PAIR_MINIMUM_WITH_G_NF',minimum_h_geom_m=h,
        h_lower_m=lower,minimum_g_nf_m=h-lower,roundoff_m=pad,fraction_of_minimum=t)


def wall_handoff_certificate(start,end,wall,lower):
    ids=wall.swept_candidates(start.center_m,end.center_m,start.radius_m+lower)
    pad=roundoff_length(start.center_m,end.center_m,start.radius_m,wall.triangles)
    if not len(ids):return True,dict(proof='NO_WALL_IN_HANDOFF_QUERY_ENVELOPE',triangle_count=0,roundoff_m=pad)
    first=capsule_triangle_many(Capsule(start.center_m,[0,0,1],start.radius_m,0.),wall.triangles[ids])
    last=capsule_triangle_many(Capsule(end.center_m,[0,0,1],end.radius_m,0.),wall.triangles[ids])
    distance=float(np.linalg.norm(end.center_m-start.center_m));minimum_bound=float('inf')
    for k,i in enumerate(ids):
        g0,g1=float(first[0][k]),float(last[0][k])
        if min(g0,g1)-lower < -pad:
            return False,dict(proof='HANDOFF_ENDPOINT_BELOW_LOWER',triangle_id=int(i),minimum_g_nf_m=min(g0,g1)-lower,roundoff_m=pad)
        if min(g0,g1)-distance>=lower-pad:
            minimum_bound=min(minimum_bound,min(g0,g1)-distance-lower);continue
        tri=wall.triangles[i];bounds=[]
        for n in [first[3][k],last[3][k]]:
            # Original sphere support and original finite triangle; subtract only the constraint gap.
            h0=float(n@(start.center_m-tri[0])-start.radius_m-np.max((tri-tri[0])@n))
            h1=float(n@(end.center_m-tri[0])-end.radius_m-np.max((tri-tri[0])@n))
            bounds.append(min(h0,h1)-lower)
        bound=max(bounds);minimum_bound=min(minimum_bound,bound)
        if bound < -pad:return False,dict(proof='HANDOFF_SWEEP_NOT_CERTIFIED',triangle_id=int(i),minimum_g_nf_bound_m=bound,roundoff_m=pad)
    return True,dict(proof='CONTINUOUS_ORIGINAL_SUPPORT_PLANE_MINUS_H_LOWER',triangle_count=len(ids),minimum_g_nf_bound_m=minimum_bound,roundoff_m=pad)


def first_handoff_event(shapes,velocities,dt,wall,policy):
    """First finite-geometry crossing along the solved straight center path.

    Distance to a convex finite triangle is convex. Its squared-distance
    derivative brackets the unique minimum (possibly an interval of minima).
    Checking that minimum catches crossings even when BOTH endpoints are clear.
    Original continuous certificates independently certify every accepted path.
    """
    from .convex_triangle import triangle_closest_many
    events=[]
    def root(gap,upper,pad,identity):
        g0,g1=gap(0.),gap(upper)
        if g1>=-pad or g0<=pad:return
        lo=0.;hi=upper;iterations=0
        for iterations in range(64):
            mid=lo+(hi-lo)/2
            if mid==lo or mid==hi:break
            g=gap(mid)
            if g>=0:lo=mid
            else:hi=mid
            if 0<=g<=pad/4:
                lo=mid;break
        value=gap(lo)
        if 0<=value<=pad:
            events.append(dict(time_from_start_s=lo,g_nf_at_root_m=value,
                roundoff_m=pad,root_iterations=iterations+1,**identity))
    if wall is not None:
        for i,s in sorted(shapes.items()):
            lower=policy.lower_handoff_gap(policy.reference_length(s))['h_lower_m']
            velocity=velocities[i];delta=dt*velocity;finish=s.center_m+delta
            ids=wall.swept_candidates(s.center_m,finish,s.radius_m+lower)
            if not len(ids) or not np.any(delta):continue
            wp,_=triangle_closest_many(s.center_m,wall.triangles[ids])
            d=s.center_m-wp;distance=np.linalg.norm(d,axis=1)
            g0=distance-s.radius_m-lower
            pad=roundoff_length(s.center_m,finish,s.radius_m,wall.triangles)
            # Supporting plane of each original convex triangle. No distance band.
            direction=np.einsum('ij,j->i',d,delta)/distance
            possible=(g0+np.minimum(direction,0)<-pad)&(g0>pad)
            for k in np.flatnonzero(possible):
                tri=wall.triangles[ids[k]:ids[k]+1]
                def value(t):
                    center=s.center_m+t*velocity;point,_=triangle_closest_many(center,tri)
                    vector=center-point[0]
                    return float(np.linalg.norm(vector)-s.radius_m-lower),float(vector@velocity)
                end_g,end_derivative=value(dt);upper=dt;minimum_iterations=0
                if end_g>=-pad:
                    if end_derivative<=0:continue
                    lo=0.;hi=dt
                    for minimum_iterations in range(64):
                        mid=lo+(hi-lo)/2
                        if mid==lo or mid==hi:break
                        if value(mid)[1]<0:lo=mid
                        else:hi=mid
                    upper=lo+(hi-lo)/2
                root(lambda t:value(t)[0],upper,pad,dict(kind='WALL',particle_id=i,
                    triangle_id=int(ids[k]),minimum_search_iterations=minimum_iterations,
                    full_candidate_endpoint_g_nf_m=end_g))
    for i,j in all_pairs(shapes):
        a,b=shapes[i],shapes[j]
        lower=policy.lower_handoff_gap(policy.reference_length(a,b))['h_lower_m']
        pad=roundoff_length(a.center_m,b.center_m,a.radius_m,b.radius_m)
        relative=a.center_m-b.center_m;dv=velocities[i]-velocities[j]
        denominator=float(dv@dv)
        upper=float(np.clip(-(relative@dv)/denominator,0.,dt)) if denominator else 0.
        def gap(t):return float(np.linalg.norm(relative+t*dv)-a.radius_m-b.radius_m-lower)
        root(gap,upper,pad,dict(kind='PAIR',particle_id=i,particle_j_id=j))
    return min(events,key=lambda x:(x['time_from_start_s'],x['particle_id'],x.get('triangle_id',-1))) if events else None


def partitioned_wall_handoff_certificate(start,end,wall,lower):
    """Compose ORIGINAL certificates on the SAME held-velocity straight path.

    Near a finite edge the endpoints' support planes can fail to prove a safe
    glancing path. The minimum-distance witness partitions the proof, not the
    dynamics. No state is accepted unless every original subcertificate passes;
    an actual crossing still fails. No altered normals, gaps, pads or floors.
    """
    from .convex_triangle import triangle_closest_many
    from .physical_time_refinement import MAX_REFINEMENT_DEPTH
    intervals=[(0.,1.)];proved=[];delta=end.center_m-start.center_m
    for _ in range(MAX_REFINEMENT_DEPTH):
        if not intervals:
            return True,dict(proof='UNION_OF_ORIGINAL_CONTINUOUS_HANDOFF_CERTIFICATES',
                proof_partition_only=True,held_velocity_path_unchanged=True,subcertificates=proved)
        lo,hi=intervals.pop(0)
        a=start.moved(start.center_m+lo*delta) if lo else start
        b=start.moved(start.center_m+hi*delta) if hi!=1 else end
        safe,record=wall_handoff_certificate(a,b,wall,lower)
        if safe:
            proved.append(dict(t0_fraction=lo,t1_fraction=hi,certificate=record));continue
        if record['proof']!='HANDOFF_SWEEP_NOT_CERTIFIED':return False,record
        tri=wall.triangles[record['triangle_id']:record['triangle_id']+1]
        left,right=lo,hi
        def value(t):
            point=start.center_m+t*delta;wp,_=triangle_closest_many(point,tri)
            d=point-wp[0]
            return float(np.linalg.norm(d)-start.radius_m-lower),float(d@delta)
        if value(lo)[1]>=0 or value(hi)[1]<=0:return False,record
        for iteration in range(64):
            mid=left+(right-left)/2
            if mid==left or mid==right:break
            if value(mid)[1]<0:left=mid
            else:right=mid
        fraction=left+(right-left)/2
        if not lo<fraction<hi or value(fraction)[0]<-record['roundoff_m']:return False,record
        intervals[0:0]=[(lo,fraction),(fraction,hi)]
    return False,dict(proof='HANDOFF_ORIGINAL_CERTIFICATE_PARTITION_LIMIT',
        certified_subinterval_count=len(proved),unproved_subinterval_count=len(intervals))
