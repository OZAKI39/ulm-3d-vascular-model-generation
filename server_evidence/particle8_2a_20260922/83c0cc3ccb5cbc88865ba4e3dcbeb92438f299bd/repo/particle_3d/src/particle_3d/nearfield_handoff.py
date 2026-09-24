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
