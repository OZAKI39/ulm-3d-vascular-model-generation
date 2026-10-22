"""Many-particle true-time stepping with simultaneous WALL / pair constraints."""
from dataclasses import dataclass, replace
import numpy as np
from .particle_shapes import Sphere, Ellipsoid, Capsule, EPS, roundoff_length
from .rbc_orientation import rotation_matrix, advance_orientation
from .pair_geometry import pair_gap
from .pair_broadphase import all_pairs, broadphase
from .kinematic_contact import ContactConstraint, project_contacts
from .physical_time_refinement import TrialNeedsSubdivision, refine_interval, swept_clearance_certificate
from .wall_gap import wall_gap, touching_contacts


class InitialParticleOverlap(ValueError):
    def __init__(self, gap):
        self.record = dict(status='INITIAL_PARTICLE_OVERLAP', pair=gap.to_dict())
        super().__init__(str(self.record))


def move_shape(shape, velocity, omega, dt):
    center = shape.center_m+dt*velocity
    if isinstance(shape, Ellipsoid):
        delta = rotation_matrix(advance_orientation([1., 0, 0, 0], omega, dt))
        return Ellipsoid(center, shape.axes_m, delta@shape.rotation)
    moved=shape.moved(center)
    if isinstance(shape,Capsule):
        # P3's constructor normalizes any new axis. Pure P4 translation must
        # retain the exact existing axis, rather than normalize it a second time.
        object.__setattr__(moved,'axis_world',shape.axis_world)
    return moved


def isolated_sphere_contact_path(shapes, free, dt):
    """Exact integration of this same frictionless velocity projection.

For two touching spheres with constant free relative velocity w, n'= (w -
(w.n)n)/R while closing. Thus c=n.w_hat obeys c'=|w|(1-c²)/R. Integrating
this ODE avoids a sequence of vanishing free-flight/chord events on a curved
contact surface. Centers are the integrated solution, never projected apart.
The actual initial separation R is retained, including its original roundoff.
"""
    i,j=sorted(shapes);a,b=shapes[i],shapes[j]
    relative=a.center_m-b.center_m;distance=np.linalg.norm(relative);n=relative/distance
    w=free[i]-free[j];speed=np.linalg.norm(w);mean=(free[i]+free[j])/2
    if speed==0 or w@n>=0:return None
    direction=w/speed;c=float(np.clip(n@direction,-1.,1.));perpendicular=n-c*direction
    transverse=np.linalg.norm(perpendicular)
    if transverse<=64*EPS:
        new_relative=relative.copy();new_relative_velocity=np.zeros(3);duration=dt
    else:
        e=perpendicular/transverse
        parameter=np.arcsinh(c/transverse);release=-distance*parameter/speed;duration=min(dt,release)
        updated_c=np.tanh(parameter+speed*duration/distance)
        updated_n=updated_c*direction+np.sqrt(max(0.,1-updated_c**2))*e
        new_relative=distance*updated_n+(dt-duration)*w
        new_relative_velocity=w-min(float(w@updated_n),0.)*updated_n if duration==dt else w
    center=(a.center_m+b.center_m)/2+mean*dt
    end={i:a.moved(center+new_relative/2),j:b.moved(center-new_relative/2)}
    velocities={i:mean+new_relative_velocity/2,j:mean-new_relative_velocity/2}
    return end,velocities,dict(method='EXACT_ISOLATED_SPHERE_PROJECTED_VELOCITY_ODE',contact_duration_s=duration,
        preserved_initial_center_separation_m=distance,position_projection=False,
        interval_average_velocities_m_s={k:(end[k].center_m-shapes[k].center_m)/dt for k in [i,j]})


def _plane_gap(a, b, normal):
    return float(normal@(a.support(-normal)-b.support(normal)))


def pair_interval_certificate(a, b, end_a, end_b, omega_a, omega_b, dt):
    """Continuous supporting-plane certificate, including rotating ellipsoids.

For an ellipsoid h(n,t)=||diag(axes) R(t)^T n||, |h''| is bounded by
(amax^2/amin+amax)*|Omega|^2. Endpoint chord minus B*dt^2/8 is thus a
rigorous support-gap lower bound. Recursive *geometric* subinterval checks
certify the same trial; an uncertified trial requests real physical splitting.
"""
    if type(a) is not type(end_a) or type(b) is not type(end_b):
        return False, 'MODE_CHANGE_REQUIRES_P3_RECOMPUTED_PAIR_PATH'
    for old, new in [(a, end_a), (b, end_b)]:
        if isinstance(old, Capsule) and (old.radius_m != new.radius_m or old.cylindrical_length_m != new.cylindrical_length_m or not np.array_equal(old.axis_world, new.axis_world)):
            return False, 'P3_CAPSULE_UPDATE_REQUIRES_RECOMPUTED_PAIR_PATH'
    g0 = pair_gap(a, b); g1 = pair_gap(end_a, end_b)
    pad = max(g0.roundoff_budget_m, g1.roundoff_budget_m)
    if min(g0.gap_m, g1.gap_m) < -pad:
        return False, 'PAIR_PENETRATING_ENDPOINT'
    angles = [np.linalg.norm(w)*dt if isinstance(s, Ellipsoid) else 0. for s, w in [(a, omega_a), (b, omega_b)]]
    movement = sum(np.linalg.norm(x.center_m-y.center_m)+x.bounding_radius_m*angle for x, y, angle in zip([a,b],[end_a,end_b],angles))
    if min(g0.gap_m, g1.gap_m) > movement+pad:
        return True, 'PAIR_SURFACE_MOTION_BOUND'
    second = sum((max(s.axes_m)**2/min(s.axes_m)+max(s.axes_m))*np.linalg.norm(w)**2
                 for s, w in [(a,omega_a),(b,omega_b)] if isinstance(s, Ellipsoid))
    velocities = [(end_a.center_m-a.center_m)/dt, (end_b.center_m-b.center_m)/dt]
    for normal in [g0.normal_j_to_i, g1.normal_j_to_i]:
        stack = [(0., dt, a, b, end_a, end_b, 0)]; certified = True
        while stack:
            t0,t1,x0,y0,x1,y1,depth = stack.pop()
            lo,hi = _plane_gap(x0,y0,normal),_plane_gap(x1,y1,normal)
            if min(lo,hi)-second*(t1-t0)**2/8 >= -pad:
                continue
            if min(lo,hi) < -pad or depth >= 12:
                certified = False; break
            middle=(t0+t1)/2
            xm=move_shape(a,velocities[0],omega_a,middle);ym=move_shape(b,velocities[1],omega_b,middle)
            stack.extend([(t0,middle,x0,y0,xm,ym,depth+1),(middle,t1,xm,ym,x1,y1,depth+1)])
        if certified:
            return True, 'CONTINUOUS_PAIR_SUPPORT_PLANE_SECOND_DERIVATIVE_BOUND'
    return False, 'PAIR_SWEEP_NOT_CERTIFIED'


@dataclass(frozen=True)
class ParticleWorldState:
    time_s: float
    shapes: dict
    velocities: dict
    omegas: dict
    pair_gaps: list
    wall_gaps: dict
    projection: dict
    boundary_events: dict
    boundary_event: str = 'ACTIVE'

    def to_dict(self):
        particles=[]
        for i in sorted(self.shapes):
            s=self.shapes[i]
            row=dict(particle_id=i,shape_mode=s.mode,center_m=s.center_m,velocity_m_s=self.velocities[i],omega_s_inv=self.omegas[i],boundary_event=self.boundary_events.get(i,'ACTIVE'))
            if isinstance(s,Sphere):row['radius_m']=s.radius_m
            elif isinstance(s,Ellipsoid):row.update(axes_m=s.axes_m,rotation=s.rotation,short_axis=s.rotation[:,2])
            else:row.update(radius_m=s.radius_m,cylindrical_length_m=s.cylindrical_length_m,axis_world=s.axis_world)
            if i in self.wall_gaps:row.update(wall_gap_m=self.wall_gaps[i].gap_m,wall_roundoff_m=self.wall_gaps[i].roundoff_m)
            particles.append(row)
        return dict(time_s=self.time_s,particles=particles,pair_gaps=[g.to_dict() for g in self.pair_gaps],projection=self.projection,boundary_event=self.boundary_event)


def initial_world(shapes, velocity_provider, wall=None):
    shapes=dict(sorted(shapes.items()));gaps=[]
    for i,j in all_pairs(shapes):
        gap=pair_gap(shapes[i],shapes[j],i,j)
        if gap.state=='PENETRATING':raise InitialParticleOverlap(gap)
        gaps.append(gap)
    wall_gaps={} if wall is None else {i:wall_gap(s,wall) for i,s in shapes.items()}
    if any(g.state=='PENETRATING' for g in wall_gaps.values()):raise ValueError('INITIAL_WALL_OVERLAP')
    velocities={};omegas={}
    for i,s in shapes.items():velocities[i],omegas[i]=map(np.asarray,velocity_provider(i,s,0.))
    return ParticleWorldState(0.,shapes,velocities,omegas,gaps,wall_gaps,{}, {})


def world_trial(velocity_provider, *, wall=None, boundary_classifier=None, on_accept=None, use_broadphase=True):
    def trial(old,dt):
        free={};omega={}
        for i,s in old.shapes.items():free[i],omega[i]=map(np.asarray,velocity_provider(i,s,old.time_s))
        pairs=broadphase(old.shapes) if use_broadphase else all_pairs(old.shapes)
        contacts=[]
        for i,j in pairs:
            gap=pair_gap(old.shapes[i],old.shapes[j],i,j)
            if gap.state=='PENETRATING':raise ValueError('Accepted state contains particle overlap')
            if gap.state=='TOUCHING':contacts.append(ContactConstraint.pair(gap))
        if wall is not None:
            for i,s in old.shapes.items():
                for gap in touching_contacts(s,wall):contacts.append(ContactConstraint.wall(i,gap))
        projected=project_contacts(old.shapes,free,omega,contacts)
        exact=None
        if wall is None and boundary_classifier is None and len(old.shapes)==2 and all(isinstance(s,Sphere) for s in old.shapes.values()) and contacts:
            exact=isolated_sphere_contact_path(old.shapes,free,dt)
        if exact is not None:
            end,velocities,integration=exact
            gaps=[pair_gap(end[i],end[j],i,j) for i,j in all_pairs(end)]
            if any(g.state=='PENETRATING' for g in gaps):raise TrialNeedsSubdivision('EXACT_SPHERE_ODE_NUMERICAL_WITNESS_FAILURE')
            record=dict(projected.record,integration=integration)
            new=ParticleWorldState(old.time_s+dt,end,velocities,projected.omegas,gaps,{},record,{},'ACTIVE')
            if on_accept is not None:on_accept(new)
            return new,'ANALYTIC_CONTINUOUS_SPHERE_CONTACT_ODE; FIXED_RADIUS; PAIRWISE_OPPOSITE_TRANSLATION'
        end={i:move_shape(s,projected.velocities[i],projected.omegas[i],dt) for i,s in old.shapes.items()}
        events={};event='ACTIVE'
        if boundary_classifier is not None:
            hits=[]
            for i,s in old.shapes.items():
                hit=boundary_classifier.first_event(s.center_m,end[i].center_m)
                if hit is not None:hits.append((hit.segment_fraction,i,hit))
            if hits:
                fraction,i,hit=min(hits,key=lambda v:(v[0],v[1]))
                if not hit.role.startswith('OUTLET_'):raise TrialNeedsSubdivision('Center boundary event '+hit.role)
                dt*=fraction;events[i]=hit.role;event=hit.role
                end={i:move_shape(s,projected.velocities[i],projected.omegas[i],dt) for i,s in old.shapes.items()}
        bounds={i:np.linalg.norm(end[i].center_m-s.center_m)+s.bounding_radius_m*np.linalg.norm(projected.omegas[i])*dt if isinstance(s,Ellipsoid) else np.linalg.norm(end[i].center_m-s.center_m) for i,s in old.shapes.items()}
        swept=broadphase(old.shapes,end_shapes=end,surface_motion_bounds=bounds) if use_broadphase else all_pairs(old.shapes)
        proofs=[]
        for i,j in swept:
            clear,proof=pair_interval_certificate(old.shapes[i],old.shapes[j],end[i],end[j],projected.omegas[i],projected.omegas[j],dt)
            if not clear:raise TrialNeedsSubdivision(proof)
            proofs.append(proof)
        gaps=[pair_gap(end[i],end[j],i,j) for i,j in all_pairs(end)]
        if any(g.state=='PENETRATING' for g in gaps):raise TrialNeedsSubdivision('PAIR_PENETRATING_TRIAL')
        wall_gaps={}
        if wall is not None:
            for i,s in end.items():
                gap=wall_gap(s,wall);wall_gaps[i]=gap
                if gap.state=='PENETRATING':raise TrialNeedsSubdivision('WALL_PENETRATING_TRIAL',gap.gap_m,gap.wall_triangle_id)
                clear,proof,triangle=swept_clearance_certificate(old.shapes[i],s,wall,rotation_angle=np.linalg.norm(projected.omegas[i])*dt if isinstance(s,Ellipsoid) else 0.)
                if not clear:raise TrialNeedsSubdivision(proof,gap.gap_m,triangle)
                proofs.append(proof)
        new=ParticleWorldState(old.time_s+dt,end,projected.velocities,projected.omegas,gaps,wall_gaps,projected.record,events,event)
        if on_accept is not None:on_accept(new)
        return new,';'.join(sorted(set(proofs))) or 'DISJOINT_SWEPT_AABB'
    return trial


def simulate_world(shapes,velocity_provider,dt,horizon,*,wall=None,boundary_classifier=None,use_broadphase=True):
    old=initial_world(shapes,velocity_provider,wall);rows=[old.to_dict()];ledger=[]
    trial=world_trial(velocity_provider,wall=wall,boundary_classifier=boundary_classifier,on_accept=lambda s:rows.append(s.to_dict()),use_broadphase=use_broadphase)
    for step in range(1,int(np.ceil(horizon/dt))+1):
        old,_=refine_interval(old,min(step*dt,horizon),trial,ledger=ledger)
        if old.boundary_event.startswith('OUTLET_'):break
    return old,rows,ledger
