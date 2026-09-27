"""Fixed P6 validation fixtures, no new physics, sampling or production choices."""
from dataclasses import replace
from itertools import combinations_with_replacement
import numpy as np
from .particle_shapes import Sphere,Ellipsoid,Capsule,EPS
from .particle4_cases import shapes_for_geometry,place_pair
from .particle3_cases import plane_triangle
from .rbc_orientation import quaternion_from_short_axis,angular_velocity,short_axis
from .wall_geometry import WallGeometry
from .lammps_state import BridgeParticle,PROPERTIES
from .lammps_neighbors import ValidationNeighborPolicy,standalone_candidates
from .particle5_motion import assemble_scene
from .resistance_solver import solve_resistance,contact_jacobian
from .pair_broadphase import all_pairs
from .pair_geometry import pair_gap
from .wall_gap import touching_contacts
from .kinematic_contact import ContactConstraint

BRANCH='dev/particle-6-lammps-bridge-20260920'


def mixed_particles():
    sphere,oblate,capsule=shapes_for_geometry()
    q=quaternion_from_short_axis([1,2,3]);retained=quaternion_from_short_axis([2,-1,3])
    shapes={11:Sphere([0,0,0],1e-6),101:Sphere([2.001e-6,0,0],1e-6),
        307:oblate.moved([10e-6,0,0]),409:capsule.moved([20e-6,0,0]),
        503:oblate.moved([30e-6,1e-6,0]),887:capsule.moved([40e-6,0,0])}
    records=[]
    for i,s in shapes.items():
        v,w=mixed_provider(i,s,0)
        records.append(BridgeParticle.from_shape(i,s,q=q if isinstance(s,Ellipsoid) else retained,
            velocity=v,omega=w))
    return records,ValidationNeighborPolicy(15e-6,.5e-6,'FIXED_MIXED_P4_METADATA_AND_NEIGHBOR_VALIDATION')


def mixed_provider(i,shape,time):
    if isinstance(shape,Sphere):return np.array([2e-7 if i==11 else -2e-7,0,0]),np.array([0.,0.,.2])
    v=np.array([1e-7,.5e-7,0])
    if isinstance(shape,Capsule):return v,np.zeros(3)
    g=np.array([[0.,.4,0],[0,0,.1],[0,0,0]])
    a,c=shape.axes_m[0],shape.axes_m[2]
    return v,angular_velocity(shape.rotation[:,2],(c*c-a*a)/(c*c+a*a),g)


def sphere_particles():
    shapes={17:Sphere([0,0,1e-6],1e-6),203:Sphere([2e-6,0,1e-6],1e-6),901:Sphere([4.001e-6,0,1.005e-6],1e-6)}
    records=[BridgeParticle.from_shape(i,s,velocity=sphere_provider(i,s,0)[0],omega=sphere_provider(i,s,0)[1]) for i,s in shapes.items()]
    return records,ValidationNeighborPolicy(8e-6,.25e-6,'P5_SPHERES_WITH_SIMULTANEOUS_WALL_CONTACT_AND_NEAR_PAIR'),WallGeometry([plane_triangle()])


def sphere_provider(i,shape,time):
    return np.array([{17:1e-6,203:0.,901:-1e-6}[i],0.,-.5e-6]),np.array([.1,.2,.3])


def six_pair_shapes():
    base=shapes_for_geometry();shapes={};rows=[]
    for k,(i,j) in enumerate(combinations_with_replacement(range(3),2)):
        a,b=place_pair(base[i],base[j],[1,0,0],0.)
        offset=np.array([40e-6*k,0,0]);ida=101+107*k;idb=ida+41
        shapes[ida]=a.moved(a.center_m+offset);shapes[idb]=b.moved(b.center_m+offset)
        rows.append(dict(ids=[ida,idb],shape_pair=[a.mode,b.mode]))
    shapes[991]=Sphere([0,12e-6,0],1e-6) # query candidates outside exact contact / near-field eligibility
    return shapes,rows,ValidationNeighborPolicy(18e-6,1e-6,'ALL_SIX_SHAPE_COMBINATIONS_WITH_FIXED_CONTACT_WITNESSES')


def state_errors(a,b,*,steps=1):
    a={p.particle_id:p for p in a};b={p.particle_id:p for p in b}
    if set(a)!=set(b):raise ValueError('Stable IDs differ')
    errors=dict(position_m=0.,quaternion=0.,velocity_m_s=0.,omega_s_inv=0.,geometry=0.,orientation_rad=0.)
    exact=True;field_errors={}
    for i in sorted(a):
        pa,pb=a[i],b[i]
        for field in ['type_code','mode_code','position']+[f for _,f,_ in PROPERTIES]:
            x=np.asarray(getattr(pa,field));y=np.asarray(getattr(pb,field))
            delta=float(np.max(np.abs(x-y)));field_errors[field]=max(field_errors.get(field,0.),delta)
            exact=exact and np.array_equal(x,y)
        for field,key in [('position','position_m'),('q','quaternion'),('velocity','velocity_m_s'),('omega','omega_s_inv')]:
            errors[key]=max(errors[key],float(np.max(np.abs(getattr(pa,field)-getattr(pb,field)))))
        for _,field,_ in PROPERTIES:
            if field not in ['type_code','mode_code','q','velocity','omega']:
                errors['geometry']=max(errors['geometry'],float(np.max(np.abs(np.asarray(getattr(pa,field))-getattr(pb,field)))))
        xa,xb=short_axis(pa.q),short_axis(pb.q)
        errors['orientation_rad']=max(errors['orientation_rad'],float(np.arctan2(np.linalg.norm(np.cross(xa,xb)),abs(xa@xb))))
    scales={key:max(float(np.max(np.abs(getattr(p,field)))) for p in list(a.values())+list(b.values())) for field,key in [('position','position_m'),('q','quaternion'),('velocity','velocity_m_s'),('omega','omega_s_inv')]}
    budgets={key:512*EPS*max(1,steps)*max(scale,np.finfo(float).tiny) for key,scale in scales.items()}
    budgets['orientation_rad']=512*EPS*max(1,steps)
    return dict(**errors,exact_equal=bool(exact),field_errors=field_errors,roundoff_budgets=budgets)


def resistance_snapshot(records,provider,mu,wall,*,stepper=None):
    shapes={p.particle_id:p.shape() for p in records};free={i:np.r_[provider(i,s,0)[0],provider(i,s,0)[1]] for i,s in shapes.items()}
    system=assemble_scene(shapes,free,mu,wall) if stepper is None else stepper._assemble(shapes,free,mu,wall)
    contacts=[]
    for i,j in all_pairs(shapes):
        g=pair_gap(shapes[i],shapes[j],i,j)
        if g.state=='TOUCHING':contacts.append(ContactConstraint.pair(g))
    if wall is not None:
        for i,s in shapes.items():contacts.extend(ContactConstraint.wall(i,g) for g in touching_contacts(s,wall))
    solved=solve_resistance(system,particles=shapes,constraints=contacts);j,ordered=contact_jacobian(system,shapes,contacts)
    return dict(R=system.matrix.toarray(),b=system.rhs,U=solved.velocity,J=j,
        contact_ids=[c.canonical_id for c in ordered],active_constraints=[c.canonical_id for c,lam in zip(ordered,solved.record.get('multipliers',[])) if lam>0],
        solver=solved.record,blocks=system.blocks)
