"""Validation-only fixtures and independent scalar references for Particle-5."""
from pathlib import Path
import numpy as np
from .particle_shapes import Sphere, unit
from .hydrodynamic_resistance import PhysicalNearField
from .resistance_assembly import assemble_resistance_system
from .resistance_solver import solve_resistance, contact_jacobian
from .particle5_motion import assemble_scene, simulate_resistance
from .particle4_motion import simulate_world
from .particle3_cases import plane_triangle
from .wall_geometry import WallGeometry
from .pair_geometry import pair_gap
from .kinematic_contact import ContactConstraint, project_contacts

P5_BRANCH = 'dev/particle-5-resistance-lubrication-20260920'
P4_COMMIT = 'a1e09cb3e2db7c91d7adc75adf40a0f939a880bf'
RATIOS = [.1, .01, .001, .0001]
ROLE = 'VALIDATION_GAPS_ONLY'


def scalar_pair_solution(ai, aj, mu, h, ui, uj):
    """Independent two scalar balance elimination; no global assembly/solve."""
    zi, zj = 6*np.pi*mu*ai, 6*np.pi*mu*aj
    reff = ai*aj/(ai+aj); zp = 6*np.pi*mu*reff**2/h
    relative = (ui-uj)/(1+zp/zi+zp/zj)
    return np.array([ui-zp/zi*relative, uj+zp/zj*relative])


def benchmarks(mu):
    stokes, walls, pairs = [], [], []
    for viscosity in [mu, .001, .01]:
        for radius in [.25e-6, 1e-6, 2e-6, 5e-6]:
            shape = Sphere(np.zeros(3),radius)
            free = np.array([.001,-.002,.003,2.,-3.,4.])
            s = assemble_resistance_system({1:shape},{1:free},mu=viscosity); v=solve_resistance(s)
            stokes.append(dict(radius_m=radius,mu_pa_s=viscosity,translation=s.self_diagonal[0],rotation=s.self_diagonal[3],
                analytic_translation=6*np.pi*viscosity*radius,analytic_rotation=8*np.pi*viscosity*radius**3,
                isolated_velocity_error=float(np.max(np.abs(v.velocity[:3]-free[:3]))),
                isolated_angular_error=float(np.max(np.abs(v.velocity[3:]-free[3:]))),solver=v.record))
    for ratio in RATIOS:
        a=1e-6; h=a*ratio; shape=Sphere([0,0,a+h],a)
        free=np.array([2e-6,3e-6,-1e-6,1.,2.,3.])
        block=PhysicalNearField(1,None,h,[0,0,1],0.,ROLE)
        s=assemble_resistance_system({1:shape},{1:free},[block],mu=mu);v=solve_resistance(s)
        analytic=free.copy();analytic[2]*=h/(a+h)
        walls.append(dict(gap_m=h,gap_ratio=ratio,relevant_radius_m=a,coefficient=s.blocks[0]['coefficient_kg_s'],
            normalized_scaling=s.blocks[0]['coefficient_kg_s']*h/(6*np.pi*mu*a*a),normal_attenuation=v.velocity[2]/free[2],
            analytic_attenuation=h/(a+h),velocity=v.velocity,free=free,analytic=analytic,
            analytic_error=float(np.max(np.abs(v.velocity[:3]-analytic[:3]))),tangent_error=float(np.linalg.norm(v.velocity[:2]-free[:2])),
            angular_error=float(np.linalg.norm(v.velocity[3:]-free[3:])),eligibility=s.blocks[0],solver=v.record))
        for ai,aj in [(1e-6,1e-6),(.5e-6,2e-6)]:
            reff=ai*aj/(ai+aj);h=ratio*reff
            shapes={7:Sphere([0,0,0],ai),19:Sphere([ai+aj+h,0,0],aj)}
            free={7:np.array([1e-6,2e-6,0,1,2,3]),19:np.array([-2e-6,2e-6,0,4,5,6])}
            spec=PhysicalNearField(7,19,h,[-1,0,0],0.,ROLE)
            s=assemble_resistance_system(shapes,free,pair_blocks=[spec],mu=mu);v=solve_resistance(s)
            analytic=scalar_pair_solution(ai,aj,mu,h,free[7][0],free[19][0])
            delta=v.velocity[0]-v.velocity[6]
            pairs.append(dict(radius_i_m=ai,radius_j_m=aj,gap_m=h,gap_ratio=ratio,relevant_radius_m=reff,
                coefficient=s.blocks[0]['coefficient_kg_s'],normalized_scaling=s.blocks[0]['coefficient_kg_s']*h/(6*np.pi*mu*reff**2),
                relative_attenuation=delta/(free[7][0]-free[19][0]),analytic=analytic,
                velocities=v.velocity[[0,6]],analytic_error=float(np.max(np.abs(v.velocity[[0,6]]-analytic))),
                action_reaction_sum=(s.blocks[0]['coefficient_kg_s']*delta*np.array([1.,-1.])).sum(),
                eligibility=s.blocks[0],solver=v.record))
    return stokes,walls,pairs


def algebraic_pair_block(normal=(1,2,3)):
    """ζtest=1 in internally normalized TEST units; never a physical block."""
    n=unit(normal);row=np.r_[n,np.zeros(3),-n,np.zeros(3)]
    return np.outer(row,row)


def sparse_scene(mu, count=12, *, candidates=True):
    shapes={};free={};walls=[]
    rng=np.random.default_rng(2026092505)
    for k in range(count):
        a=(1+.1*(k%3))*1e-6
        if k%2:
            previous=shapes[k-1];center=previous.center_m+np.array([previous.radius_m+a+1e-9,0,0])
        else:
            center=np.array([12e-6*(k//2),0,3e-6])
        shapes[k]=Sphere(center,a);free[k]=rng.normal(size=6)*np.array([1e-5]*3+[10.]*3)
        if k%3==0:
            walls.append(PhysicalNearField(k,None,.005*a,[0,0,1],0.))
    from .particle5_motion import validation_candidates,pair_spec
    from .pair_broadphase import all_pairs
    candidate_pairs=validation_candidates(shapes) if candidates else all_pairs(shapes)
    system=assemble_resistance_system(shapes,free,walls,[pair_spec(shapes,i,j) for i,j in candidate_pairs],mu=mu)
    return shapes,system,candidate_pairs


def contact_cases(mu):
    rows=[]
    # Two touching unequal spheres, each at a distinct small positive plane gap.
    ai,aj=1e-6,2e-6;zi=ai*1.001;zj=aj*1.005
    dx=np.sqrt((ai+aj)**2-(zj-zi)**2)
    two={1:Sphere([0,0,zi],ai),2:Sphere([dx,0,zj],aj)}
    three={1:Sphere([0,0,10e-6],.5e-6),2:Sphere([1.5e-6,0,10e-6],1e-6),3:Sphere([4e-6,0,10e-6],1.5e-6)}
    for name,shapes in [('two_near_wall',two),('three_spheres',three)]:
        free={i:np.array([(-1)**(i+1)*i*1e-6,0,-i*1e-6,0,0,0]) for i in shapes}
        walls=[]
        if name=='two_near_wall':
            walls=[PhysicalNearField(i,None,s.center_m[2]-s.radius_m,[0,0,1],0.) for i,s in shapes.items()]
        system=assemble_resistance_system(shapes,free,walls,mu=mu)
        gaps=[pair_gap(shapes[i],shapes[j],i,j) for i in shapes for j in shapes if i<j]
        contacts=[ContactConstraint.pair(g) for g in gaps if g.state=='TOUCHING']
        p5=solve_resistance(system,particles=shapes,constraints=contacts)
        v={i:p5.unconstrained[6*k:6*k+3] for k,i in enumerate(system.ids)}
        w={i:p5.unconstrained[6*k+3:6*k+6] for k,i in enumerate(system.ids)}
        p4=project_contacts(shapes,v,w,contacts)
        u4=np.concatenate([np.r_[p4.velocities[i],p4.omegas[i]] for i in system.ids])
        j,_=contact_jacobian(system,shapes,contacts);delta=u4-p5.unconstrained
        rows.append(dict(case=name,comparison_reference='SAME_UNCONSTRAINED_HYDRO_VELOCITY',ids=system.ids,
            centers_m=[shapes[i].center_m for i in system.ids],radii_m=[shapes[i].radius_m for i in system.ids],
            hydro=p5.unconstrained,p4=u4,p5=p5.velocity,p4_normal_speeds=j@u4,p5_normal_speeds=j@p5.velocity,
            p4_resistance_objective=float(.5*delta@(system.matrix@delta)),p5_resistance_objective=p5.record['resistance_correction_objective'],
            velocity_difference=float(np.linalg.norm(u4-p5.velocity)),p4_audit=p4.record,p5_audit=p5.record))
    return rows


def synthetic_motion(mu, kind, divisor=1, *, p4=False):
    a=1e-6
    if kind=='wall':
        wall=WallGeometry([plane_triangle()]);shapes={1:Sphere([0,0,1.1*a],a)}
        provider=lambda i,s,t:(np.array([0,0,-a]),np.zeros(3))
        dt,horizon=.1/divisor,4.
    else:
        wall=None;shapes={1:Sphere([0,0,0],a),2:Sphere([2.05*a,0,0],a)}
        provider=lambda i,s,t:(np.array([a if i==1 else -a,0,0]),np.zeros(3))
        dt,horizon=.025/divisor,1.
    if p4:
        end,states,ledger=simulate_world(shapes,provider,dt,horizon,wall=wall)
    else:
        end,states,ledger=simulate_resistance(shapes,provider,mu,dt,horizon,wall=wall,role=ROLE)
    summary=dict(kind=kind,model='P4' if p4 else 'P5',dt_s=dt,horizon_s=horizon,
        final_time_s=end.time_s,time_coverage_error_s=abs(sum(r['dt_s'] for r in ledger)-horizon),
        role='VALIDATION_ONLY',contact_events=sum(s['projection'].get('contact_count',0)>0 for s in states),
        minimum_gap_m=min(min(p['wall_gap_m'] for p in s['particles']) if kind=='wall' else s['pair_gaps'][0]['gap_m'] for s in states),
        final_centers_m=[end.shapes[i].center_m for i in sorted(end.shapes)])
    return summary,states,ledger
