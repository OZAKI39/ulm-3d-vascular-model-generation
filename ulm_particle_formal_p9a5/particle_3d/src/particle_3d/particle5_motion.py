"""P5 sphere dynamics on immutable P3/P4 geometry and physical-time certificates."""
import numpy as np
from .hydrodynamic_resistance import PhysicalNearField, VALIDATION_NEARFIELD_RATIO_MAX, evaluate_block, require_sphere
from .resistance_assembly import assemble_resistance_system
from .resistance_solver import solve_resistance, ResistanceSystemIllConditioned
from .pair_broadphase import all_pairs, broadphase
from .pair_geometry import pair_gap
from .wall_gap import wall_gap, touching_contacts
from .kinematic_contact import ContactConstraint, MultiContactInfeasible
from .particle4_motion import initial_world, ParticleWorldState, pair_interval_certificate
from .physical_time_refinement import refine_interval, TrialNeedsSubdivision, swept_clearance_certificate


def pair_spec(shapes, i, j, *, role='VALIDATION_ONLY'):
    gap = pair_gap(shapes[i], shapes[j], i, j)
    return PhysicalNearField(i, j, gap.gap_m, gap.normal_j_to_i, gap.roundoff_budget_m, role)


def validation_candidates(shapes):
    # P4 AABB reuse. Query dilation covers h<=.01*Reff since Reff<=min(ai,aj).
    # This changes no particle geometry; it is not a production neighbor skin.
    return broadphase(shapes, surface_motion_bounds={i: VALIDATION_NEARFIELD_RATIO_MAX*s.radius_m for i,s in shapes.items()})


def assemble_scene(shapes, free, mu, wall=None, *, use_broadphase=True, role='VALIDATION_ONLY'):
    for shape in shapes.values():
        require_sphere(shape)
    pairs = validation_candidates(shapes) if use_broadphase and role == 'VALIDATION_ONLY' else all_pairs(shapes)
    walls = []
    if wall is not None:
        for i, s in shapes.items():
            g = wall_gap(s, wall)
            walls.append(PhysicalNearField(i, None, g.gap_m, g.normal_inward, g.roundoff_m, role))
    return assemble_resistance_system(shapes, free, walls, [pair_spec(shapes,i,j,role=role) for i,j in pairs], mu=mu)


def resistance_trial(velocity_provider, mu, *, wall=None, boundary_classifier=None, role='VALIDATION_ONLY', on_accept=None):
    def trial(old, dt):
        free = {i: np.concatenate(velocity_provider(i,s,old.time_s)) for i,s in old.shapes.items()}
        system = assemble_scene(old.shapes, free, mu, wall, role=role)
        contacts = []
        pair_records = []
        for i, j in all_pairs(old.shapes):
            gap = pair_gap(old.shapes[i], old.shapes[j], i, j)
            if gap.state == 'PENETRATING':
                raise ValueError('Accepted state contains overlap')
            if gap.state == 'TOUCHING':
                contacts.append(ContactConstraint.pair(gap))
            pair_records.append(evaluate_block(pair_spec(old.shapes,i,j,role=role),old.shapes,mu)[2])
        if wall is not None:
            for i,s in old.shapes.items():
                contacts.extend(ContactConstraint.wall(i,g) for g in touching_contacts(s,wall))
        try:
            solved = solve_resistance(system, particles=old.shapes, constraints=contacts)
        except (ResistanceSystemIllConditioned, MultiContactInfeasible) as error:
            raise TrialNeedsSubdivision(str(error)) from error
        v = {i: solved.velocity[6*k:6*k+3] for k,i in enumerate(system.ids)}
        w = {i: solved.velocity[6*k+3:6*k+6] for k,i in enumerate(system.ids)}
        end = {i:s.moved(s.center_m+dt*v[i]) for i,s in old.shapes.items()}
        events, event = {}, 'ACTIVE'
        if boundary_classifier is not None:
            hits = []
            for i,s in old.shapes.items():
                hit = boundary_classifier.first_event(s.center_m,end[i].center_m)
                if hit is not None:
                    hits.append((hit.segment_fraction,i,hit))
            if hits:
                fraction,i,hit = min(hits,key=lambda x:(x[0],x[1]))
                if not hit.role.startswith('OUTLET_'):
                    raise TrialNeedsSubdivision('Center boundary event '+hit.role)
                dt *= fraction; events[i] = hit.role; event = hit.role
                end = {i:s.moved(s.center_m+dt*v[i]) for i,s in old.shapes.items()}
        proofs = []
        for i,j in all_pairs(end):
            clear, proof = pair_interval_certificate(old.shapes[i],old.shapes[j],end[i],end[j],w[i],w[j],dt)
            if not clear:
                raise TrialNeedsSubdivision(proof)
            proofs.append(proof)
        gaps = [pair_gap(end[i],end[j],i,j) for i,j in all_pairs(end)]
        wall_gaps = {}
        if wall is not None:
            for i,s in end.items():
                gap = wall_gap(s,wall); wall_gaps[i] = gap
                if gap.state == 'PENETRATING':
                    raise TrialNeedsSubdivision('WALL_PENETRATING_TRIAL',gap.gap_m,gap.wall_triangle_id)
                clear,proof,triangle = swept_clearance_certificate(old.shapes[i],s,wall,rotation_angle=0.)
                if not clear:
                    raise TrialNeedsSubdivision(proof,gap.gap_m,triangle)
                proofs.append(proof)
        audit = dict(solved.record, potential_blocks=system.blocks, pair_eligibility=pair_records,
            free_generalized_velocities={i:free[i] for i in system.ids},
            hydro_generalized_velocity=solved.unconstrained,
            evaluated_at_time_s=old.time_s, eligibility_role=role,
            query_bound_role='VALIDATION_QUERY_BOUND_NOT_PHYSICAL_MOTION')
        new = ParticleWorldState(old.time_s+dt,end,v,w,gaps,wall_gaps,audit,events,event)
        if on_accept is not None:
            on_accept(new)
        return new, ';'.join(sorted(set(proofs))) or 'NO_GEOMETRIC_OBSTACLES'
    return trial


def simulate_resistance(shapes, velocity_provider, mu, dt, horizon, *, wall=None, boundary_classifier=None, role='VALIDATION_ONLY', max_depth=48):
    if not np.isfinite([dt,horizon]).all() or min(dt,horizon) <= 0:
        raise ValueError('Positive finite validation dt and horizon required')
    old = initial_world(shapes,velocity_provider,wall)
    rows, ledger = [old.to_dict()], []
    trial = resistance_trial(velocity_provider,mu,wall=wall,boundary_classifier=boundary_classifier,
                             role=role,on_accept=lambda s:rows.append(s.to_dict()))
    for step in range(1,int(np.ceil(horizon/dt))+1):
        old,_ = refine_interval(old,min(step*dt,horizon),trial,ledger=ledger,max_depth=max_depth)
        if old.boundary_event.startswith('OUTLET_'):
            break
    return old,rows,ledger
