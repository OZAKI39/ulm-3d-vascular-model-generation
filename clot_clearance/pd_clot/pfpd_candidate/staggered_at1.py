"""Early AT1 alternating minimization using the existing ADR callback API.

All attempted iterates are saved. Only simultaneously converged load steps
are accepted. High damage stops before installing the proposed phase state.
"""
from pathlib import Path
import json
import numpy as np
from .at1_evolution import at1_constrained_update, critical_driver, audit_update, HighDamageGate, NAME
from .adr_adapter_at1 import fixture, boundary_view, pair_weights, equilibrate_tensile, tensile_domain
from .baqp import evaluate


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')


def energy_parts(family, state, s_pf_bond, yc):
    crack = float(np.dot(pair_weights(family), 2*yc*s_pf_bond))
    return dict(elastic_J=state.energy, crack_J=crack, total_J=state.energy+crack)


def run_path(config, strains, directory):
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=False)
    write_json(out/'CONFIG.json', config)
    cloud, family, material, top = fixture(config)
    x = family.X.copy(); s_pf_bond = np.zeros(len(family.pairs))
    accepted_x, accepted_s = x.copy(), s_pf_bond.copy()
    yc = critical_driver(config['phase']['Gc_demo_J_per_m2'], family.horizon, config['phase']['kernel'])
    tol = config['tolerances']; ec = config['mechanical_equilibrium']
    weights = pair_weights(family); accepted = []; iterations = []
    work = 0.; status = 'CONVERGED'; event = None
    np.savez_compressed(out/'reference.npz', X=family.X, pairs=family.pairs, volume=family.volume,
                        weights=weights, top=top, fixed=cloud.fixed)
    for load_step, strain in enumerate(strains):
        step = out/f'load_{load_step:03d}'; step.mkdir()
        previous_energy = None; converged_checks = 0; load_converged = False
        view = boundary_view(cloud, top, strain)
        x[view.fixed] = view.X[view.fixed]
        try:
            for iteration in range(config['staggered']['maximum_iterations']):
                stem = step/f'iteration_{iteration:03d}'
                mechanical = equilibrate_tensile(cloud, family, material, top, x, s_pf_bond, strain, config)
                write_json(str(stem)+'_mechanical.json', dict(summary=mechanical.summary,
                                                            checks=mechanical.checks, history=mechanical.history))
                x = mechanical.x.copy()
                np.savez_compressed(str(stem)+'_mechanical.npz', x=x, s_pf_bond=s_pf_bond)
                if not mechanical.converged:
                    raise RuntimeError('MECHANICAL_EQUILIBRIUM_FAILED_BEFORE_PHASE_UPDATE')
                before = evaluate(family, x, material, g=(1-s_pf_bond)**2,
                                  h=np.ones(len(family.pairs)))
                min_green = tensile_domain(before, tol['tensile_min_green_strain'])
                # Save Y and old s BEFORE attempting the irreversible update.
                np.savez_compressed(str(stem)+'_driver.npz', Y=before.W0, x=x, s_old=s_pf_bond)
                update = at1_constrained_update(s_pf_bond, before.W0, yc, config['phase']['sc'])
                audit = audit_update(update, before.W0, yc, kkt_tolerance=tol['kkt_scaled'],
                                     energy_tolerance=tol['phase_energy_relative'])
                if not audit['passed']:
                    raise RuntimeError('AT1_KKT_OR_BOND_ENERGY_AUDIT_FAILED')
                after = evaluate(family, x, material, g=(1-update.s_new)**2,
                                 h=np.ones(len(family.pairs)))
                changes = {name:float(np.max(abs(getattr(after,name)-getattr(before,name))))
                           for name in ['K','Fbar','Ftilde','shape_i','shape_j']}
                if any(not np.array_equal(getattr(after,k), getattr(before,k)) for k in changes):
                    raise RuntimeError('KINEMATIC_IDENTITY_VIOLATION')
                energy_before = energy_parts(family, before, s_pf_bond, yc)
                energy_after = energy_parts(family, after, update.s_new, yc)
                delta_energy = energy_after['total_J']-energy_before['total_J']
                budget = tol['phase_energy_absolute_J'] + tol['phase_energy_relative']*abs(energy_before['total_J'])
                if delta_energy > budget:
                    raise RuntimeError('FIXED_X_TOTAL_ENERGY_INCREASE')
                phase_step = float(np.max(abs(update.s_new-s_pf_bond)))
                force_residual = float(np.linalg.norm(after.force[~cloud.fixed])/ec['force_floor_N'])
                energy_change = (None if previous_energy is None else
                    abs(energy_after['total_J']-previous_energy)/max(abs(previous_energy), ec['energy_floor_J']))
                okay = (force_residual <= ec['relative_force_residual_tolerance'] and
                        phase_step <= tol['phase_step'] and energy_change is not None and
                        energy_change <= tol['stagger_energy_relative'])
                row = dict(load_step=load_step, strain=strain, iteration=iteration,
                           mechanical_converged=mechanical.converged,
                           mechanical_iterations=mechanical.summary['iterations'],
                           pre_phase_force_residual=mechanical.summary['final']['relative_force_residual'],
                           relative_force_residual=force_residual, maximum_delta_s=phase_step,
                           relative_stagger_energy_change=energy_change, phase_energy_change_J=delta_energy,
                           energy_before=energy_before, energy_after=energy_after,
                           maximum_s=float(update.s_new.max()), mean_s=float(np.average(update.s_new, weights=weights)),
                           min_green_strain=min_green, kinematic_changes=changes, audit=audit,
                           simultaneous_criteria_passed=bool(okay))
                iterations.append(row); write_json(str(stem)+'_phase.json', row)
                np.savez_compressed(str(stem)+'_phase.npz', x=x, Y=before.W0, weights=weights,
                                    **{k:getattr(update,k) for k in update.__dataclass_fields__})
                s_pf_bond = update.s_new.copy()
                previous_energy = energy_after['total_J']
                converged_checks = converged_checks+1 if okay else 0
                if converged_checks > config['staggered'].get('extra_converged_iterations', 0):
                    load_converged = True
                    displacement = float(strain*np.ptp(family.X[:,2]))
                    reaction = float(-after.force[top,2].sum())
                    if accepted:
                        work += .5*(reaction+accepted[-1]['reaction_N'])*(displacement-accepted[-1]['displacement_m'])
                    accepted_row = dict(load_step=load_step, strain=strain, displacement_m=displacement,
                                        reaction_N=reaction, maximum_s=float(s_pf_bond.max()),
                                        mean_s=float(np.average(s_pf_bond,weights=weights)),
                                        external_work_J=work, **energy_after,
                                        energy_balance_residual_J=energy_after['total_J']-work,
                                        relative_force_residual=force_residual,
                                        maximum_delta_s=phase_step, stagger_iterations=iteration+1,
                                        relative_stagger_energy_change=energy_change,
                                        minimum_phase_increment=float((s_pf_bond-accepted_s).min()))
                    accepted.append(accepted_row)
                    accepted_x, accepted_s = x.copy(), s_pf_bond.copy()
                    write_json(step/'ACCEPTED.json', accepted_row)
                    np.savez_compressed(step/'accepted_state.npz', x=x, s_pf_bond=s_pf_bond, Y=after.W0)
                    break
            if not load_converged:
                raise RuntimeError('STAGGERED_MAXIMUM_ITERATIONS_EXCEEDED')
        except HighDamageGate as exc:
            status = 'HIGH_DAMAGE_GATE_REACHED'
            event = dict(status=status, reason=str(exc), load_step=load_step, strain=strain,
                         sc=exc.cap, proposed_maximum_s=float(exc.update.s_new.max()),
                         last_accepted_maximum_s=float(accepted_s.max()),
                         proposed_state_accepted=False, continued_after_gate=False)
            write_json(step/'HIGH_DAMAGE_GATE.json', event)
            np.savez_compressed(step/'high_damage_state.npz', x_equilibrated_attempt=x,
                                x_last_accepted=accepted_x, s_last_accepted=accepted_s,
                                Y=before.W0, **{k:getattr(exc.update,k) for k in exc.update.__dataclass_fields__})
            break
        except (ValueError, FloatingPointError, RuntimeError, np.linalg.LinAlgError) as exc:
            status = 'FAILED'
            event = dict(status=status, reason=str(exc), exception=type(exc).__name__,
                         load_step=load_step, strain=strain, phase_update_after_failure=False)
            write_json(step/'FAILED.json', event); write_json(out/'FAILED.json', event)
            np.savez_compressed(step/'failed_state.npz', x_current=x, s_current=s_pf_bond,
                                x_last_accepted=accepted_x, s_last_accepted=accepted_s)
            break
    write_json(out/'LOAD_PATH.json', accepted); write_json(out/'STAGGERED_HISTORY.json', iterations)
    summary = dict(candidate=NAME, status=status, event=event, accepted_steps=len(accepted),
                   requested_steps=len(strains), Yc_Pa=yc, maximum_accepted_s=float(accepted_s.max()),
                   maximum_phase_energy_increase_J=max([0.]+[r['phase_energy_change_J'] for r in iterations]),
                   maximum_kinematic_change=max([0.]+[max(r['kinematic_changes'].values()) for r in iterations]),
                   minimum_accepted_phase_increment=min([0.]+[r['minimum_phase_increment'] for r in accepted]),
                   maximum_accepted_force_residual=max([0.]+[r['relative_force_residual'] for r in accepted]),
                   maximum_energy_balance_residual_J=max([0.]+[abs(r['energy_balance_residual_J']) for r in accepted]),
                   energy_balance_note='Trapezoidal external work diagnostic; not a closure gate or fitted correction',
                   high_damage_coupling='NOT_IMPLEMENTED', real_CFD='NOT_USED', physical_calibration='NOT_PERFORMED')
    write_json(out/'SUMMARY.json', summary)
    return summary
