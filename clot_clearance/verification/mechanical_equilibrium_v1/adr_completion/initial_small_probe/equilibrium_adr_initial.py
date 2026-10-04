"""Kilic & Madenci (2010), DOI 10.1016/j.tafmec.2010.08.001, eqs 11–19.

Additive continuation of equilibrium.py: that historical implementation and
its tests remain unchanged. ADR is a NUMERICAL EQUILIBRATION DEVICE.
The NOSB adaptation of eq14 uses the existing total-force numerical tangent;
it does not use the paper's bond-model-specific eq15. See METHOD_AUDIT.md.
"""
from dataclasses import dataclass
import numpy as np
from . import equilibrium as legacy
from .erosion import SettlingMonitor

make_fixture = legacy.make_fixture
FrozenForceEvaluator = legacy.FrozenForceEvaluator


@dataclass
class ADRResult(legacy.EquilibriumResult):
    fictitious_diagonal: np.ndarray | None = None
    pseudo_half_velocity: np.ndarray | None = None


def tangent_row_sum(cloud, x, evaluate, perturbation):
    """Absolute rows of -dR_free/dx_free; stream columns without dense storage."""
    free = ~cloud.fixed
    indices = np.flatnonzero(np.repeat(free, 3))
    rows = np.zeros(len(indices))
    for index in indices:
        plus, minus = x.copy(), x.copy()
        plus.flat[index] += perturbation
        minus.flat[index] -= perturbation
        column = (evaluate(plus).residual[free].ravel() -
                  evaluate(minus).residual[free].ravel()) / (2 * perturbation)
        rows += abs(column)
    if not np.isfinite(rows).all() or np.any(rows <= 0):
        raise FloatingPointError('Nonpositive/nonfinite tangent row sum')
    return rows.reshape((-1, 3))


def adaptive_damping(u, residual, previous_residual, mass, half_velocity, dtau):
    """Eqs18–19, zero-velocity rule from paragraph immediately after eq19."""
    local = np.divide(-(residual - previous_residual) / mass, dtau * half_velocity,
                      out=np.zeros_like(u), where=half_velocity != 0)
    denominator = float(np.sum(u*u))
    quotient = float(np.sum(u*u*local)/denominator) if denominator else 0.
    if not np.isfinite(quotient) or quotient < 0:
        raise FloatingPointError(f'Undefined ADR Rayleigh quotient: {quotient}')
    return 2*np.sqrt(quotient)


def equilibrate(cloud, x0, v0, integrity, config, force_callback=None,
                method=None, verification_only=False):
    method = method or config['mechanical_equilibrium']['method']
    if method != 'adaptive_dynamic_relaxation':
        return legacy.equilibrate(cloud, x0, v0, integrity, config, force_callback,
                                  method, verification_only)
    c = config; ec = c['mechanical_equilibrium']; ac = ec['adr']
    # Reuse all state/tolerance guards, replacing only historical D=0 restriction.
    legacy.validate_options(cloud, x0, v0, np.ones(len(cloud.pairs)), c,
                            'legacy_fixed_damping', verification_only)
    g = np.array(integrity, float, copy=True)
    if g.shape != (len(cloud.pairs),) or not np.isfinite(g).all() or np.any((g <= 0) | (g > 1)):
        raise ValueError('ADR development path requires intact attached topology, 0 < integrity <= 1')
    g.setflags(write=False)
    for key in ['pseudo_time_step', 'mass_row_sum_factor', 'finite_difference_step_over_spacing',
                'pseudo_kinetic_to_stored_tolerance']:
        if not np.isfinite(ac[key]) or ac[key] <= 0: raise ValueError(f'Invalid ADR {key}')
    if ac['mass_row_sum_factor'] <= 1: raise ValueError('Strict eq14 stability margin required')
    extra = ac['additional_iterations_after_convergence']
    if isinstance(extra, bool) or int(extra) != extra or extra < 0: raise ValueError('Invalid extra iteration count')
    if np.any(v0 != 0): raise ValueError('ADR starts with zero pseudo velocity, not physical incoming velocity')
    evaluate = force_callback or FrozenForceEvaluator(cloud, g, c)
    x = np.array(x0, copy=True); free = ~cloud.fixed
    dtau = ac['pseudo_time_step']; perturbation = cloud.spacing*ac['finite_difference_step_over_spacing']
    row_sum = tangent_row_sum(cloud, x, evaluate, perturbation)
    mass = ac['mass_row_sum_factor']*.25*dtau**2*row_sum
    half = np.zeros_like(x); state = evaluate(x); previous = None
    history = []; checks = []; consecutive = 0; first_pass = None; failure = None
    monitor = SettlingMonitor(g > 0, ec['driver_absolute_floor'])

    def capture(iteration, damping):
        row = legacy.diagnostics(cloud, x, state, c)
        pseudo = .5*float(np.sum(mass*half[free]**2))
        row.update(equilibrium_iteration=iteration, equilibrium_pseudo_time=iteration*dtau,
                   mechanical_time_s=None, phase='FROZEN_DAMAGE_EQUILIBRATION',
                   maximum_damage=float((1-g).max()), driver_drift=None,
                   adr_damping_coefficient=float(damping), adr_pseudo_kinetic_measure=pseudo,
                   adr_pseudo_kinetic_to_stored_ratio=pseudo/max(abs(state.stored_energy_J), ec['energy_floor_J']),
                   represented_ultrasound_cycles=0.)
        history.append(row)
        # The old monitor is used ONLY for strain ranges, not its time/energy fields.
        monitor.add(iteration, iteration, state.epsilon_plus, 0.)
        return row

    capture(0, 0.)
    for iteration in range(1, ec['maximum_iterations']+1):
        residual = state.residual[free]
        if iteration == 1:
            damping = 0.; half[free] = dtau*residual/mass/2  # (13)
        else:
            try:
                damping = adaptive_damping((x-cloud.X)[free], residual, previous,
                                           mass, half[free], dtau)
            except FloatingPointError as exc:
                failure = str(exc); break
            half[free] = ((2-damping*dtau)*half[free]+2*dtau*residual/mass)/(2+damping*dtau)  # (11)
        previous = residual.copy()
        candidate = x + dtau*half  # (12)
        candidate[cloud.fixed] = cloud.X[cloud.fixed]
        try:
            next_state = evaluate(candidate)
        except (FloatingPointError, ValueError) as exc:
            failure = str(exc); break
        x = candidate; state = next_state
        row = capture(iteration, damping)
        if iteration % ec['check_interval'] == 0:
            drift = monitor.summary(ec['driver_drift_tolerance'])['driver_relative_drift']
            row['driver_drift'] = drift
            okay = (row['relative_force_residual'] <= ec['relative_force_residual_tolerance'] and
                    drift <= ec['driver_drift_tolerance'] and
                    row['adr_pseudo_kinetic_to_stored_ratio'] <= ac['pseudo_kinetic_to_stored_tolerance'])
            consecutive = consecutive+1 if okay else 0
            checks.append(dict(iteration=iteration, relative_force_residual=row['relative_force_residual'],
                               driver_drift=drift, pseudo_kinetic_ratio=row['adr_pseudo_kinetic_to_stored_ratio'],
                               all_criteria_passed=bool(okay), consecutive_passes=consecutive))
            if consecutive >= ec['minimum_consecutive_converged_checks'] and first_pass is None:
                first_pass = iteration
            monitor = SettlingMonitor(g > 0, ec['driver_absolute_floor'])
            monitor.add(iteration, iteration, state.epsilon_plus, 0.)
            if first_pass is not None and iteration >= first_pass+extra and consecutive >= ec['minimum_consecutive_converged_checks']:
                break
    bound_ratio = None
    if ac['final_mass_bound_audit']:
        end_sum = tangent_row_sum(cloud, x, evaluate, perturbation/2)
        bound_ratio = float(np.max(.25*dtau**2*end_sum/mass))
        if bound_ratio > 1.: failure = 'Endpoint tangent violates eq14 fictitious mass bound'
    passed = bool(failure is None and checks and checks[-1]['iteration'] == history[-1]['equilibrium_iteration']
                  and consecutive >= ec['minimum_consecutive_converged_checks'])
    summary = dict(status='EQUILIBRIUM_CRITERIA_PASSED' if passed else 'INCOMPLETE_MECHANICAL_EQUILIBRIUM',
                   method=method, converged=passed, iterations=history[-1]['equilibrium_iteration'],
                   first_converged_iteration=first_pass, mechanical_time_s=None,
                   equilibrium_pseudo_time=history[-1]['equilibrium_pseudo_time'], erosion_exposure_time_s=0.,
                   represented_ultrasound_cycles=0., damage_advanced=False, fixed_topology=True,
                   physical_density_unchanged=True, physical_density_kg_m3=c['clot']['density_kg_m3'],
                   fictitious_mass_min=float(mass.min()), fictitious_mass_max=float(mass.max()),
                   fictitious_mass_units='N * pseudo_time^2 / m; not physical kg',
                   endpoint_mass_bound_ratio=bound_ratio, failure_reason=failure, final=history[-1],
                   algorithm_source='Kilic & Madenci (2010), eq11-14,18-19; numerical NOSB tangent for eq14',
                   convergence_scope='Free force + pseudo kinetic ratio + tensile driver drift; numerical device only')
    full_mass = np.zeros_like(x); full_mass[free] = mass
    return ADRResult(x, None, method, passed, history, checks, summary, full_mass, half)
