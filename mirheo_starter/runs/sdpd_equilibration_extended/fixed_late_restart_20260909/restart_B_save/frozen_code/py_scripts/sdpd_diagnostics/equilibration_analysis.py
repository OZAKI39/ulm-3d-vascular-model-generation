"""Fixed-window screening, with trend, sampling and thermal agreement separated."""
import math
import numpy as np
from py_scripts.fluid_physics.analysis import block_statistics, correlation_time_samples, t95


def empty_result(status='READY_TO_RUN_AWAITING_BUDGET', reason='NOT_RUN_AWAITING_AUTHORIZATION'):
    return {
        'summary': {'status':status, 'execution_status':reason, 'program_completed':None,
                    'liquid_stationary':None, 'stable_temperature_matches_target':None,
                    'temperature_mean_star':None, 'temperature_mean_K':None,
                    'temperature_CI95_star':None, 'actual_steps':None,
                    'actual_time_star':None, 'actual_physical_time_s':None,
                    'selection':None, 'all_liquid_properties_qualified':False,
                    'viscosity_new_measurement':None, 'EOS_new_measurement':None,
                    'human_review':'PENDING', 'reason':reason},
        'stationarity': {'status':reason, 'temperature':None, 'pressure':None,
                         'structure':None, 'momentum':None, 'conservation':None},
        'sampling': {'status':reason, 'total_samples':None, 'formal_samples':None,
                     'block_count':None, 'block_samples':None, 'discarded_tail_samples':None,
                     'measured_stationary_ACF_samples':None, 'CI_basis':None},
        'statistics':None}


def completed(execution, completion, plan):
    return (execution.get('status') == 'COMPLETED' and not execution.get('timeout')
            and not execution.get('remaining_own_group_processes')
            and completion.get('status') == 'COMPLETED_PLANNED_STEPS'
            and completion.get('steps') == plan['steps']
            and math.isclose(completion.get('time_star',-1), plan['duration_star'], abs_tol=1e-12))


def metric_rules(plan):
    s, c = plan['parameters'], plan['criteria']
    n, target = s['task']['n_star'], s['kBT_star']
    rules = {
        'kBT_COM_star': (target, c['temperature_two_half_drift_relative'], c['temperature_quarter_range_relative']),
        'pressure_star': (n*target, c['pressure_drift_thermal_pressure_units'], c['pressure_quarter_range_thermal_pressure_units']),
        'kernel_number_mean': (n, c['kernel_mean_drift_over_global_n'], 2*c['kernel_mean_drift_over_global_n']),
        'kernel_number_variance': (n*n, c['kernel_variance_drift_over_global_n_squared'], 2*c['kernel_variance_drift_over_global_n_squared'])}
    for q in ['kernel_q10','kernel_q50','kernel_q90']:
        rules[q] = (n, c['kernel_quantile_drift_over_global_n'], 2*c['kernel_quantile_drift_over_global_n'])
    for axis in 'xyz':
        rules['mean_v'+axis] = (math.sqrt(target/s['m_star']), c['COM_velocity_drift_over_thermal_speed'], 2*c['COM_velocity_drift_over_thermal_speed'])
    return rules


def sufficient_statistics(t, x, minimum, block_samples, min_blocks):
    r = block_statistics(t, x, min_duration=minimum, min_blocks=min_blocks, fixed_block_samples=block_samples)
    if r['status'] != 'SUFFICIENT':
        r['ci95_halfwidth'] = None
        r['ci_method'] = 'No qualified CI: too few complete independent blocks.'
    return r


def analyze_series(data, plan, execution, completion, structure):
    """Synthetic callers must keep their evidence outside real result packages."""
    r = empty_result('STATIONARITY_OR_SAMPLING_INCONCLUSIVE','NOT_EVALUATED')
    summary, checks, audit = r['summary'], r['stationarity'], r['sampling']
    s, c = plan['parameters'], plan['criteria']
    steps = int(completion.get('steps',0))
    summary.update(actual_steps=steps, actual_time_star=steps*s['dt_star'],
                   actual_physical_time_s=steps*s['dt_star']*s['locked_units']['t0'],
                   program_completed=execution.get('exit_code') == 0)
    if not completed(execution, completion, plan):
        partial = (execution.get('timeout') or completion.get('status') == 'COOPERATIVE_PARTIAL')
        summary.update(status='STATIONARITY_OR_SAMPLING_INCONCLUSIVE' if partial else 'NUMERICAL_OR_RUNTIME_FAILURE',
                       execution_status='PARTIAL' if partial else 'FAILED',
                       reason='Required continuous physical interval was not completed; exit code zero is insufficient.')
        return r
    summary['execution_status'] = 'COMPLETED_FULL_EXPERIMENT'
    try:
        t = np.asarray(data['time_star'],float)
        n = np.asarray(data['N'],float)
        if (len(t) < 2 or np.any(np.diff(t) <= 0)
                or not np.allclose(np.asarray(data['step'])*s['dt_star'],t,atol=1e-12,rtol=0)
                or not np.allclose(np.diff(t),plan['sample_interval_star'],atol=1e-12,rtol=0)
                or t[0] != 0 or not math.isclose(t[-1],plan['duration_star'],abs_tol=1e-12)):
            raise ValueError('STEP_TIME_OR_CADENCE_MISMATCH')
        mass = np.asarray(data['total_mass_star'],float)
        if not np.all(n == s['expected_N']) or not np.all(mass == s['expected_N']*s['m_star']):
            raise ValueError('PARTICLE_COUNT_OR_MASS_NOT_CONSERVED')
        checks['conservation'] = {'status':'PASS', 'N':s['expected_N'], 'total_mass_star':float(mass[0]),
                                  'scope':'Constructed periodic global density is not evidence of local structural equilibrium.'}
        rules = metric_rules(plan)
        if any(not np.isfinite(np.asarray(data[k])[1:]).all() for k in rules):
            raise ValueError('NONFINITE_MEASUREMENT')
        if not np.allclose(np.asarray(data['density_time_star'])[1:],t[1:]-s['dt_star'],atol=1e-12,rtol=0):
            raise ValueError('DENSITY_SAMPLING_PHASE_MISMATCH')
        lo, hi = plan['formal_window_star']
        mid = (lo+hi)/2
        selected = (t > lo+1e-12) & (t <= hi+1e-12)
        ts = t[selected]
        audit.update(total_samples=len(t), startup_samples=int(np.sum(t <= lo+1e-12)),
                     formal_samples=len(ts), formal_window_star=[lo,hi],
                     endpoint_convention='Right endpoints: (0.30,0.40]; fixed before execution.',
                     minimum_observation_duration_star=c['min_observation_duration_star'],
                     measured_stationary_ACF_samples=None,
                     historical_ACF_role='Old drifting 82.44-sample estimate is planning context only; never copied into new measured ACF.')
        if hi-lo+1e-12 < c['min_observation_duration_star'] or len(ts) < 16:
            summary['reason'] = 'Insufficient complete formal interval.'
            return r
        checks.update(status='POINT_TRENDS_FIRST', structure=structure)
        trends = {}
        for key, (scale, limit, quarter_limit) in rules.items():
            x = np.asarray(data[key],float)[selected]
            a, b = x[ts <= mid+1e-12], x[ts > mid+1e-12]
            quarter_values = [x[(ts > left+1e-12) & (ts <= right+1e-12)]
                              for left, right in zip(np.linspace(lo,hi,5)[:-1],np.linspace(lo,hi,5)[1:])]
            quarters = [float(q.mean()) for q in quarter_values]
            residual_spread = float(np.sqrt(np.mean(np.concatenate([q-q.mean() for q in quarter_values])**2))/scale)
            differences = np.diff(quarters)
            drift = float(abs(b.mean()-a.mean())/scale)
            monotonic = bool(np.all(differences > 0) or np.all(differences < 0))
            trends[key] = {'normalization_scale_star':scale, 'two_half_limit':limit,
                           'two_half_point_drift':drift, 'quarter_means':quarters,
                           'quarter_range':float(np.ptp(quarters)/scale), 'quarter_range_limit':quarter_limit,
                           'within_quarter_RMS_over_scale':residual_spread,
                           'trend_noise_rule':'Require a coherent material half-window change larger than 3 times the within-quarter RMS; no independent-sample assumption or steady-state CI.',
                           'coherent_drift':bool(monotonic and drift > c['coherent_drift_over_within_quarter_rms_min']*residual_spread and (drift > limit or np.ptp(quarters)/scale > quarter_limit)),
                           'drift_CI95_halfwidth':None, 'drift_upper_bound':None,
                           'status':'POINT_SCREEN_ONLY'}
        checks['metrics'] = trends
        if structure.get('measurement_status') != 'PASS':
            raise ValueError('STATE_SNAPSHOT_OR_CHANNEL_AUDIT_FAILED')
        if any(v['coherent_drift'] for v in trends.values()) or structure.get('status') == 'TRANSIENT_PERSISTS':
            summary.update(status='TRANSIENT_PERSISTS', liquid_stationary=False,
                           reason='Material coherent trend remains in fixed quarters/halves or fixed structural snapshots. No steady-state ACF or stable-temperature bias reported.')
            checks['status'] = 'TRANSIENT_PERSISTS'
            # A conditional diagnostic, never an automatic extension or convergence guarantee.
            q = trends['kBT_COM_star']['quarter_means']
            last_slope = (q[-1]-q[-2])/((hi-lo)/4)
            audit['next_observation_condition'] = {'terminal_temperature_slope_per_star':last_slope,
                'additional_time_to_target_at_constant_last_slope_star': max(0.,(s['kBT_star']-q[-1])/last_slope) if last_slope else None,
                'warning':'Linear terminal-slope extrapolation only; it can miss a biased asymptote and is not a convergence forecast.'}
            return r
        if (any(v['two_half_point_drift'] > v['two_half_limit'] or v['quarter_range'] > v['quarter_range_limit'] for v in trends.values())
                or structure.get('status') != 'PASS'):
            summary['reason'] = 'Fixed point trends or structural evidence are inconclusive; steady-state ACF not asserted.'
            checks['status'] = 'INCONCLUSIVE'
            return r
        minimum = plan['minimum_block_duration_star']
        spacing = plan['sample_interval_star']
        tau = {}
        for key in rules:
            x = np.asarray(data[key],float)[selected]
            tau[key] = {'full':correlation_time_samples(x),
                        'first_half':correlation_time_samples(x[ts <= mid+1e-12]),
                        'second_half':correlation_time_samples(x[ts > mid+1e-12])}
        block = max(math.ceil(minimum/spacing), math.ceil(c['ACF_multiplier']*max(v for x in tau.values() for v in x.values())))
        statistics = {}
        for key, (scale, limit, _) in rules.items():
            x = np.asarray(data[key],float)[selected]
            stats = sufficient_statistics(ts,x,minimum,block,c['min_blocks'])
            halves = []
            for mask in [ts <= mid+1e-12, ts > mid+1e-12]:
                halves.append(sufficient_statistics(ts[mask],x[mask],minimum,block,c['min_half_blocks']))
            stats['fixed_halves'] = halves
            statistics[key] = stats
            if stats['status'] == 'SUFFICIENT' and all(h['status'] == 'SUFFICIENT' for h in halves):
                delta = abs(halves[1]['mean']-halves[0]['mean'])/scale
                # Each mean and its CI use identical complete blocks in its fixed half.
                error = sum(h['ci95_halfwidth'] for h in halves)/scale
                trends[key].update(two_half_complete_block_drift=delta, drift_CI95_halfwidth=error,
                                   drift_upper_bound=delta+error, status='PASS' if delta+error <= limit else 'INCONCLUSIVE',
                                   CI_basis='Sum of the two half-mean 95% t halfwidths; conservative interval, no independence assumption or joint coverage claim.')
            else:
                trends[key]['status'] = 'INSUFFICIENT_BLOCKS'
        r['statistics'] = statistics
        temp = statistics['kBT_COM_star']
        sufficient = all(v['status'] == 'SUFFICIENT' and all(h['status'] == 'SUFFICIENT' for h in v['fixed_halves']) for v in statistics.values())
        audit.update(status='SUFFICIENT' if sufficient else 'INSUFFICIENT_BLOCKS',
                     block_samples=block, block_count=temp['block_count'],
                     discarded_tail_samples=temp['discarded_tail_samples'],
                     estimator_sample_count=temp.get('estimator_sample_count'),
                     block_duration_star=temp.get('block_duration_star'),
                     measured_stationary_ACF_samples=tau,
                     ACF_scope='Estimated only after fixed point-trend screens; conditional on stationarity, not independent proof.',
                     CI_basis=temp['ci_method'], statistics=statistics,
                     additional_formal_samples_condition=max(0,c['min_blocks']*block-len(ts),2*c['min_half_blocks']*block-len(ts)),
                     additional_duration_star_condition=max(0,c['min_blocks']*block-len(ts),2*c['min_half_blocks']*block-len(ts))*spacing)
        if not sufficient or any(v['status'] != 'PASS' for v in trends.values()):
            summary['reason'] = 'Point trends are small, but independent blocks or drift uncertainty do not establish stationarity.'
            checks['status'] = 'INCONCLUSIVE'
            return r
        checks.update(status='PASS', temperature=trends['kBT_COM_star'], pressure=trends['pressure_star'],
                      momentum={k:trends[k] for k in trends if k.startswith('mean_v')})
        mean, hw, target = temp['mean'], temp['ci95_halfwidth'], s['kBT_star']
        error = abs(mean/target-1)
        summary.update(liquid_stationary=True, temperature_mean_star=mean,
                       temperature_mean_K=mean/target*s['target_temperature_K'],
                       temperature_CI95_star=[mean-hw,mean+hw],
                       temperature_CI95_K=[(mean-hw)/target*s['target_temperature_K'],(mean+hw)/target*s['target_temperature_K']],
                       thermal_error_including_uncertainty=error+hw/target)
        if error+hw/target <= c['temperature_relative_error_including_CI']:
            summary.update(status='EQUILIBRIUM_SCREEN_PASS', stable_temperature_matches_target=True,
                           reason='Stationarity and the 2% thermal screen pass for this one box, parameter set, dt and interval only.')
        elif error-hw/target > c['temperature_relative_error_including_CI']:
            summary.update(status='STATIONARY_BUT_THERMALLY_BIASED', stable_temperature_matches_target=False,
                           reason='Stable mean and its CI lie outside the 2% band. Next: unforced dt/dt2 comparison of stable states.')
        else:
            summary.update(reason='Stationary evidence exists, but thermal uncertainty overlaps the 2% decision boundary.')
        return r
    except (ValueError,KeyError,IndexError,FloatingPointError) as exc:
        summary.update(status='MEASUREMENT_OR_ANALYSIS_ERROR',reason=str(exc),liquid_stationary=None)
        checks['status'] = 'MEASUREMENT_OR_ANALYSIS_ERROR'
        return r
