"""Variable-N experiment I/O; every physical step calls unchanged WaterSequential.step."""
import time
import traceback
import resource
import numpy as np
import h5py
from mevo_mb.fast_sim.water_sequential import WaterSequential
from mevo_mb.fast_sim.nearfield import pair_reference
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.control.scheduled_pzt import Schedule
from scripts.straight_tube_exp1.f10s.common import context
from scripts.straight_tube_exp1.f10s.run import shorten
from scripts.straight_tube_exp1.f10z.run import Monitor
from .common import *

COLUMNS = ['accepted', 'residual_m', 'overlap_m', 'max_iteration_m', 'total_contact_m', 'wall_m',
    'full_step_s', 'solver_s', 'resistance_s', 'iterations', 'physical_path_sum_m',
    'contact_displacement_sum_m', 'primal_m', 'dual_N', 'force_residual_N',
    'complementarity_m', 'complementarity_N_m', 'actual_complementarity_m',
    'active_count', 'min_pressure_Pa', 'max_pressure_Pa', 'max_total_force_N', 'force_scale_N']

class Campaign:
    def __init__(self, out):
        verify()
        self.started = time.perf_counter(); self.out = out; self.cfg = read(CONFIG)
        self.p = policy(); self.c = context(); self.ident = identity(); self.records = []
        self.bases = {'coarse': self.c['basis'], 'fine': CachedWaterBasis(OLD/'cache/fine.h5')}
        with h5py.File(BASE/'INITIAL_STATES.h5') as h: self.states = {int(k): h[k][:] for k in h}
        self.schedules = {k: Schedule.from_record(v) for k, v in read(BASE/'RUN_SCHEDULES.json').items()}

    def engine(self, n, cache='coarse'):
        return WaterSequential(self.states[n], self.c['physics'], self.c['response'], self.bases[cache], self.c['cavity'], policy=self.p)

    def checks(self):
        evidence = {}
        for n in self.cfg['N']:
            sim = self.engine(n); x = self.states[n]; sim.control(self.schedules['PAIR'].segments[0].event())
            f = sim.fields(x)
            event = self.schedules['PAIR'].segments[0].event()
            pb, gb = self.bases['coarse'].query_basis(x)
            u = self.bases['coarse'].vector(event['amplitude_Pa'], event['phase_rad'], event['enabled'])
            pressure = u @ pb; gradient = np.einsum('i,inj->nj', u, gb)
            volume = self.c['response'].query(pressure)
            refs = dict(pressure_Pa=pressure, primary_force_N=-.5*np.real(volume[:, 0, None]*gradient.conj()),
                        secondary_force_N=pair_reference(x, volume, pressure, self.c['physics'], sim.sim.chi_zero))
            errors = {k: float(abs(f[k]-v).max()) for k, v in refs.items()}
            passed = all(errors[k] <= max(float(abs(v).max())*1e-9, 1e-24) for k, v in refs.items())
            wall = sim.wall_check(x)
            passed = passed and wall == dict(correction_m=0., error_bits=0) and not sim.sim.flow_enabled
            evidence[str(n)] = dict(pass_checks=bool(passed), errors=errors, wall=wall,
                pressure_range_Pa=[float(abs(pressure).min()), float(abs(pressure).max())],
                net_secondary_force_N=f['secondary_force_N'].sum(0), N=n)
        write(self.out/'CHECKS.json', evidence)
        assert all(r['pass_checks'] for r in evidence.values()), evidence

    def run(self, n, action, duration, variant='DT1'):
        if len(self.records) >= self.cfg['max_trajectories']:
            raise RuntimeError('Registered trajectory limit')
        dt = .5e-6 if variant == 'DT05' else 1e-6
        cache = 'fine' if variant == 'FINE' else 'coarse'
        name = f'N{n}_{action}_{round(duration*1e6)}us_{variant}'
        folder = self.out/'runs'/name; folder.mkdir(parents=True, exist_ok=False)
        schedule = shorten(self.schedules[action], duration)
        write(folder/'INPUT.json', dict(N=n, schedule=schedule.record(), dt_s=dt, cache=cache,
            identity=self.ident, preregistration_sha256=sha(BASE/'PREREGISTRATION.json')))
        begin = time.perf_counter(); x = [self.states[n].copy()]; times = [0.]
        rows, sampled, failed = [], [], []; last = None; reasons = []; active = None
        try:
            sim = self.engine(n, cache)
            for j, (t, h, segment) in enumerate(schedule.steps(dt)):
                if time.perf_counter()-begin > self.cfg['per_run_budget_s']:
                    reasons = ['PER_RUN_WALL_BUDGET']; break
                if time.perf_counter()-self.started > self.cfg['remote_budget_s']:
                    reasons = ['GLOBAL_WALL_BUDGET']; break
                if active != segment.segment_id:
                    sim.control(segment.event()); active = segment.segment_id
                last = sim.step(h); r = last['result']; f = last['fields']
                failed.append(','.join(last['reasons']))
                rows.append([last['accepted'], r['mobility_equation_residual_m'], r['geometric_overlap_m'].max(initial=0),
                    r['max_single_iteration_m'], r['max_total_contact_correction_m'], last['wall_correction_m'],
                    last['full_step_s'], r['total_s'], last['resistance_s'], r['qp_iterations'],
                    np.linalg.norm(r['prediction_m']-x[-1], axis=1).sum(), np.linalg.norm(r['contact_displacement_m'], axis=1).sum(),
                    r['primal_violation_m'], r['dual_negative_N'], r['force_equation_residual_N'],
                    r['complementarity_equivalent_gap_m'], r['complementarity_N_m'], r['actual_complementarity_equivalent_gap_m'],
                    r['active_count'], abs(f['pressure_Pa']).min(), abs(f['pressure_Pa']).max(),
                    np.linalg.norm(f['total_force_N'], axis=1).max(), r['force_scale_N']])
                if j % max(1, round(1e-4/dt)) == 0 or not last['accepted']:
                    sampled.append(dict(time_s=t, positions_m=x[-1].copy(), **f,
                        contact_force_N=r['contact_force_N'], pair_contact_reaction_N=r['pair_force_N']))
                if not last['accepted']:
                    reasons = last['reasons']; break
                x.append(sim.x.copy()); times.append(t+h)
                if j and j % 500 == 0:
                    print(name, 'progress_ms', times[-1]*1e3, 'elapsed_s', time.perf_counter()-begin, flush=True)
        except Exception as ex:
            reasons = ['EXCEPTION: '+repr(ex)]
            (folder/'exception.txt').write_text(traceback.format_exc())
        integration_s = time.perf_counter()-begin
        if last is not None: write(folder/'LAST_STEP.json', last)
        rows = np.asarray(rows).reshape(-1, len(COLUMNS)); x = np.asarray(x)
        with h5py.File(folder/'trajectory.h5', 'w') as h:
            for k, v in dict(positions_m=x, time_s=times, particle_id=np.arange(n), step_metrics=rows).items():
                h.create_dataset(k, data=v, compression='gzip')
            h['step_metrics'].attrs['columns'] = ','.join(COLUMNS)
            h['step_failed_gates'] = np.array(failed, dtype=h5py.string_dtype())
            group = h.create_group('force_samples')
            if sampled:
                for k in sampled[0]: group.create_dataset(k, data=np.asarray([v[k] for v in sampled]), compression='gzip')
            if reasons and last is not None: h['failure_proposal_m'] = last['result']['positions_m']
            h.attrs.update(flow='STRICT_ZERO', N=n, execution_commit=self.ident['execution_commit'],
                schedule_json=__import__('json').dumps(schedule.record()), training_ready=False,
                position_output='Every accepted step; no interpolation, padding, projection or clustering simplification')
        result = dict(name=name, N=n, action=action, variant=variant, dt_s=dt, cache=cache,
            status='PASS' if not reasons and abs(times[-1]-duration)<1e-13 else 'FAIL', reasons=reasons,
            requested_s=duration, completed_s=times[-1], accepted_steps=len(x)-1, attempted_steps=len(rows),
            displacement_m=x[-1].mean(0)-x[0].mean(0), wall_s=time.perf_counter()-begin, integration_s=integration_s,
            mean_step_s=float(rows[:, 6].mean()) if len(rows) else None,
            solver_s=float(rows[:, 7].sum()), resistance_s=float(rows[:, 8].sum()),
            rss_peak_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            trajectory_sha256=sha(folder/'trajectory.h5'), labels=LABELS, training_ready=False)
        result['max_metrics'] = {k: float(rows[:, i].max(initial=0)) for i, k in enumerate(COLUMNS) if i not in [0, 19]}
        result['pressure_min_Pa'] = float(rows[:, 19].min()) if len(rows) else None
        write(folder/'RESULT.json', result); self.records.append(result); write(self.out/'RUNS.json', self.records)
        print(name, result['status'], 'completed_ms', times[-1]*1e3, 'wall_s', result['wall_s'], reasons, flush=True)
        return result

def main():
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument('--label', required=True); a = ap.parse_args()
    out = BASE/a.label; out.mkdir(parents=True, exist_ok=False)
    begin = time.perf_counter(); c = None; monitor = Monitor(); monitor.thread.start(); status = 'EXCEPTION'
    try:
        c = Campaign(out); c.checks()
        r = c.run(32, 'PAIR', c.cfg['probe_s'])
        with h5py.File(out/'runs'/r['name']/'trajectory.h5') as h, h5py.File(BASE/'N32_REGRESSION.h5') as old:
            error = float(abs(h['positions_m'][:]-old['positions_m'][:]).max())
            matched = bool(np.allclose(h['time_s'][:], old['time_s'][:], atol=1e-14, rtol=0))
        write(out/'N32_REGRESSION.json', dict(max_coordinate_error_m=error, same_times=matched,
            passed=bool(r['status']=='PASS' and error<1e-18 and matched)))
        assert r['status']=='PASS' and error<1e-18 and matched, 'Historical 32-particle regression failed'
        probes = {}
        for n in [64, 128]:
            c.run(n, 'OFF', c.cfg['probe_s'])
            probes[n] = c.run(n, 'PAIR', c.cfg['probe_s'])
        decisions = []
        for n in [64, 128]:
            probe = probes[n]; chosen = None
            if probe['status']=='PASS':
                for duration in c.cfg['formal_windows_descending_s']:
                    single = probe['mean_step_s']*duration/c.cfg['dt_s']
                    projected = 1.35*single*5 + 15
                    eligible = projected<c.cfg['per_N_projected_budget_s'] and 1.35*single*2<210 and projected<c.cfg['remote_budget_s']-(time.perf_counter()-c.started)
                    decisions.append(dict(N=n, window_s=duration, projected_wall_s=projected, eligible=eligible))
                    if eligible: chosen=duration; break
            if chosen is None:
                decisions.append(dict(N=n, status='FORMAL_NOT_RUN', reason='Failed probe or registered cost limit'))
            write(out/'DECISIONS.json', decisions)
            if chosen is not None:
                off = c.run(n, 'OFF', chosen)
                for variant in ['DT1', 'DT05', 'FINE']:
                    r = c.run(n, 'PAIR', chosen, variant)
                    if r['status']!='PASS':
                        decisions.append(dict(N=n, status='REFINEMENTS_STOPPED', reason=r['reasons']))
                        write(out/'DECISIONS.json', decisions); break
        status = 'COMPLETED'
    except Exception:
        (out/'EXCEPTION.txt').write_text(traceback.format_exc()); raise
    finally:
        monitor.done.set(); monitor.thread.join(4)
        write(out/'resources.json', monitor.rows)
        write(out/'CAMPAIGN.json', dict(status=status, wall_s=time.perf_counter()-begin,
            identity=identity(), runs=c.records if c else [], CPU_user_s=resource.getrusage(resource.RUSAGE_SELF).ru_utime,
            CPU_system_s=resource.getrusage(resource.RUSAGE_SELF).ru_stime))

if __name__ == '__main__': main()
