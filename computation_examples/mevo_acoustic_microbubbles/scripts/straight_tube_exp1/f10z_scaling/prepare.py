"""One preregistered realization per N, constant nominal number density."""
import numpy as np
import h5py
from mevo_mb.coupling.pdms_cavity import apparatus
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from scripts.straight_tube_exp1.f10z.prepare import cloud, action
from .common import *
from .metrics import series

def main():
    BASE.mkdir(parents=True, exist_ok=True)
    if (BASE / 'PREREGISTRATION.json').exists():
        raise RuntimeError('Immutable preregistration already exists')
    cfg = dict(N=[32, 64, 128], seeds={'32': 18002, '64': 18064, '128': 18128},
        radius_m=R, frequency_Hz=2.5e6, flow='STRICT_ZERO',
        source_limit_Pa_equivalent=1500, local_nonzero_pressure_range_Pa=[0.1, 40000],
        half_width_rule_m='18e-6*(N/32)**(1/3)', minimum_initial_separation_m=5e-6,
        probe_s=1e-4, formal_windows_descending_s=[.002, .001, .0005],
        dt_s=1e-6, half_dt_s=.5e-6, per_run_budget_s=240, remote_budget_s=1200,
        per_N_projected_budget_s=520, max_trajectories=13,
        admission='Original step gates unchanged; full-window dt/2 and 1.5um cache comparisons, fixed IDs and edge events',
        cost_rule='Probe step time * window/dt; 1.35 safety factor, five equivalent runs plus 15s setup; <=520s per N, dt/2 <210s and remaining global budget. Else try next shorter registered window. Failure retained, no retry with looser tolerances.',
        actions=['OFF', 'PAIR'], labels=LABELS, training_ready=False,
        interpretation='Exploratory one realization per N, fixed nominal density, not a pure causal N comparison or statistically representative ensemble',
        model='Unchanged WaterSequential.step; original GPU forces/walls and CPU SQP/resistance')
    write(CONFIG, cfg)
    states = {}
    with h5py.File(HISTORICAL / 'INITIAL_STATES.h5') as h:
        states[32] = h['medium_train'][:]
    for n in [64, 128]:
        x = cloud(cfg['seeds'][str(n)], 18e-6 * (n / 32)**(1/3), 5e-6, n=n)
        states[n] = x - x.mean(0)
    cavity, _ = apparatus(read(OLD / 'CONFIG_SNAPSHOT.json')['geometry'])
    caches = {k: CachedWaterBasis(OLD / 'cache' / (k + '.h5')) for k in ['coarse', 'fine']}
    audit = {}
    with h5py.File(BASE / 'INITIAL_STATES.h5', 'w') as h:
        for n, x in states.items():
            m = series(np.array([0.]), x[None])
            wall = np.array([cavity.closest_wall(a)[0] for a in x]) - R
            assert wall.min() > 0 and m['min_distance_m'][0] >= 5e-6
            for cache in caches.values():
                coord = (x - cache.origin) / cache.spacing
                assert np.all(coord >= 1) and np.all(coord < np.array(cache.shape)-2)
            h.create_dataset(str(n), data=x)
            audit[str(n)] = dict(N=n, seed=cfg['seeds'][str(n)], half_width_m=18e-6*(n/32)**(1/3),
                number_density_m3=n/(36e-6*(n/32)**(1/3))**3,
                min_distance_m=m['min_distance_m'][0], Rg_m=m['rg_m'][0], components=int(m['components'][0]),
                min_wall_clearance_m=wall.min(), max_abs_coordinate_m=abs(x).max(), legal=True)
    write(BASE / 'INITIAL_STATE_AUDIT.json', audit)
    write(BASE / 'RUN_SCHEDULES.json', {a: action(a, .002).record() for a in cfg['actions']})
    # Regression fixture is an exact prefix of a completed historical trajectory.
    ref = next(HISTORICAL.glob('*/runs/medium_train__PAIR__2000us__DT1/trajectory.h5'))
    with h5py.File(ref) as h, h5py.File(BASE / 'N32_REGRESSION.h5', 'w') as out:
        t = h['time_s'][:]; keep = t <= cfg['probe_s'] + 1e-13
        out['time_s'] = t[keep]; out['positions_m'] = h['positions_m'][:][keep]
        out.attrs['historical_sha256'] = sha(ref)
    inputs = {str((OLD/p).relative_to(ROOT)): sha(OLD/p) for p in INPUTS}
    for p in [CONFIG, *[BASE/n for n in ['INITIAL_STATES.h5', 'INITIAL_STATE_AUDIT.json', 'RUN_SCHEDULES.json', 'N32_REGRESSION.h5']]]:
        inputs[str(p.relative_to(ROOT))] = sha(p)
    # All source files, not only the historically protected subset, stay unchanged.
    protected = {str(p.relative_to(ROOT)): sha(p) for p in (ROOT/'src/mevo_mb').rglob('*.py')}
    protected.update({p: sha(ROOT/p) for p in PROTECTED})
    write(BASE/'PREREGISTRATION.json', dict(protocol=cfg, pre_results=True,
        baseline_commit=identity()['execution_commit'], input_sha256=inputs, protected_source_sha256=protected))
    print(audit)

if __name__ == '__main__': main()
