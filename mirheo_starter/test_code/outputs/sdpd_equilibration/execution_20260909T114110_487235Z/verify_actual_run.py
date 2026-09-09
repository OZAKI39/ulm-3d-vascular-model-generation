"""Independent CPU cross-check of this authorized run; never launches a worker."""
import csv
import datetime
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path('/home/lzy/projects/mirheo_starter')
AUDIT = Path(__file__).resolve().parent
RAW = ROOT / 'runs/sdpd_equilibration/unforced_plateau_20260909/longer_unforced_thermal_plateau'


def read(path):
    return json.loads(path.read_text())


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    spec = read(RAW / 'actual_parameters.json')
    plan = read(RAW / 'run_plan.json')
    execution = read(RAW / 'execution.json')
    completion = read(RAW / 'worker_completion.json')
    observer = read(RAW / 'equilibration_observer.json')
    with (RAW / 'raw_statistics.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    data = {k: np.array([float(r[k]) if r[k] else np.nan for r in rows]) for k in rows[0]}
    checks = {}
    checks['registered_parameters'] = spec == plan['parameters']
    checks['approved_plan'] = plan['plan_sha256'] == read(AUDIT / 'user_confirmation.json')['scope']['plan_sha256']
    checks['task_deadline'] = execution['elapsed_monotonic_s'] <= 600
    checks['full_completion'] = (completion['steps'] == 400000 and completion['status'] == 'COMPLETED_PLANNED_STEPS'
                                 and execution['status'] == 'COMPLETED' and not execution['timeout'])
    checks['continuous_cadence'] = bool(np.array_equal(data['step'], np.arange(0, completion['steps'] + 1, 200)))
    checks['continuous_time'] = bool(np.allclose(data['time_star'], data['step'] * 1e-6, atol=1e-12, rtol=0))
    checks['one_initialization'] = observer['continuous_coordinator_count'] == observer['velocity_initialization_count'] == 1
    checks['particle_and_mass_conservation'] = bool(np.all(data['N'] == 4096) and np.all(data['total_mass_star'] == 4096))
    checks['no_drive_or_extra_thermostat'] = (spec['task']['force_star'] == 0 and not observer['additional_thermostat']
                                            and not observer['velocity_rescaling_during_run'])
    initial_ids = np.load(RAW / 'initial_ids.npy', allow_pickle=False)
    snapshot_rows = []
    mass, volume = spec['m_star'], float(np.prod(spec['domain_star']))
    for f in observer['snapshots']:
        prefix = RAW / f"snapshot_{f['index']:02d}"
        v = np.asarray(np.load(str(prefix) + '_post_velocities.npy', allow_pickle=False), dtype=np.float64)
        ids = np.load(str(prefix) + '_ids.npy', allow_pickle=False)
        stress = np.asarray(np.load(str(prefix) + '_stresses.npy', allow_pickle=False), dtype=np.float64)
        centered = v - v.mean(axis=0)
        thermal = mass * float(np.sum(centered * centered)) / (3 * (len(v) - 1))
        raw_thermal = mass * float(np.sum(v * v)) / (3 * len(v))
        pressure = float(np.sum(stress[:, [0, 3, 5]])) / (3 * volume) + mass * float(np.sum(centered * centered)) / (3 * volume)
        i = int(np.flatnonzero(data['step'] == f['step'])[0])
        deltas = {'COM_temperature_star': abs(thermal - data['kBT_COM_star'][i]),
                  'raw_temperature_star': abs(raw_thermal - data['kBT_raw_star'][i]),
                  'pressure_star': abs(pressure - data['pressure_star'][i])}
        pre_v = np.asarray(np.load(str(prefix) + '_pre_velocities.npy', allow_pickle=False), float)
        force = np.asarray(np.load(str(prefix) + '_pre_forces.npy', allow_pickle=False), float)
        pre_x = np.asarray(np.load(str(prefix) + '_positions.npy', allow_pickle=False), float)
        post_x = np.asarray(np.load(str(prefix) + '_post_positions.npy', allow_pickle=False), float)
        kick = v - pre_v - spec['dt_star'] * force / mass
        drift = post_x - pre_x - spec['dt_star'] * v
        domain = np.asarray(spec['domain_star'])
        drift -= domain * np.rint(drift / domain)
        phase = {'kick_max_absolute': float(np.max(abs(kick))), 'drift_PBC_max_absolute': float(np.max(abs(drift)))}
        ok = bool(max(deltas.values()) < 1e-9 and np.array_equal(np.sort(ids), np.sort(initial_ids)))
        ok = ok and phase['kick_max_absolute'] <= plan['criteria']['phase_kick_absolute_tolerance'] and phase['drift_PBC_max_absolute'] <= plan['criteria']['phase_drift_absolute_tolerance']
        snapshot_rows.append({'step': f['step'], 'passed': ok, 'absolute_recomputation_errors': deltas, 'same_phase_recomputation': phase})
    checks['snapshot_temperature_pressure_ID_recomputation'] = all(r['passed'] for r in snapshot_rows)
    from py_scripts.fluid_physics.runner import group_members
    process = read(RAW / 'process_identity.json')
    remaining = group_members(process['pgid'])
    checks['no_remaining_task_processes'] = not remaining and not execution['remaining_own_group_processes']
    before = read(AUDIT / 'authorized_preflight.json')['budget']['budget_snapshot']
    checks['old_ledgers_unchanged'] = all(digest(Path(name)) == h for name, h in before['ledger_sha256'].items())
    lo, hi = plan['formal_window_star']
    formal = (data['time_star'] > lo + 1e-12) & (data['time_star'] <= hi + 1e-12)
    descriptive = {}
    for k in ['kBT_COM_star', 'temperature_K', 'pressure_pa', 'kernel_number_mean', 'kernel_number_variance', 'mean_vx', 'mean_vy', 'mean_vz']:
        values = data[k][formal]
        descriptive[k] = {'sample_count': int(len(values)), 'mean': float(values.mean()) if len(values) else None,
                          'min': float(values.min()) if len(values) else None, 'max': float(values.max()) if len(values) else None}
    result = {'status': 'PASS' if all(checks.values()) else 'REVIEW_REQUIRED',
              'created_at': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'checks': checks,
              'raw_directory': str(RAW), 'raw_manifest_sha256': digest(RAW / 'output_sha256.json'),
              'snapshot_recomputation': snapshot_rows, 'remaining_own_processes': remaining,
              'startup_peak_temperature_ratio': float(np.max(data['kBT_COM_star']) / spec['kBT_star']),
              'formal_window_descriptive_only': descriptive,
              'descriptive_scope': 'Fixed-window raw descriptions, without stationarity claim or independent-sample CI. Formal conclusions come from the frozen analysis.'}
    (AUDIT / 'raw_validation_final.json').write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'status': result['status'], 'checks': checks}, ensure_ascii=False))
    if not all(checks.values()):
        raise SystemExit(1)


if __name__ == '__main__':
    main()
