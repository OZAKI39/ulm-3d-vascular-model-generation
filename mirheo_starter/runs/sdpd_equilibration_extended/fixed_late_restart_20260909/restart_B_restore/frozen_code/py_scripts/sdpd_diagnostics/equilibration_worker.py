"""Thin observation extension of the existing native continuous SDPD worker.

Importing this module never initializes CUDA. Only the authorized MPI launcher
calls main. Diagnostic snapshots are not restart/checkpoint files.
"""
import csv
import json
import time
from pathlib import Path
import numpy as np
from py_scripts.fluid_physics.analysis import native_pressure
from py_scripts.fluid_physics.common import sha256_file, write_json
from .observations import velocity_statistics

RAW_COLUMNS = [
    'step', 'time_star', 'physical_time_s', 'density_time_star', 'pressure_force_time_star',
    'N', 'mass_star', 'total_mass_star', 'global_rho_star', 'global_rho_kg_m3',
    'kBT_raw_star', 'kBT_COM_star', 'temperature_K', 'component_kBT_x', 'component_kBT_y',
    'component_kBT_z', 'mean_vx', 'mean_vy', 'mean_vz', 'momentum_x', 'momentum_y', 'momentum_z',
    'pressure_star', 'pressure_pa', 'virial_pressure_star', 'kinetic_pressure_star',
    'kernel_number_mean', 'kernel_number_variance', 'kernel_number_min', 'kernel_number_max',
    'kernel_q10', 'kernel_q50', 'kernel_q90', 'kernel_cv', 'max_speed_star',
    'compute_wall_s', 'sampling_wall_s', 'chunk_compute_wall_s']
PRESSURE_DEFINITION = ('sum(stresses.xx + stresses.yy + stresses.zz)/(3V) + '
    'm*sum(|v-COM|^2)/(3V); compression positive, no second EOS term. '
    'Persistent saved_stresses from beforeIntegration; native Stress layout [xx,xy,xz,yy,yz,zz]. Force/density are pre-kick '
    '(step-1)*dt; momenta are post-kick step*dt. This is the existing discrete '
    'virial definition, with the phase difference explicitly retained.')


def mechanical_pressure(stress, vel, mass, volume, dt, step):
    stress, vel = np.asarray(stress, float), np.asarray(vel, float)
    if stress.shape != (len(vel), 6) or len(vel) <= 1 or not np.isfinite(stress).all():
        raise ValueError('INVALID_NATIVE_STRESS_LAYOUT_OR_VALUES')
    n = len(vel)
    raw = mass * np.sum(vel**2) / (3*n)
    mean = vel.mean(axis=0)
    stats = np.array([(step*dt, n, raw, *(mass*mean))],
                     dtype=[(k, float) for k in ['time','num_particles','kBT','vx','vy','vz']])
    virial_sum = float(stress[:, [0, 3, 5]].sum()/3)
    virial = np.array([((step-1)*dt, virial_sum)], dtype=[('time',float),('pressure',float)])
    p, _ = native_pressure(stats, virial, mass, volume, dt)
    return float(p[0]), virial_sum/volume


class EquilibrationObserver:
    def __init__(self, spec, directory, pv, channel):
        if spec['task']['force_star'] != 0 or spec['candidate']['method'] != 'SDPD':
            raise ValueError('UNFORCED_SDPD_ONLY')
        self.spec, self.d, self.pv, self.channel = spec, Path(directory), pv, channel
        self.volume = float(np.prod(spec['domain_star']))
        self.initial_ids = np.asarray(pv.get_indices(), np.int64)
        self.n = len(self.initial_ids)
        if self.n != spec['expected_N'] or len(np.unique(self.initial_ids)) != self.n:
            raise ValueError('INITIAL_PARTICLE_COUNT_OR_ID_ERROR')
        np.save(self.d/'initial_ids.npy', self.initial_ids, allow_pickle=False)
        identity = {name: sha256_file(self.d/name) for name in
                    ['initial_positions_global.npy','initial_velocities.npy','initial_ids.npy','initial_state.json']}
        write_json(self.d/'initial_state_identity.json', {
            'step':0, 'time_star':0., 'physical_time_s':0., 'sha256':identity,
            'continuity':'Fresh Uniform + one Gaussian draw and one COM subtraction; no restart.',
            'snapshot_is_checkpoint':False, 'native_random_seed':'Pinned rank + hash(pv); StepRandomGen retained by the same interaction.',
            'velocity_seed':spec['initial_velocity_seed']})
        self.stream = (self.d/'raw_statistics.csv').open('x', buffering=1)
        self.writer = csv.DictWriter(self.stream, fieldnames=RAW_COLUMNS)
        self.writer.writeheader()
        vel = np.load(self.d/'initial_velocities.npy', allow_pickle=False)
        row = self.velocity_row(0, vel)
        row.update(compute_wall_s=0., sampling_wall_s=0., chunk_compute_wall_s=0.)
        self.writer.writerow(row)  # t=0 pressure/kernel density genuinely unmeasured.
        self.snapshots = []

    def velocity_row(self, step, vel):
        s, mass = self.spec, self.spec['m_star']
        q = velocity_statistics(vel, mass)
        if len(vel) != self.n or not np.isfinite(vel).all():
            raise ValueError('PARTICLE_COUNT_OR_VELOCITY_ERROR')
        row = dict(step=step, time_star=step*s['dt_star'], physical_time_s=step*s['dt_star']*s['locked_units']['t0'],
                   N=self.n, mass_star=mass, total_mass_star=self.n*mass,
                   global_rho_star=self.n*mass/self.volume,
                   global_rho_kg_m3=self.n*mass/self.volume*s['locked_units']['si_per_star']['mass_density'],
                   kBT_raw_star=float(mass*np.sum(vel**2)/(3*self.n)), kBT_COM_star=q['COM_kBT'],
                   temperature_K=q['COM_kBT']/s['kBT_star']*s['target_temperature_K'], max_speed_star=q['max_speed'])
        for j, axis in enumerate('xyz'):
            row['component_kBT_'+axis] = q['component_kBT'][j]
            row['mean_v'+axis] = q['COM_velocity'][j]
            row['momentum_'+axis] = q['COM_velocity'][j]*self.n*mass
        return row

    def sample(self, step, vel, density, compute_s, sample_s, sample_started, chunk_s):
        row = self.velocity_row(step, vel)
        den = np.asarray(density, float)
        if len(den) != self.n or not np.isfinite(den).all() or np.min(den) <= 0:
            raise ValueError('INVALID_NATIVE_DENSITY_FIELD')
        p, virial = mechanical_pressure(self.channel(self.pv,'saved_stresses'), vel,
                                        self.spec['m_star'], self.volume, self.spec['dt_star'], step)
        q10, q50, q90 = np.quantile(den, [.1,.5,.9])
        row.update(density_time_star=(step-1)*self.spec['dt_star'],
                   pressure_force_time_star=(step-1)*self.spec['dt_star'], pressure_star=p,
                   pressure_pa=p*self.spec['locked_units']['si_per_star']['pressure'],
                   virial_pressure_star=virial, kinetic_pressure_star=p-virial,
                   kernel_number_mean=float(den.mean()), kernel_number_variance=float(den.var()),
                   kernel_number_min=float(den.min()), kernel_number_max=float(den.max()),
                   kernel_q10=float(q10), kernel_q50=float(q50), kernel_q90=float(q90),
                   kernel_cv=float(den.std()/den.mean()), compute_wall_s=compute_s,
                   sampling_wall_s=sample_s+time.monotonic()-sample_started, chunk_compute_wall_s=chunk_s)
        if not all(np.isfinite(v) for v in row.values()):
            raise ValueError('NONFINITE_OBSERVATION')
        self.writer.writerow(row)

    def snapshot(self, index, step):
        ids = np.asarray(self.pv.get_indices(), np.int64)
        if not np.array_equal(np.sort(ids), np.sort(self.initial_ids)):
            raise ValueError('PARTICLE_IDENTITIES_NOT_CONSERVED')
        path = self.d/f'snapshot_{index:02d}_ids.npy'
        np.save(path, ids, allow_pickle=False)
        stress_path=self.d/f'snapshot_{index:02d}_stresses.npy'
        np.save(stress_path,self.channel(self.pv,'saved_stresses'),allow_pickle=False)
        self.snapshots.append({'index':index,'step':step,'ids_sha256':sha256_file(path),
                               'stresses_sha256':sha256_file(stress_path),
                               'same_phase_positions_and_density':True,'snapshot_is_checkpoint':False})

    def finish(self, completion):
        self.stream.close()
        write_json(self.d/'equilibration_observer.json', {'pressure_definition':PRESSURE_DEFINITION,
            'snapshots':self.snapshots, 'initial_N':self.n, 'actual_steps':completion['steps'],
            'physical_time_s':completion['time_star']*self.spec['locked_units']['t0'],
            'continuous_coordinator_count':1, 'velocity_initialization_count':1,
            'additional_thermostat':False, 'velocity_rescaling_during_run':False,
            'sampling_interval_star':self.spec['snapshot_every']*self.spec['dt_star'],
            'difference_from_previous_observer':'Same 200-step cadence; one additional persistent pre-integration stress saver (device copy each step), six fixed snapshots with IDs and stress/velocity/density summaries. CPU snapshot pair audits run after GPU exit. Added cost is unmeasured and covered only by planning margins.',
            'stop_reason':completion['status'], 'restart_capability':'NOT_VERIFIED_NOT_USED'})


def main():
    from .gpu_worker import main as native_main
    native_main(observer_factory=EquilibrationObserver)


if __name__ == '__main__':
    main()
