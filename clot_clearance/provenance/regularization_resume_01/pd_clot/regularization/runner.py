"""Extension of fragmentation.run with measured energy and consistent fluid data.

The legacy runner remains frozen. Integration, surface extraction, supported
NOSB mechanics and graph lineage are reused without prescribed separation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil
import time
import traceback
import warnings
import numpy as np
from ..geometry import make_cloud
from ..mechanics import timestep_limit
from ..damage import particle_damage
from ..fragment_mechanics import prepare_supported_shape, evaluate_supported
from ..fragment_topology import FragmentTracker, surface_particles
from ..fragment_fluid import relaxation_half_step
from ..output import write_json
from .damage import EnergyBondDamage
from .energy import EnergyLedger, kinetic, potential
from .flow import ConsistentFragmentFluid
from .diagnostics import classify, localization, ClassifiedClearance


def run(c, output, calibration, quiet=False, provider=None):
    out = Path(output)
    out.mkdir(parents=True, exist_ok=False)
    write_json(out/'CONFIG.json', c)
    cal = json.loads(Path(calibration).read_text()) if isinstance(calibration, (str, Path)) else calibration
    write_json(out/'CALIBRATION.json', cal)
    cloud = make_cloud(c['clot'])
    x, v = cloud.X.copy(), np.zeros_like(cloud.X)
    dt = c['simulation']['dt_s']
    limit = timestep_limit(cloud, c)
    if not 0 < dt <= limit:
        raise ValueError(f'Unsafe explicit timestep {dt}; estimate {limit}')
    dc, sim = c['damage'], c['simulation']
    period = 1/sim['representative_frequency_Hz']
    steps_per_cycle = int(round(period/dt))
    if abs(steps_per_cycle*dt-period) > 1e-12:
        raise ValueError('dt must divide the representative cycle')
    damage = EnergyBondDamage(cloud, c, cal)
    fluid = ConsistentFragmentFluid(c, provider=provider)
    tracker = FragmentTracker(cloud, c['transport']['x_clearance_m'])
    clearance = ClassifiedClearance(len(x))
    mass = cloud.volume*c['clot']['density_kg_m3']
    shape = prepare_supported_shape(cloud, damage.g, c['safety'])
    attached, components, stats = tracker.update(damage.g, x, v, 0)
    exposed, _ = surface_particles(cloud, damage.g, c['surface'])
    initial_exposed = exposed.copy()
    states, history, component_history = [], [], []
    N, t, tick = 0, 0., 0
    impulse = np.zeros_like(x)
    events = {name: None for name in ['first_bond_failure', 'first_detachment', 'first_resolved_detachment', 'first_clearance_crossing', 'first_reduced_rank']}
    root = Path(__file__).resolve().parents[2]
    source_files = [*root.glob('pd_clot/*.py'), *root.glob('pd_clot/regularization/*.py'), *root.glob('vendor/*.py')]
    provenance = out/'source_snapshot'
    for p in source_files:
        dst = provenance/p.relative_to(root)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, dst)
    write_json(out/'IDENTITY.json', dict(model='NOSB_PD_energy_regularization_development',
               no_prescribed_cleavage=True, initial_velocity='zero', forced_failure_verification=False,
               config_sha256=hashlib.sha256(json.dumps(c, sort_keys=True).encode()).hexdigest(),
               calibration_sha256=hashlib.sha256(json.dumps(cal, sort_keys=True).encode()).hexdigest(),
               source_sha256={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files},
               particle_count=len(x), bond_count=len(damage.g), dt_estimate_s=limit,
               load_and_transport_share_one_field=True, source_center_m=c['streaming']['bubble_center_m']))
    (out/'vtk').mkdir()
    start = time.perf_counter()
    minimum_J, maximum_acceleration = 1., 0.

    def forces():
        mech = evaluate_supported(cloud, x, damage.g, c['material'], c['safety'], shape)
        fs, normals, traction = fluid.surface_force(cloud, x, damage.g, exposed, attached, t)
        acc = (mech[0]+fs)/mass[:, None]
        acc[cloud.fixed] = 0.
        if not np.isfinite(acc).all():
            raise FloatingPointError('Nonfinite acceleration')
        return mech, acc, fs, normals, traction

    initial_mech = forces()[0]
    energy = EnergyLedger(sum(potential(cloud, initial_mech)))

    def apply_damping(decay):
        nonlocal v
        before = kinetic(mass, v)
        v *= decay
        energy.damping_dissipation_J += before-kinetic(mass, v)

    def apply_transport(free, half_dt):
        nonlocal v
        if not c['transport']['enabled']:
            return
        if c['transport']['mode'] != 'relaxation':
            from .transport import fragment_drag_half_step
            updated, drag_rows = fragment_drag_half_step(cloud, x, v, mass, components, category,
                                                         fluid.flow_provider, t, half_dt, c)
            transport_records.extend(dict(cycles=N, time_s=t, **r) for r in drag_rows)
        else:
            updated = relaxation_half_step(v, fluid.velocity(x, t), free, half_dt, c['transport']['tau_h_s'])
        energy.transport_work_J += kinetic(mass, updated)-kinetic(mass, v)
        impulse[:] += mass[:, None]*(updated-v)
        v = updated

    transport_records = []
    category, components, resolution = classify(cloud, damage.g, shape[3], attached, tracker.labels, components, c['fragment_resolution'])

    def record(step, recent, amplitude):
        import pyvista as pv
        mech, _, _, normals, traction = forces()
        D, Dmax = particle_damage(cloud, damage.g)
        elastic, stabilization = potential(cloud, mech)
        energy_row = energy.row(elastic, stabilization, kinetic(mass, v))
        if (abs(energy_row['relative_numerical_energy_residual']) > c['energy_audit']['maximum_relative_residual']
                and abs(energy_row['numerical_energy_residual_J']) > c['energy_audit']['minimum_absolute_tolerance_J']):
            raise FloatingPointError(f'Unexplained energy residual: {energy_row}')
        clearance_stats = clearance.update(tracker.cleared, category, cloud.volume)
        row = dict(macro_step=step, represented_cycles=N, mechanical_time_s=t,
                   mean_damage=float(D.mean()), maximum_damage=float(D.max()),
                   active_bond_count=int(damage.active.sum()), broken_bond_count=int((~damage.active).sum()),
                   broken_bond_fraction=float(np.mean(~damage.active)), particle_count=len(x), total_volume_m3=float(cloud.volume.sum()),
                   maximum_displacement_m=float(np.linalg.norm(x-cloud.X, axis=1).max()),
                   rank_counts=[int(np.sum(shape[3] == r)) for r in range(4)],
                   surface_particle_count=int(exposed.sum()), new_surface_particle_count=int(np.sum(exposed & ~initial_exposed)),
                   maximum_cycle_stretch_amplitude=float(amplitude.max()),
                   maximum_cycle_energy_driver_J=float(damage.last_driver_J.max()),
                   base_max_error_m=float(np.abs(x[cloud.fixed]-cloud.X[cloud.fixed]).max()),
                   **stats, **resolution, **clearance_stats, **energy_row,
                   **localization(cloud, D, damage.g, attached, np.array(c['streaming']['bubble_center_m']), c['localization']['influence_radius_m']))
        history.append(row)
        component_history.append(dict(cycles=N, mechanical_time_s=t, components=components))
        state = dict(x=x.copy(), v=v.copy(), integrity=damage.g.copy(), bond_D=damage.D.copy(), bond_active=damage.active.copy(),
                     damage=D.copy(), fragment_id=tracker.labels.copy(), attached=attached.copy(), surface=exposed.copy(),
                     rank=shape[3].copy(), cleared=tracker.cleared.copy(), hydro_impulse=impulse.copy(), cycles=N,
                     recent_broken=recent.copy(), amplitude=amplitude.copy(), resolution_class=category.copy(),
                     bond_cycle_energy_J=damage.last_driver_J.copy(), bond_dissipated_energy_J=damage.dissipated_J.copy())
        states.append(state)
        poly = pv.PolyData(x)
        arrays = dict(reference_position_m=cloud.X, displacement=x-cloud.X, velocity=v, damage=D, max_bond_damage=Dmax,
                      fragment_id=tracker.labels, attached_to_base=attached.astype(np.uint8), is_surface_particle=exposed.astype(np.uint8),
                      FIXED_BASE=cloud.fixed.astype(np.uint8), cleared_material=tracker.cleared.astype(np.uint8),
                      deformation_rank=shape[3], J_or_intrinsic_stretch=np.nan_to_num(mech[2], nan=0.),
                      strain_energy_density=mech[3], stabilization_energy_density=mech[4], particle_volume_m3=cloud.volume,
                      fluid_velocity_m_s=fluid.velocity(x, t), surface_normal=normals, surface_traction_Pa=traction,
                      resolution_class=category)
        for key, val in arrays.items():
            poly[key] = val
        poly.save(out/f'vtk/particles_{step:04d}.vtp')
        bonds = pv.PolyData(x, lines=np.column_stack((np.full(len(cloud.pairs), 2), cloud.pairs)).ravel())
        for key, val in dict(integrity=damage.g, accumulated_D=damage.D, active=damage.active.astype(np.uint8),
                             recently_broken=recent.astype(np.uint8), cycle_energy_J=damage.last_driver_J,
                             cumulative_damage_dissipation_J=damage.dissipated_J).items():
            bonds.cell_data[key] = val
        bonds.save(out/f'vtk/bonds_{step:04d}.vtp')
        write_json(out/'PROGRESS.json', row)
        if resolution['fragment_resolution_status'] == 'FRAGMENT_RESOLUTION_INADEQUATE':
            warnings.warn('FRAGMENT_RESOLUTION_INADEQUATE: detached singleton or low-rank fraction exceeds the configured quality flag', RuntimeWarning)
        if not quiet:
            print(json.dumps({k: row[k] for k in ['represented_cycles', 'mean_damage', 'broken_bond_fraction', 'detached_volume_fraction',
                                                 'resolved_fragment_volume_fraction', 'P_singleton', 'P_lowrank',
                                                 'relative_numerical_energy_residual', 'fragment_resolution_status']}), flush=True)

    try:
        record(0, np.zeros(len(damage.g), bool), np.zeros(len(damage.g)))
        with (out/'failure_ledger.jsonl').open('w') as legacy_ledger, (out/'fracture_energy_ledger.jsonl').open('w') as fracture_ledger:
            for macro in range(1, sim['number_of_macro_steps']+1):
                lo, hi = np.full(len(damage.g), np.inf), np.full(len(damage.g), -np.inf)
                impulse[:] = 0
                warmup = sim.get('warmup_cycles', 1) if macro == 1 else 0
                nsteps = steps_per_cycle*(sim['representative_cycles']+warmup)
                mech, acc, fs, _, _ = forces()
                free = (~attached) & (~cloud.fixed)
                decay = np.exp(-np.where(attached, sim['damping_s_inv'], sim['detached_damping_s_inv'])[:, None]*dt/2)
                for k in range(nsteps):
                    apply_damping(decay)
                    apply_transport(free, dt/2)
                    v += .5*dt*acc
                    dx = dt*v
                    x += dx
                    x[cloud.fixed] = cloud.X[cloud.fixed]
                    v[cloud.fixed] = 0
                    tick += 1
                    t = tick*dt
                    if np.linalg.norm(x-cloud.X, axis=1).max() > c['safety']['maximum_displacement_m']:
                        raise FloatingPointError('Configured displacement guard exceeded')
                    if np.linalg.norm(x[:, 1:3], axis=1).max() > c['pipe']['radius_m']*c['safety']['maximum_radius_ratio']:
                        raise FloatingPointError('Particle left pipe cross-section; no wall-contact model')
                    old_fs = fs
                    mech, acc, fs, _, _ = forces()
                    energy.external_work_J += float(np.sum(.5*(old_fs+fs)*dx))
                    peak = float(np.linalg.norm(acc, axis=1).max())
                    maximum_acceleration = max(maximum_acceleration, peak)
                    if peak > c['safety']['maximum_acceleration_m_s2']:
                        raise FloatingPointError('Acceleration guard exceeded')
                    if np.any(shape[3] > 0):
                        minimum_J = min(minimum_J, float(mech[2][shape[3] > 0].min()))
                    v += .5*dt*acc
                    apply_transport(free, dt/2)
                    apply_damping(decay)
                    v[cloud.fixed] = 0
                    if k >= warmup*steps_per_cycle:
                        stretch = np.linalg.norm(x[cloud.pairs[:, 1]]-x[cloud.pairs[:, 0]], axis=1)/cloud.length
                        lo, hi = np.minimum(lo, stretch), np.maximum(hi, stretch)
                amplitude = .5*(hi-lo)
                before_energy = sum(potential(cloud, mech))
                before_rank = shape[3].copy()
                ids, previous_D, increment = damage.advance(amplitude, dc['DeltaN'])
                N += dc['DeltaN']
                shape = prepare_supported_shape(cloud, damage.g, c['safety'])
                new_mech = evaluate_supported(cloud, x, damage.g, c['material'], c['safety'], shape)
                damage.dissipated_J += energy.damage_drop(before_energy, sum(potential(cloud, new_mech)),
                                                         (damage.D-previous_D)*damage.last_driver_J)
                recent = np.zeros(len(damage.g), bool)
                recent[ids] = True
                for bid in ids:
                    pair = cloud.pairs[bid]
                    old_row = dict(cycles=N, bond_id=int(bid), particle_ids=pair.tolist(),
                                   Q=float(amplitude[bid]), D_before=float(previous_D[bid]), D_after=float(damage.D[bid]),
                                   increment=float(increment[bid]), D_break=1.)
                    legacy_ledger.write(json.dumps(old_row)+'\n')
                    fracture_ledger.write(json.dumps(dict(failure_cycle=N, bond_id=int(bid), particle_ids=pair.tolist(),
                        pre_failure_damage=float(previous_D[bid]), local_cyclic_energy_driver_J=float(damage.last_driver_J[bid]),
                        accumulated_dissipated_energy_J=float(damage.dissipated_J[bid]), particle_spacing_m=cloud.spacing,
                        horizon_m=cloud.horizon, local_support_rank_before=before_rank[pair].tolist(),
                        local_support_rank_after=shape[3][pair].tolist(), calibrated_critical_energy_J=float(damage.critical_J[bid])))+'\n')
                legacy_ledger.flush()
                fracture_ledger.flush()
                attached, components, stats = tracker.update(damage.g, x, v, N)
                category, components, resolution = classify(cloud, damage.g, shape[3], attached, tracker.labels, components, c['fragment_resolution'])
                exposed, _ = surface_particles(cloud, damage.g, c['surface'])
                conditions = dict(first_bond_failure=bool(np.any(~damage.active)), first_detachment=bool(np.any(~attached)),
                                  first_resolved_detachment=bool(np.any(category == 1)),
                                  first_clearance_crossing=bool(tracker.cleared.any()), first_reduced_rank=bool(np.any(shape[3] < 3)))
                for name, condition in conditions.items():
                    if condition and events[name] is None:
                        events[name] = dict(cycles=N, macro_step=macro, mechanical_time_s=t)
                record(macro, recent, amplitude)
        summary = dict(status='COMPLETED_REGULARIZATION_DEVELOPMENT_RUN', label=c['verification_label'],
                       elapsed_wall_s=time.perf_counter()-start, events=events, final=history[-1],
                       minimum_J_or_intrinsic_stretch=minimum_J, maximum_acceleration_m_s2=maximum_acceleration,
                       all_particles_retained=True, forced_failure_verification=False, D_break=1.,
                       physical_predictivity_claim=False, source_center_m=c['streaming']['bubble_center_m'],
                       load_and_transport_share_one_field=True, fracture_resolution_status=resolution['fragment_resolution_status'],
                       energy_time_convention='Actual representative mechanical cycles only; not fictitious work multiplied by DeltaN',
                       calibration_scope=cal['calibration_scope'])
        write_json(out/'SUMMARY.json', summary)
        write_json(out/'history.json', history)
        write_json(out/'components.json', component_history)
        write_json(out/'lineage.json', tracker.lineage)
        write_json(out/'clearance_ledger.json', tracker.crossings)
        write_json(out/'transport_diagnostics.json', transport_records)
        np.savez_compressed(out/'states.npz', X=cloud.X, pairs=cloud.pairs, volume=cloud.volume, fixed=cloud.fixed,
                            **{key: np.array([state[key] for state in states]) for key in states[0]})
        with (out/'history.csv').open('w', newline='') as f:
            keys = [k for k, val in history[0].items() if not isinstance(val, list)]
            writer = csv.DictWriter(f, keys, extrasaction='ignore')
            writer.writeheader()
            writer.writerows(history)
        for kind in ['particles', 'bonds']:
            entries = '\n'.join(f'<DataSet timestep="{row["represented_cycles"]}" file="vtk/{kind}_{i:04d}.vtp"/>' for i, row in enumerate(history))
            (out/f'{kind}.pvd').write_text('<?xml version="1.0"?><VTKFile type="Collection" version="0.1"><Collection>'+entries+'</Collection></VTKFile>')
        return summary
    except Exception as error:
        write_json(out/'FAILED.json', dict(error=repr(error), cycles=N, mechanical_time_s=t, events=events, traceback=traceback.format_exc()))
        np.savez_compressed(out/'failed_state.npz', X=cloud.X, x=x, v=v, pairs=cloud.pairs, integrity=damage.g, bond_D=damage.D)
        write_json(out/'history_partial.json', history)
        raise


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--calibration', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    run(json.loads(args.config.read_text()), args.output, args.calibration)


if __name__ == '__main__':
    main()
