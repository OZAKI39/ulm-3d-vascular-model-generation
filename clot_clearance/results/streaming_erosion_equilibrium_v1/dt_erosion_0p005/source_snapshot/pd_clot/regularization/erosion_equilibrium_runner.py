"""Equilibrium-gated slow exposure; no fictitious time in treatment clocks."""
import csv, gzip, hashlib, json, shutil, time, resource
from pathlib import Path
import numpy as np
from .equilibrium_adr import make_fixture, FrozenForceEvaluator, equilibrate
from .erosion import StreamingErosionDamage
from .erosion_runner import aggregate_bond_values
from ..damage import particle_damage, fragments

ROOT = Path(__file__).resolve().parents[2]


def write(path, value):
    with Path(path).open('x') as f: json.dump(value, f, indent=2, allow_nan=False)


def snapshot(out, c):
    write(out/'CONFIG.json', c)
    hashes = {}
    for p in sorted([*ROOT.glob('pd_clot/*.py'), *ROOT.glob('pd_clot/regularization/*.py')]):
        rel = p.relative_to(ROOT); dest = out/'source_snapshot'/rel
        dest.parent.mkdir(parents=True, exist_ok=True); shutil.copy2(p, dest)
        hashes[str(rel)] = hashlib.sha256(p.read_bytes()).hexdigest()
    write(out/'IDENTITY.json', dict(config_sha256=hashlib.sha256(json.dumps(c, sort_keys=True).encode()).hexdigest(),
                                  source_sha256=hashes, random_seed=None, randomness_used=False,
                                  time_convention='Exposure advances ONLY with accepted damage update; ADR has no physical time',
                                  physical_kinetic_energy='Not computed by quasi-static method'))


def save_equilibrium(out, cloud, g, c, result, exposure=0., phase='PRE_RELAXATION'):
    out.mkdir(parents=True, exist_ok=False)
    state = FrozenForceEvaluator(cloud, g, c)(result.x)
    np.savez_compressed(out/'state.npz', X=cloud.X, x=result.x, pairs=cloud.pairs, fixed=cloud.fixed,
                        volume=cloud.volume, integrity=g, bond_D=1-g,
                        epsilon_plus=state.epsilon_plus, internal_force_N=state.internal_force_N,
                        external_force_N=state.external_force_N, surface_traction_Pa=state.surface_traction_Pa,
                        fictitious_diagonal=result.fictitious_diagonal, pseudo_half_velocity=result.pseudo_half_velocity,
                        erosion_exposure_time_s=exposure)
    rows = [dict(h, phase=phase, erosion_exposure_time_s=exposure) for h in result.history]
    with gzip.open(out/'history.csv.gz', 'wt', newline='') as f:
        w = csv.DictWriter(f, list(rows[0])); w.writeheader(); w.writerows(rows)
    columns = {}
    for key in rows[0]:
        values = [h[key] for h in rows]
        if all(v is None for v in values): continue
        columns[key] = np.array(values) if isinstance(next(v for v in values if v is not None), str) else np.array([np.nan if v is None else v for v in values], float)
    np.savez_compressed(out/'history.npz', **columns)
    write(out/'checks.json', result.checks)
    write(out/'SUMMARY.json', dict(result.summary, phase=phase, erosion_exposure_time_s=exposure,
                                  final=dict(result.summary['final'], phase=phase, erosion_exposure_time_s=exposure)))
    return state


def run(c, output, solver=equilibrate):
    """An incomplete pre/update/post solve is saved and cannot advance damage."""
    if not c['damage']['enabled'] or c['damage']['mode'] != 'STREAMING_EROSION':
        raise ValueError('Enabled independent erosion mode required')
    if c['mechanical_equilibrium']['method'] != 'adaptive_dynamic_relaxation':
        raise ValueError('This runner requires the verified ADR interface')
    out = Path(output); out.mkdir(parents=True, exist_ok=False); snapshot(out, c)
    start = time.perf_counter(); cpu = time.process_time()
    cloud = make_fixture(c); x = cloud.X.copy(); zero = np.zeros_like(x)
    damage = StreamingErosionDamage(len(cloud.pairs), c['streaming_erosion'])
    de = c['streaming_erosion']['erosion_exposure_step_s']; count = c['streaming_erosion']['number_of_erosion_steps']
    if not np.isfinite(de) or de <= 0 or isinstance(count, bool) or int(count) != count or count < 1:
        raise ValueError('Positive slow step and integer count required')
    exposure = 0.; states = []; history = []; jumps = []; equilibrations = []; failure = None

    def record(force, label, converged):
        pd, _ = particle_damage(cloud, damage.g)
        _, _, topology = fragments(cloud, damage.g, 0.)
        rel = np.linalg.norm(force.residual[~cloud.fixed])/max(np.linalg.norm(force.external_force_N[~cloud.fixed]), c['mechanical_equilibrium']['force_floor_N'])
        row = dict(erosion_exposure_time_s=exposure, phase=label, equilibrium_passed=converged,
                   mean_damage=float(pd.mean()), maximum_damage=float(pd.max()), maximum_bond_damage=float(damage.D.max()),
                   maximum_epsilon_plus=float(force.epsilon_plus.max()), mean_epsilon_plus=float(force.epsilon_plus.mean()),
                   maximum_displacement_m=float(np.linalg.norm(x-cloud.X, axis=1).max()), stored_energy_J=force.stored_energy_J,
                   relative_force_residual=float(rel), broken_bond_count=int((~damage.active).sum()),
                   detached_volume_fraction=topology['detached_volume_fraction'])
        history.append(row)
        states.append(dict(x=x.copy(), bond_D=damage.D.copy(), integrity=damage.g.copy(), damage=pd,
                           epsilon_plus=force.epsilon_plus.copy(), particle_driver=aggregate_bond_values(cloud, force.epsilon_plus),
                           erosion_exposure_time_s=exposure, equilibrium_passed=converged))

    # n=0 starts from undeformed D=0; n=count settles final D without another exposure update.
    for n in range(int(count)+1):
        phase = 'PRE_RELAXATION' if n == 0 else 'EROSION_STEP'
        before_D = damage.D.copy(); before_g = damage.g.copy()
        result = solver(cloud, x, zero, damage.g, c)
        np.testing.assert_array_equal(damage.D, before_D); np.testing.assert_array_equal(damage.g, before_g)
        x = result.x.copy()
        force = save_equilibrium(out/f'equilibrium_{n:03d}', cloud, damage.g, c, result, exposure, phase)
        equilibrations.append(dict(index=n, exposure_s=exposure, phase=phase, converged=result.converged,
                                   iterations=result.summary['iterations'], relative_force_residual=result.summary['final']['relative_force_residual']))
        record(force, phase, result.converged)
        if not result.converged:
            failure = 'INCOMPLETE_MECHANICAL_EQUILIBRIUM'; break
        if n == count: break
        U_before = force.stored_energy_J; before_x = x.copy(); driver = force.epsilon_plus.copy()
        failed, old, increment = damage.advance(driver, de)
        exposure = (n+1)*de
        # The frozen mechanics supports the intended attached development case.
        # Topology change is retained, then stopped; no detached pseudo-dynamics.
        from ..fragment_mechanics import prepare_supported_shape
        from .fast_mechanics import evaluate_supported
        from .energy import potential
        shape = prepare_supported_shape(cloud, damage.g, c['safety'])
        after = evaluate_supported(cloud, x, damage.g, c['material'], c['safety'], shape)
        U_after = sum(potential(cloud, after))
        np.testing.assert_array_equal(x, before_x)
        jump = dict(update=n+1, erosion_exposure_start_s=n*de, erosion_exposure_end_s=exposure,
                    equilibrium_passed_before_update=True, equilibrium_relative_force_residual=result.summary['final']['relative_force_residual'],
                    U_before_damage_J=U_before, U_after_damage_J=U_after, delta_U_damage_J=U_before-U_after,
                    fixed_coordinates_at_damage_jump=True, coordinates_sha256=hashlib.sha256(x.tobytes()).hexdigest(),
                    maximum_damage_increment=float(increment.max()), broken_bond_count=int(len(failed)))
        jumps.append(jump)
        jump_dir = out/f'damage_jump_{n+1:03d}'; jump_dir.mkdir()
        np.savez_compressed(jump_dir/'state.npz', x=x, D_before=old, D_after=damage.D,
                            epsilon_plus_used=driver, rate_s_inv=damage.last_rate_s_inv)
        write(jump_dir/'AUDIT.json', jump)
        print(json.dumps(dict(update=n+1, exposure_s=exposure, iterations=result.summary['iterations'],
                              relative_force_residual=result.summary['final']['relative_force_residual'], delta_U_damage_J=U_before-U_after)), flush=True)
        if len(failed):
            failure = 'INCOMPLETE_TOPOLOGY_OUTSIDE_EQUILIBRIUM_SCOPE'
            np.savez_compressed(out/'stopped_state.npz', x=x, bond_D=damage.D, integrity=damage.g, erosion_exposure_time_s=exposure)
            break
    energy_negative = any(j['delta_U_damage_J'] < 0 for j in jumps)
    summary = dict(status=failure or 'COMPLETED_EQUILIBRATED_EROSION_DEVELOPMENT_RUN',
                   final=history[-1], accepted_damage_updates=len(jumps), requested_damage_updates=count,
                   saved_exposure_s=exposure, pre_relaxation=equilibrations[0], equilibrations=equilibrations,
                   signed_damage_energy_drop_J=sum(j['delta_U_damage_J'] for j in jumps),
                   damage_energy_increase_count=sum(j['delta_U_damage_J'] < 0 for j in jumps),
                   damage_energy_status='ENERGY_DISSIPATION_NOT_ESTABLISHED' if energy_negative else 'NO_NEGATIVE_DAMAGE_JUMP_IN_THIS_RUN',
                   physical_calibration='NOT PERFORMED', mesh_objectivity='STILL OPEN', wall_attachment='STILL OPEN',
                   real_microbubble_CFD='NOT USED', pseudo_time_used_as_exposure=False, represented_cycles_used=False,
                   mechanical_energy_budget='Only fixed-coordinate signed damage jumps; ADR pseudo-work is not physical work',
                   wall_time_s=time.perf_counter()-start, cpu_time_s=time.process_time()-cpu,
                   peak_process_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    np.savez_compressed(out/'states.npz', X=cloud.X, pairs=cloud.pairs, fixed=cloud.fixed, volume=cloud.volume,
                        **{k:np.array([s[k] for s in states]) for k in states[0]})
    write(out/'history.json', history); write(out/'damage_jumps.json', jumps); write(out/'SUMMARY.json', summary)
    return summary
