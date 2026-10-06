"""Frozen intact normal/shear benchmarks; serial memory-bounded execution."""
import argparse
import csv
import gc
import hashlib
import json
import resource
import threading
import time
from pathlib import Path
import numpy as np
import psutil

from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.nonlocal_boundary import BlockBoundary, AffineSurrogate, block_family, convert_volume_force
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import FACE_NAMES
from pd_clot.pfpd_candidate.nonlocal_boundary.volume_force_conversion import stress_channel_load
from pd_clot.pfpd_candidate.nonlocal_boundary.volume_dirichlet_reference import constrained_support, impose
from pd_clot.pfpd_candidate.nonlocal_boundary.diagnostics import load_checks, support_profile, fields
from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_adapter import solve

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'verification/nonlocal_boundary_formulation_v1'
RESULTS = ROOT/'results/nonlocal_boundary_formulation_v1'


def plain(x):
    if isinstance(x, np.ndarray): return x.tolist()
    if isinstance(x, np.generic): return x.item()
    raise TypeError(type(x))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, default=plain, allow_nan=False)+'\n')


def load(path): return json.loads(Path(path).read_text())


def csvsave(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def intact_material(): return Material(model='svk', mu=1000., bulk=666.6666666666666)


def surrogate(kind):
    F = np.eye(3)
    if kind == 'normal': F[2, 2] += .003
    elif kind == 'shear': F[0, 2] = .003
    else: raise ValueError(kind)
    return AffineSurrogate(F, intact_material())


def integrity():
    hashes = load(OUT/'PROTECTED_INPUT_HASHES.json')
    different = [p for p, value in hashes.items() if hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != value]
    return dict(protected_files=len(hashes), differences=different, PASS=not different)


def guard_memory():
    mem = psutil.virtual_memory()
    if mem.available < .2*mem.total: raise MemoryError('20% system RAM headroom guard')


def benchmark_reference(kind, cfg, spec=None, ordering=False):
    spec = cfg['refinement'][0] if spec is None else spec
    guard_memory(); g, f = block_family(spec['n'], cfg['L_m'], spec['delta_m']); u = surrogate(kind)
    x = u.position(f.X); st = evaluate(f, x, u.material)
    converted = convert_volume_force(g, f, u, st)
    checks = load_checks(g, f.X, converted.nodal_force_N, u)
    profile, quantiles = support_profile(g, f, converted.nodal_force_N, spec['h_m'])
    metrics, plus = fields(g, f, st, u)
    reference = dict(kind=kind, geometry=spec, response=u.response(), face_targets=u.exact_face_integrals(g),
                     B1_conversion=checks, B1_support_quantiles=quantiles, affine_fields=metrics,
                     support_nodes=int(converted.support.sum()), direct_horizon_nodes=int(converted.seed.sum()),
                     out_of_support_residual_N=float(np.linalg.norm(converted.omitted_internal_force_N)),
                     construction_identity_max_N=float(np.max(abs(st.force+converted.nodal_force_N))),
                     identity_is_scientific_gate=False, volume_sum_ratio=float(f.volume.sum()/g.volume),
                     affine_energy_ratio=st.energy/(u.response()['W']*g.volume))
    P = u.response()['P']; channels = {}
    for name, Pc in [('normal', np.diag(np.diag(P))), ('tangential', P-np.diag(np.diag(P)))]:
        if np.linalg.norm(Pc) == 0:
            channels[name] = dict(status='EXACT_ZERO_TARGET_AND_CHANNEL')
            continue
        channel = stress_channel_load(f, st, Pc)
        channels[name] = load_checks(g, f.X, channel, u, P_override=Pc)
    reference['traction_channels'] = channels
    # Report a meaningful nonzero cut-load diagnostic, without pretending that
    # a volumetric load has a unique face association.
    halves = []
    for j, face in enumerate(FACE_NAMES):
        axis, sign = j//2, -1 if j % 2 == 0 else 1
        mask = sign*(f.X[:, axis]-g.center[axis]) > 0
        actual = converted.nodal_force_N[mask].sum(axis=0)
        target = u.exact_face_integrals(g)[j]['force']
        halves.append(dict(face=face, half_space_force_N=actual, face_target_N=target,
                           stable_scaled_difference=float(np.linalg.norm(actual-target)/checks['force_scale_N']),
                           interpretation='Diagnostic only: half-space force is not assigned surface traction'))
    reference['half_space_vs_face_diagnostic'] = halves
    arrays = dict(X=f.X, volume=f.volume, x=x, force_internal=st.force, load=converted.nodal_force_N,
                  support=converted.support, direct_region=converted.seed, Fbar=st.Fbar,
                  Psi_plus_nodal=plus, distance=g.distance(f.X))
    if ordering:
        # New independent label order; a new Family preserves physical geometry.
        rng = np.random.default_rng(9271); perm = rng.permutation(len(f.X)); inv = np.argsort(perm)
        pairs = np.sort(inv[f.pairs], axis=1); pairs = pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))]
        pf = Family.from_arrays(f.X[perm], f.volume[perm], pairs, f.horizon)
        ps = evaluate(pf, u.position(pf.X), u.material)
        error = float(np.linalg.norm(ps.force[inv]-st.force)/max(np.linalg.norm(st.force), 1e-30))
        reference['permutation'] = dict(seed=9271, relative_force_L2_error=error,
                                       PASS=error <= cfg['permutation_gate'])
        del pf, ps
    del st; gc.collect()
    return reference, arrays, profile


def field_profiles(g, f, x, state, external, fixed, u, nodal_plus, h):
    Pscale = np.linalg.norm(u.response()['P'])*max(g.lengths)**2
    residual = state.force+external; residual[fixed] = 0.
    depth = np.minimum(f.X[:, 2]-g.lower[2], g.upper[2]-f.X[:, 2]); layer = np.rint(depth/h).astype(int)
    rows = []
    for k in np.unique(layer):
        mask = layer == k
        active = mask & ~fixed
        rows.append(dict(depth_h=int(k), depth_m=float(k*h), depth_over_delta=float(k*h/f.horizon),
                         Fzz_mean=float(np.mean(state.Fbar[mask, 2, 2])), Fxz_mean=float(np.mean(state.Fbar[mask, 0, 2])),
                         Psi_plus_mean_Pa=float(nodal_plus[mask].mean()), Psi_plus_max_Pa=float(nodal_plus[mask].max()),
                         force_residual_RMS_scaled=float(np.linalg.norm(residual[mask])/np.sqrt(mask.sum())/Pscale),
                         free_nodes=int(active.sum()),
                         free_force_residual_RMS_scaled=float(np.linalg.norm(residual[active])/np.sqrt(active.sum())/Pscale) if active.any() else None))
    # Centerline first interior increment on both prescribed planes.
    jumps = []
    for sign in [-1, 1]:
        z0 = g.lower[2] if sign == -1 else g.upper[2]
        targets = [np.array([*g.center[:2], z0]), np.array([*g.center[:2], z0-sign*h])]
        ids = [int(np.argmin(np.linalg.norm(f.X-t, axis=1))) for t in targets]
        disp = x-f.X
        jumps.append(dict(side='bottom' if sign == -1 else 'top',
                          displacement_jump_m=disp[ids[1]]-disp[ids[0]],
                          surrogate_jump_m=u.displacement(f.X)[ids[1]]-u.displacement(f.X)[ids[0]]))
    return rows, jumps


def solve_candidate(kind, candidate, cfg):
    path = OUT/f'{kind}_{candidate}.json'
    if path.exists(): return load(path)
    guard_memory(); spec = cfg['refinement'][0]
    g, f = block_family(spec['n'], cfg['L_m'], spec['delta_m']); u = surrogate(kind)
    prescribed = u.position(f.X); fixed = g.dirichlet_mask(f.X)
    if candidate == 'B0':
        external = u.surface_nodal_load(g, f.X, g.neumann_faces); initial = prescribed.copy()
    elif candidate == 'B1':
        with np.load(OUT/f'{kind.upper()}_AFFINE.npz') as saved: external = saved['load'].copy()
        # Nonzero smooth perturbation, zero on all faces. Recovery is separate
        # from the exact cancellation identity used to construct the load.
        q = (f.X-g.lower)/g.lengths
        perturb = .02*.003*cfg['L_m']*np.prod(np.sin(np.pi*q), axis=1)
        initial = prescribed+perturb[:, None]*np.array([1., -.4, .7])
    elif candidate == 'B2':
        fixed, prescribed = constrained_support(g, f, u)
        external = np.zeros_like(f.X); initial = prescribed.copy()
    else: raise ValueError(candidate)
    scale = np.linalg.norm(u.response()['P'])*cfg['L_m']**2
    print('START', kind, candidate, 'free_dofs', 3*int((~fixed).sum()), flush=True)
    x, info = solve(f, u.material, prescribed, fixed, initial, external, spec['h_m'], scale, cfg)
    state = evaluate(f, x, u.material); metrics, plus = fields(g, f, state, u)
    total_external = external.copy(); total_external[fixed] = -state.force[fixed]
    checks = load_checks(g, f.X, total_external, u)
    profiles, jumps = field_profiles(g, f, x, state, external, fixed, u, plus, spec['h_m'])
    reactions = -(state.force+external)*fixed[:, None]
    row = dict(kind=kind, candidate=candidate, solve=info, fields=metrics, load_checks=checks,
               prescribed_load_checks=load_checks(g, f.X, external, u) if candidate == 'B1' else None,
               constraint_nodes=int(fixed.sum()), constrained_displacement_error_m=float(abs(x[fixed]-prescribed[fixed]).max()),
               relative_displacement_error=float(np.linalg.norm(x-prescribed)/max(np.linalg.norm(prescribed-f.X), 1e-30)),
               first_interior_jumps=jumps, reaction_top_N=reactions[g.face_mask(f.X, 'z+')].sum(axis=0),
               reaction_bottom_N=reactions[g.face_mask(f.X, 'z-')].sum(axis=0),
               interior_force_residual_L2_scaled=float(np.linalg.norm((state.force+external)[~fixed])/scale),
               exterior_nodes=0, all_energy_is_physical_domain=True,
               label='LOCALITY-STRONG REFERENCE' if candidate == 'B2' else candidate)
    save(path, row)
    np.savez_compressed(OUT/f'{kind}_{candidate}_STATE.npz', X=f.X, x=x, volume=f.volume,
                        external=external, total_external=total_external, reaction=reactions, fixed=fixed,
                        force_internal=state.force, Fbar=state.Fbar, Psi_plus_nodal=plus, initial=initial)
    csvsave(OUT/f'{kind}_{candidate}_FIELD_PROFILE.csv', profiles)
    print('DONE', kind, candidate, 'Fcore', metrics['F_error_core_RMS'], 'VW', checks['virtual_work_error'], flush=True)
    del state, f; gc.collect(); guard_memory()
    return row


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--stage', choices=['references', 'solves', 'refinement', 'all'], default='all')
    args = parser.parse_args(); cfg = load(OUT/'PREDECLARED_CONFIG.json')
    assert integrity()['PASS']; started = time.monotonic(); stop = threading.Event()
    resources = dict(peak_RSS_bytes=0, minimum_available_RAM_bytes=psutil.virtual_memory().available,
                     CPU_workers=1, BLAS_threads=1, GPU='GPU_NOT_USED_BY_DESIGN')
    def monitor():
        p = psutil.Process()
        while not stop.is_set():
            resources['peak_RSS_bytes'] = max(resources['peak_RSS_bytes'], p.memory_info().rss)
            resources['minimum_available_RAM_bytes'] = min(resources['minimum_available_RAM_bytes'], psutil.virtual_memory().available)
            stop.wait(.05)
    t = threading.Thread(target=monitor, daemon=True); t.start()
    try:
        if args.stage in ['references', 'all']:
            for kind in ['normal', 'shear']:
                if kind == 'shear':
                    normal = load(OUT/'NORMAL_REFERENCE.json')['B1_conversion']
                    if not (normal['resultant_PASS'] and normal['virtual_work_PASS']):
                        save(OUT/'SHEAR_BENCHMARK.json', dict(status='NOT_RUN_NORMAL_GATE')); break
                ref, a, profile = benchmark_reference(kind, cfg, ordering=True)
                save(OUT/f'{kind.upper()}_REFERENCE.json', ref)
                np.savez_compressed(OUT/f'{kind.upper()}_AFFINE.npz', **a)
                csvsave(OUT/f'{kind.upper()}_BOUNDARY_FORCE_SUPPORT.csv', profile)
                print('REFERENCE', kind, ref['B1_conversion']['resultant_error'], ref['B1_conversion']['moment_error'], ref['B1_conversion']['virtual_work_error'], flush=True)
                del a; gc.collect()
        if args.stage in ['solves', 'all']:
            for kind in ['normal', 'shear']:
                if not (OUT/f'{kind.upper()}_REFERENCE.json').exists(): continue
                for candidate in ['B0', 'B1', 'B2']: solve_candidate(kind, candidate, cfg)
        if args.stage in ['refinement', 'all']:
            if all((OUT/f'{kind.upper()}_REFERENCE.json').exists() and
                   (OUT/f'{kind}_B1.json').exists() and load(OUT/f'{kind}_B1.json')['solve']['converged'] and
                   load(OUT/f'{kind}_B1.json')['fields']['F_error_max'] <= cfg['field_accuracy_gate_relative_to_imposed_gradient'] and
                   all(load(OUT/f'{kind}_B1.json')['load_checks'][k] for k in ['resultant_PASS', 'moment_PASS', 'virtual_work_PASS']) and
                   all(load(OUT/f'{kind.upper()}_REFERENCE.json')['B1_conversion'][k]
                       for k in ['resultant_PASS', 'moment_PASS', 'virtual_work_PASS']) for kind in ['normal', 'shear']):
                rows = []
                for spec in cfg['refinement'][1:]:
                    for kind in ['normal', 'shear']:
                        ref, a, profile = benchmark_reference(kind, cfg, spec)
                        save(OUT/f"REFINEMENT_{spec['name']}_{kind}.json", ref)
                        csvsave(OUT/f"REFINEMENT_{spec['name']}_{kind}_SUPPORT.csv", profile)
                        rows.append(ref); del a; gc.collect()
                save(OUT/'REFINEMENT.json', rows)
                save(OUT/'REFINEMENT_STATUS.json', dict(status='COMPLETED_AFFINE_DIAGNOSTIC_ONLY'))
            else:
                save(OUT/'REFINEMENT_STATUS.json', dict(status='NOT_RUN_MAIN_B1_ACCEPTANCE_NOT_MET'))
    finally:
        resources['wall_seconds'] = time.monotonic()-started
        resources['peak_RSS_bytes'] = max(resources['peak_RSS_bytes'], resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
        resources['minimum_available_RAM_fraction'] = resources['minimum_available_RAM_bytes']/psutil.virtual_memory().total
        stop.set(); t.join(); save(OUT/f'RESOURCES_{args.stage}.json', resources)
        save(OUT/'HISTORICAL_INTEGRITY.json', integrity())


if __name__ == '__main__': main()
